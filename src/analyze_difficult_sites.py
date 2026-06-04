from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

PROCESSED_DIR = Path("data/processed")
OUT_DIR = Path("artifacts/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

IN_PATH = PROCESSED_DIR / "rbd_dms_substitutions_delta_esm_v0.csv"

PRED_PATH = OUT_DIR / "difficult_site_predictions.csv"
SITE_ERROR_PATH = OUT_DIR / "difficult_site_error_summary.csv"
TOP_SITE_PATH = OUT_DIR / "top_difficult_sites.csv"

RANDOM_STATE = 42
N_REPEATS = 5
TEST_FRACTION = 0.2
FEATURE_SET = "mutation_structure_deltaesm_background"

print("Reading delta ESM dataset...")
df = pd.read_csv(IN_PATH)
print("\nDataset shape:")
print(df.shape)

df["mutation_id"] = (
    df["wildtype"].astype(str)
    + df["site"].astype(int).astype(str)
    + df["mutant"].astype(str)
)

esm_ref_cols = [
    c for c in df.columns
    if c.startswith("esm2_") and c.replace("esm2_", "").isdigit()
]

delta_esm_cols = [
    c for c in df.columns
    if c.startswith("delta_esm2_") and c.replace("delta_esm2_", "").isdigit()
]

print("Reference ESM columns:", len(esm_ref_cols))
print("Delta ESM columns:", len(delta_esm_cols))

if len(delta_esm_cols) == 0:
    raise RuntimeError("No delta ESM columns found. Run run_delta_esm_validation.py first.")

hydropathy = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

volume = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
    "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7,
    "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
    "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}

charge = {
    "A": 0, "R": 1, "N": 0, "D": -1, "C": 0,
    "Q": 0, "E": -1, "G": 0, "H": 0.1, "I": 0,
    "L": 0, "K": 1, "M": 0, "F": 0, "P": 0,
    "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0,
}


def map_prop(series, prop_dict):
    return series.map(prop_dict).astype(float)


def build_features(data, feature_set):
    x = data.copy()

    x["site"] = pd.to_numeric(x["site"], errors="coerce")
    x["site_from_rbd_start"] = x["site"] - 331

    x["wt_hydropathy"] = map_prop(x["wildtype"], hydropathy)
    x["mut_hydropathy"] = map_prop(x["mutant"], hydropathy)
    x["delta_hydropathy"] = x["mut_hydropathy"] - x["wt_hydropathy"]

    x["wt_volume"] = map_prop(x["wildtype"], volume)
    x["mut_volume"] = map_prop(x["mutant"], volume)
    x["delta_volume"] = x["mut_volume"] - x["wt_volume"]

    x["wt_charge"] = map_prop(x["wildtype"], charge)
    x["mut_charge"] = map_prop(x["mutant"], charge)
    x["delta_charge"] = x["mut_charge"] - x["wt_charge"]

    x["is_to_proline"] = (x["mutant"] == "P").astype(int)
    x["is_from_proline"] = (x["wildtype"] == "P").astype(int)
    x["is_to_cysteine"] = (x["mutant"] == "C").astype(int)
    x["is_from_cysteine"] = (x["wildtype"] == "C").astype(int)

    numeric_cols = [
        "site", "site_from_rbd_start",
        "wt_hydropathy", "mut_hydropathy", "delta_hydropathy",
        "wt_volume", "mut_volume", "delta_volume",
        "wt_charge", "mut_charge", "delta_charge",
        "is_to_proline", "is_from_proline", "is_to_cysteine", "is_from_cysteine",
    ]

    structural_cols = [
        "is_rbm_rough",
        "min_dist_to_ace2_angstrom",
        "is_ace2_interface_5A",
        "n_ace2_contact_residues_5A",
        "structure_site_available",
    ]

    if "structure" in feature_set:
        numeric_cols += structural_cols

    if "refesm" in feature_set:
        numeric_cols += esm_ref_cols

    if "deltaesm" in feature_set:
        numeric_cols += delta_esm_cols
        numeric_cols += ["delta_esm2_l2", "mutant_wt_esm2_cosine"]

    cat_cols = ["wildtype", "mutant"]

    if "background" in feature_set:
        cat_cols += ["source", "background"]

    X_num = x[numeric_cols]
    X_cat = pd.get_dummies(x[cat_cols], drop_first=False)
    X = pd.concat([X_num, X_cat], axis=1).fillna(0)

    duplicate_cols = X.columns[X.columns.duplicated()].tolist()
    if duplicate_cols:
        raise RuntimeError(f"Duplicate feature columns found: {duplicate_cols[:20]}")

    return X


def regression_metrics(y_true, y_pred):
    return {
        "R2": r2_score(y_true, y_pred),
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
        "Pearson": pd.Series(y_true).corr(pd.Series(y_pred), method="pearson"),
        "Spearman": pd.Series(y_true).corr(pd.Series(y_pred), method="spearman"),
    }


print("\nBuilding feature matrix...")
X = build_features(df, FEATURE_SET)
print("Feature matrix shape:")
print(X.shape)

targets = ["ace2_binding_bg_z", "expression_bg_z"]
all_predictions = []
all_metrics = []

unique_sites = np.array(sorted(df["site"].dropna().unique()))
n_test_sites = int(len(unique_sites) * TEST_FRACTION)

for repeat in range(1, N_REPEATS + 1):
    seed = RANDOM_STATE + repeat
    rng = np.random.default_rng(seed)

    shuffled_sites = unique_sites.copy()
    rng.shuffle(shuffled_sites)
    heldout_sites = set(shuffled_sites[:n_test_sites])

    test_idx = df[df["site"].isin(heldout_sites)].index
    train_idx = df[~df["site"].isin(heldout_sites)].index

    print("\n" + "=" * 80)
    print(f"Repeat {repeat}/{N_REPEATS}")
    print(f"Held-out sites: {len(heldout_sites)}")
    print(f"Train rows: {len(train_idx)}")
    print(f"Test rows: {len(test_idx)}")
    print("=" * 80)

    for target in targets:
        print(f"\nTraining RandomForest | target={target}")
        y = pd.to_numeric(df[target], errors="coerce")

        X_train = X.loc[train_idx]
        X_test = X.loc[test_idx]
        y_train = y.loc[train_idx]
        y_test = y.loc[test_idx]

        model = RandomForestRegressor(
            n_estimators=80,
            max_depth=16,
            min_samples_leaf=2,
            random_state=seed,
            n_jobs=-1,
        )

        model.fit(X_train, y_train)
        pred = model.predict(X_test)

        metrics = regression_metrics(y_test, pred)
        metrics.update({
            "repeat": repeat,
            "seed": seed,
            "target": target,
            "feature_set": FEATURE_SET,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "n_test_sites": len(heldout_sites),
        })
        all_metrics.append(metrics)

        print("Metrics:")
        print(metrics)

        pred_df = df.loc[test_idx, [
            "source", "background", "site", "wildtype", "mutant", "mutation", "mutation_id",
            "ace2_binding_bg_z", "expression_bg_z",
            "is_rbm_rough", "is_ace2_interface_5A", "min_dist_to_ace2_angstrom",
            "n_ace2_contact_residues_5A", "structure_site_available",
        ]].copy()

        pred_df["repeat"] = repeat
        pred_df["seed"] = seed
        pred_df["target"] = target
        pred_df["y_true"] = y_test.values
        pred_df["y_pred"] = pred
        pred_df["error"] = pred_df["y_pred"] - pred_df["y_true"]
        pred_df["abs_error"] = pred_df["error"].abs()
        pred_df["squared_error"] = pred_df["error"] ** 2

        all_predictions.append(pred_df)

predictions = pd.concat(all_predictions, ignore_index=True)
predictions.to_csv(PRED_PATH, index=False, encoding="utf-8-sig")

metrics_df = pd.DataFrame(all_metrics)
print("\nRepeated site-held-out metrics:")
print(metrics_df.to_string(index=False))

site_error = (
    predictions
    .groupby(["target", "site"])
    .agg(
        n_rows=("abs_error", "count"),
        mean_abs_error=("abs_error", "mean"),
        median_abs_error=("abs_error", "median"),
        rmse=("squared_error", lambda x: np.sqrt(np.mean(x))),
        mean_error=("error", "mean"),
        y_true_mean=("y_true", "mean"),
        y_true_std=("y_true", "std"),
        y_pred_mean=("y_pred", "mean"),
        is_rbm_rough=("is_rbm_rough", "max"),
        is_ace2_interface_5A=("is_ace2_interface_5A", "max"),
        min_dist_to_ace2_angstrom=("min_dist_to_ace2_angstrom", "min"),
        n_ace2_contact_residues_5A=("n_ace2_contact_residues_5A", "max"),
        structure_site_available=("structure_site_available", "max"),
    )
    .reset_index()
)

site_error["difficulty_rank_in_target"] = (
    site_error.groupby("target")["mean_abs_error"].rank(method="dense", ascending=False)
)
site_error.to_csv(SITE_ERROR_PATH, index=False, encoding="utf-8-sig")

ace2 = site_error[site_error["target"] == "ace2_binding_bg_z"].copy()
expr = site_error[site_error["target"] == "expression_bg_z"].copy()

ace2 = ace2.rename(columns={
    "mean_abs_error": "ace2_mae",
    "rmse": "ace2_rmse",
    "mean_error": "ace2_bias",
    "difficulty_rank_in_target": "ace2_rank",
})

expr = expr.rename(columns={
    "mean_abs_error": "expression_mae",
    "rmse": "expression_rmse",
    "mean_error": "expression_bias",
    "difficulty_rank_in_target": "expression_rank",
})

keep_common = [
    "site", "is_rbm_rough", "is_ace2_interface_5A", "min_dist_to_ace2_angstrom",
    "n_ace2_contact_residues_5A", "structure_site_available",
]

combined = ace2[
    keep_common + ["ace2_mae", "ace2_rmse", "ace2_bias", "ace2_rank"]
].merge(
    expr[["site", "expression_mae", "expression_rmse", "expression_bias", "expression_rank"]],
    on="site",
    how="outer",
)

combined["combined_mae"] = combined[["ace2_mae", "expression_mae"]].mean(axis=1)
combined["combined_rank_score"] = combined[["ace2_rank", "expression_rank"]].mean(axis=1)
combined = combined.sort_values(["combined_mae", "ace2_mae", "expression_mae"], ascending=False)
combined.to_csv(TOP_SITE_PATH, index=False, encoding="utf-8-sig")

print("\nSaved prediction-level errors:")
print(PRED_PATH)
print("\nSaved site-level error summary:")
print(SITE_ERROR_PATH)
print("\nSaved top difficult sites:")
print(TOP_SITE_PATH)

print("\nTop 30 difficult sites by combined MAE:")
print(
    combined[[
        "site", "combined_mae", "ace2_mae", "expression_mae", "ace2_bias", "expression_bias",
        "is_rbm_rough", "is_ace2_interface_5A", "min_dist_to_ace2_angstrom",
        "n_ace2_contact_residues_5A",
    ]].head(30).to_string(index=False)
)

print("\nTop 20 difficult sites for ACE2 binding:")
print(
    site_error[site_error["target"] == "ace2_binding_bg_z"]
    .sort_values("mean_abs_error", ascending=False)[[
        "site", "mean_abs_error", "rmse", "mean_error", "is_rbm_rough",
        "is_ace2_interface_5A", "min_dist_to_ace2_angstrom", "n_ace2_contact_residues_5A",
    ]]
    .head(20)
    .to_string(index=False)
)

print("\nTop 20 difficult sites for expression:")
print(
    site_error[site_error["target"] == "expression_bg_z"]
    .sort_values("mean_abs_error", ascending=False)[[
        "site", "mean_abs_error", "rmse", "mean_error", "is_rbm_rough",
        "is_ace2_interface_5A", "min_dist_to_ace2_angstrom", "n_ace2_contact_residues_5A",
    ]]
    .head(20)
    .to_string(index=False)
)
