from pathlib import Path
import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

PROCESSED_DIR = Path("data/processed")
OUT_DIR = Path("artifacts/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

in_path = PROCESSED_DIR / "rbd_dms_substitutions_v0.csv"
normalized_path = PROCESSED_DIR / "rbd_dms_substitutions_normalized_v0.csv"
metrics_path = OUT_DIR / "rbd_baseline_normalized_model_metrics.csv"
norm_summary_path = OUT_DIR / "rbd_normalization_summary.csv"

RANDOM_STATE = 42

print("Reading substitution dataset...")
df = pd.read_csv(in_path)

print("\nOriginal shape:")
print(df.shape)

group_cols = ["source", "background"]

def group_zscore(s):
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std()
    if pd.isna(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean()) / sd

for col in ["ace2_binding", "expression"]:
    df[f"{col}_bg_z"] = df.groupby(group_cols)[col].transform(group_zscore)
    df[f"{col}_bg_pct"] = df.groupby(group_cols)[col].rank(pct=True)

df = df.dropna(subset=["ace2_binding_bg_z", "expression_bg_z"]).copy()

df.to_csv(normalized_path, index=False, encoding="utf-8-sig")

norm_summary = (
    df.groupby(group_cols)
    .agg(
        n_rows=("mutation", "count"),
        ace2_raw_mean=("ace2_binding", "mean"),
        ace2_raw_std=("ace2_binding", "std"),
        expression_raw_mean=("expression", "mean"),
        expression_raw_std=("expression", "std"),
        ace2_z_mean=("ace2_binding_bg_z", "mean"),
        ace2_z_std=("ace2_binding_bg_z", "std"),
        expression_z_mean=("expression_bg_z", "mean"),
        expression_z_std=("expression_bg_z", "std"),
    )
    .reset_index()
)

norm_summary.to_csv(norm_summary_path, index=False, encoding="utf-8-sig")

print("\nSaved normalized dataset:")
print(normalized_path)

print("\nSaved normalization summary:")
print(norm_summary_path)

print("\nNormalization summary:")
print(norm_summary.to_string(index=False))

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

def build_features(data, include_background=False):
    x = data.copy()

    x["site"] = pd.to_numeric(x["site"], errors="coerce")
    x["site_from_rbd_start"] = x["site"] - 331
    x["is_rbm_rough"] = ((x["site"] >= 438) & (x["site"] <= 506)).astype(int)

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
        "site",
        "site_from_rbd_start",
        "is_rbm_rough",
        "wt_hydropathy",
        "mut_hydropathy",
        "delta_hydropathy",
        "wt_volume",
        "mut_volume",
        "delta_volume",
        "wt_charge",
        "mut_charge",
        "delta_charge",
        "is_to_proline",
        "is_from_proline",
        "is_to_cysteine",
        "is_from_cysteine",
    ]

    cat_cols = ["wildtype", "mutant"]

    if include_background:
        cat_cols += ["source", "background"]

    X_num = x[numeric_cols]
    X_cat = pd.get_dummies(x[cat_cols], drop_first=False)

    X = pd.concat([X_num, X_cat], axis=1)
    X = X.fillna(0)

    return X

def regression_metrics(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    return {
        "R2": r2_score(y_true, y_pred),
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
        "Pearson": pd.Series(y_true).corr(pd.Series(y_pred), method="pearson"),
        "Spearman": pd.Series(y_true).corr(pd.Series(y_pred), method="spearman"),
    }

def run_split(df, split_name, train_idx, test_idx, target, feature_set_name, include_background):
    X = build_features(df, include_background=include_background)
    y = pd.to_numeric(df[target], errors="coerce")

    X_train = X.loc[train_idx]
    X_test = X.loc[test_idx]
    y_train = y.loc[train_idx]
    y_test = y.loc[test_idx]

    models = {
        "Ridge": make_pipeline(StandardScaler(), Ridge(alpha=1.0)),
        "RandomForest": RandomForestRegressor(
            n_estimators=120,
            max_depth=18,
            min_samples_leaf=2,
            random_state=RANDOM_STATE,
            n_jobs=-1,
        ),
    }

    rows = []

    for model_name, model in models.items():
        print(
            f"\nTraining {model_name} | "
            f"features={feature_set_name} | "
            f"split={split_name} | target={target}"
        )

        model.fit(X_train, y_train)
        pred = model.predict(X_test)

        m = regression_metrics(y_test, pred)

        rows.append({
            "feature_set": feature_set_name,
            "split": split_name,
            "target": target,
            "model": model_name,
            **m,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
        })

    return rows

rng = np.random.default_rng(RANDOM_STATE)

all_indices = np.array(df.index)
rng.shuffle(all_indices)

n_test = int(len(all_indices) * 0.2)
random_test_idx = all_indices[:n_test]
random_train_idx = all_indices[n_test:]

unique_sites = df["site"].dropna().unique()
rng.shuffle(unique_sites)

n_site_test = int(len(unique_sites) * 0.2)
heldout_sites = set(unique_sites[:n_site_test])

site_test_idx = df[df["site"].isin(heldout_sites)].index
site_train_idx = df[~df["site"].isin(heldout_sites)].index

print("\nRandom split:")
print(f"train={len(random_train_idx)}, test={len(random_test_idx)}")

print("\nSite-held-out split:")
print(f"held-out sites={len(heldout_sites)}")
print(f"train={len(site_train_idx)}, test={len(site_test_idx)}")

targets = [
    "ace2_binding_bg_z",
    "expression_bg_z",
]

feature_sets = {
    "mutation_only": False,
    "with_background": True,
}

all_metrics = []

for feature_set_name, include_background in feature_sets.items():
    for target in targets:
        all_metrics.extend(
            run_split(
                df,
                split_name="random_split",
                train_idx=random_train_idx,
                test_idx=random_test_idx,
                target=target,
                feature_set_name=feature_set_name,
                include_background=include_background,
            )
        )

        all_metrics.extend(
            run_split(
                df,
                split_name="site_held_out",
                train_idx=site_train_idx,
                test_idx=site_test_idx,
                target=target,
                feature_set_name=feature_set_name,
                include_background=include_background,
            )
        )

metrics = pd.DataFrame(all_metrics)
metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")

print("\nSaved normalized model metrics:")
print(metrics_path)

print("\nNormalized model metrics:")
print(metrics.to_string(index=False))
