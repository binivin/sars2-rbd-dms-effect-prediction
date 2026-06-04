from pathlib import Path
import numpy as np
import pandas as pd
import torch

from transformers import AutoTokenizer, AutoModel
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

PROCESSED_DIR = Path("data/processed")
OUT_DIR = Path("artifacts/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)

IN_PATH = PROCESSED_DIR / "rbd_dms_substitutions_esm_v0.csv"

DELTA_EMB_PATH = OUT_DIR / "rbd_delta_esm2_embeddings.csv"
DELTA_DATASET_PATH = PROCESSED_DIR / "rbd_dms_substitutions_delta_esm_v0.csv"

SITE_DETAIL_PATH = OUT_DIR / "delta_esm_repeated_site_heldout_detailed_metrics.csv"
SITE_SUMMARY_PATH = OUT_DIR / "delta_esm_repeated_site_heldout_summary_metrics.csv"

MUT_DETAIL_PATH = OUT_DIR / "delta_esm_mutation_heldout_detailed_metrics.csv"
MUT_SUMMARY_PATH = OUT_DIR / "delta_esm_mutation_heldout_summary_metrics.csv"

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
RANDOM_STATE = 42
N_REPEATS = 5
TEST_FRACTION = 0.2
BATCH_SIZE = 16

print("Reading ESM-merged dataset...")
df = pd.read_csv(IN_PATH)

print("\nInput shape:")
print(df.shape)

esm_ref_cols = [c for c in df.columns if c.startswith("esm2_")]
print("Reference ESM columns:", len(esm_ref_cols))

if len(esm_ref_cols) == 0:
    raise RuntimeError("No reference ESM columns found. Run run_esm_embedding_baseline.py first.")

df["mutation_id"] = (
    df["wildtype"].astype(str)
    + df["site"].astype(int).astype(str)
    + df["mutant"].astype(str)
)

df["background_mutation_id"] = (
    df["background"].astype(str)
    + "_"
    + df["mutation_id"].astype(str)
)

site_order = sorted(df["site"].dropna().astype(int).unique())
site_to_idx = {site: i for i, site in enumerate(site_order)}

global_site_wt = (
    df.groupby("site")["wildtype"]
    .agg(lambda x: x.value_counts().index[0])
    .to_dict()
)

background_sequences = {}

for bg in sorted(df["background"].dropna().astype(str).unique()):
    bg_df = df[df["background"].astype(str) == bg]

    bg_site_wt = (
        bg_df.groupby("site")["wildtype"]
        .agg(lambda x: x.value_counts().index[0])
        .to_dict()
    )

    seq_chars = []
    for site in site_order:
        aa = bg_site_wt.get(site, global_site_wt.get(site))
        if aa is None or len(str(aa)) != 1:
            raise RuntimeError(f"Could not infer amino acid for background={bg}, site={site}")
        seq_chars.append(str(aa))

    background_sequences[bg] = "".join(seq_chars)

print("\nBackground sequences reconstructed:")
for bg, seq in background_sequences.items():
    print(bg, len(seq), seq[:30])

if DELTA_EMB_PATH.exists():
    print("\nDelta ESM embedding file already exists. Loading:")
    print(DELTA_EMB_PATH)
    delta_emb = pd.read_csv(DELTA_EMB_PATH)

else:
    print("\nLoading ESM-2 model:")
    print(MODEL_NAME)

    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("Device:", device)

    model = model.to(device)
    model.eval()

    def embed_sequences(sequences, batch_size=BATCH_SIZE):
        all_embeddings = []

        for start in range(0, len(sequences), batch_size):
            batch = sequences[start:start + batch_size]

            inputs = tokenizer(
                batch,
                return_tensors="pt",
                padding=True,
                truncation=False,
            )
            inputs = {k: v.to(device) for k, v in inputs.items()}

            with torch.no_grad():
                outputs = model(**inputs)

            hidden = outputs.last_hidden_state.detach().cpu().numpy()
            all_embeddings.append(hidden)

            print(f"Embedded {min(start + batch_size, len(sequences))} / {len(sequences)}")

        return all_embeddings

    print("\nEmbedding background wild-type sequences...")
    bg_names = list(background_sequences.keys())
    bg_seqs = [background_sequences[bg] for bg in bg_names]

    bg_batches = embed_sequences(bg_seqs, batch_size=BATCH_SIZE)

    bg_embedding_by_name = {}
    cursor = 0

    for batch_hidden in bg_batches:
        for i in range(batch_hidden.shape[0]):
            bg = bg_names[cursor]
            seq_len = len(background_sequences[bg])

            bg_embedding_by_name[bg] = batch_hidden[i, 1:1 + seq_len, :]
            cursor += 1

    unique_rows = (
        df[
            [
                "background",
                "site",
                "wildtype",
                "mutant",
                "mutation_id",
                "background_mutation_id",
            ]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    print("\nUnique background-specific mutations:")
    print(len(unique_rows))

    mutant_sequences = []
    mutation_positions = []

    mismatch_count = 0

    for _, row in unique_rows.iterrows():
        bg = str(row["background"])
        site = int(row["site"])
        wt = str(row["wildtype"])
        mut = str(row["mutant"])

        seq = list(background_sequences[bg])
        pos = site_to_idx[site]

        if seq[pos] != wt:
            mismatch_count += 1

        seq[pos] = mut

        mutant_sequences.append("".join(seq))
        mutation_positions.append(pos)

    print("WT/background sequence mismatch count:", mismatch_count)

    print("\nEmbedding mutant sequences...")
    mutant_batches = embed_sequences(mutant_sequences, batch_size=BATCH_SIZE)

    delta_rows = []
    cursor = 0

    for batch_hidden in mutant_batches:
        for i in range(batch_hidden.shape[0]):
            row = unique_rows.iloc[cursor]
            bg = str(row["background"])
            pos = mutation_positions[cursor]

            mut_emb = batch_hidden[i, 1 + pos, :]
            wt_emb = bg_embedding_by_name[bg][pos, :]

            delta = mut_emb - wt_emb

            delta_l2 = float(np.linalg.norm(delta))
            wt_norm = float(np.linalg.norm(wt_emb))
            mut_norm = float(np.linalg.norm(mut_emb))

            cosine = float(
                np.dot(mut_emb, wt_emb) / ((mut_norm * wt_norm) + 1e-8)
            )

            out = {
                "background": row["background"],
                "site": int(row["site"]),
                "wildtype": row["wildtype"],
                "mutant": row["mutant"],
                "mutation_id": row["mutation_id"],
                "background_mutation_id": row["background_mutation_id"],
                "delta_esm2_l2": delta_l2,
                "mutant_wt_esm2_cosine": cosine,
            }

            for j, val in enumerate(delta):
                out[f"delta_esm2_{j:03d}"] = float(val)

            delta_rows.append(out)
            cursor += 1

    delta_emb = pd.DataFrame(delta_rows)
    delta_emb.to_csv(DELTA_EMB_PATH, index=False, encoding="utf-8-sig")

    print("\nSaved delta ESM embeddings:")
    print(DELTA_EMB_PATH)

merge_cols = [
    "background",
    "site",
    "wildtype",
    "mutant",
    "mutation_id",
    "background_mutation_id",
]

merged = df.merge(delta_emb, on=merge_cols, how="left")
merged.to_csv(DELTA_DATASET_PATH, index=False, encoding="utf-8-sig")

print("\nSaved delta ESM dataset:")
print(DELTA_DATASET_PATH)

print("Merged shape:")
print(merged.shape)

delta_cols = [
    c for c in merged.columns
    if c.startswith("delta_esm2_")
    and c.replace("delta_esm2_", "").isdigit()
]

print("Delta ESM vector columns:", len(delta_cols))

if len(delta_cols) == 0:
    raise RuntimeError("No delta ESM vector columns found.")

if merged[delta_cols].isna().any().any():
    raise RuntimeError("Some delta ESM vector features are missing after merge.")

if merged["delta_esm2_l2"].isna().any():
    raise RuntimeError("Some delta_esm2_l2 values are missing after merge.")

if merged["mutant_wt_esm2_cosine"].isna().any():
    raise RuntimeError("Some mutant_wt_esm2_cosine values are missing after merge.")

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
        "site",
        "site_from_rbd_start",
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
        numeric_cols += delta_cols
        numeric_cols += ["delta_esm2_l2", "mutant_wt_esm2_cosine"]

    cat_cols = ["wildtype", "mutant"]

    if "background" in feature_set:
        cat_cols += ["source", "background"]

    X_num = x[numeric_cols]
    X_cat = pd.get_dummies(x[cat_cols], drop_first=False)

    X = pd.concat([X_num, X_cat], axis=1)
    X = X.fillna(0)

    duplicate_cols = X.columns[X.columns.duplicated()].tolist()
    if duplicate_cols:
        raise RuntimeError(f"Duplicate feature columns found: {duplicate_cols[:20]}")

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


def fit_eval(df, train_idx, test_idx, target, feature_set, seed):
    X = build_features(df, feature_set)
    y = pd.to_numeric(df[target], errors="coerce")

    X_train = X.loc[train_idx]
    X_test = X.loc[test_idx]
    y_train = y.loc[train_idx]
    y_test = y.loc[test_idx]

    models = {
        "Ridge": make_pipeline(StandardScaler(), Ridge(alpha=1.0)),
        "RandomForest": RandomForestRegressor(
            n_estimators=40,
            max_depth=16,
            min_samples_leaf=2,
            random_state=seed,
            n_jobs=-1,
        ),
    }

    rows = []

    for model_name, model in models.items():
        print(f"{model_name} | {feature_set} | {target}")

        model.fit(X_train, y_train)
        pred = model.predict(X_test)

        rows.append({
            "feature_set": feature_set,
            "target": target,
            "model": model_name,
            **regression_metrics(y_test, pred),
            "n_train": len(train_idx),
            "n_test": len(test_idx),
        })

    return rows


feature_sets = [
    "mutation_only",
    "mutation_structure",
    "mutation_structure_refesm_background",
    "mutation_structure_deltaesm_background",
    "mutation_structure_refesm_deltaesm_background",
]

targets = [
    "ace2_binding_bg_z",
    "expression_bg_z",
]

all_site_rows = []

unique_sites = np.array(sorted(merged["site"].dropna().unique()))
n_test_sites = int(len(unique_sites) * TEST_FRACTION)

for repeat in range(1, N_REPEATS + 1):
    seed = RANDOM_STATE + repeat
    rng = np.random.default_rng(seed)

    shuffled = unique_sites.copy()
    rng.shuffle(shuffled)

    heldout_sites = set(shuffled[:n_test_sites])

    test_idx = merged[merged["site"].isin(heldout_sites)].index
    train_idx = merged[~merged["site"].isin(heldout_sites)].index

    print("\n" + "=" * 80)
    print(f"Repeated site-held-out {repeat}/{N_REPEATS}")
    print(f"Held-out sites: {len(heldout_sites)}")
    print(f"Train rows: {len(train_idx)}")
    print(f"Test rows: {len(test_idx)}")
    print("=" * 80)

    for fs in feature_sets:
        for target in targets:
            rows = fit_eval(merged, train_idx, test_idx, target, fs, seed)
            for r in rows:
                r["repeat"] = repeat
                r["seed"] = seed
                r["validation"] = "site_heldout"
            all_site_rows.extend(rows)

site_detail = pd.DataFrame(all_site_rows)
site_detail.to_csv(SITE_DETAIL_PATH, index=False, encoding="utf-8-sig")

site_summary = (
    site_detail
    .groupby(["feature_set", "target", "model"])
    .agg(
        R2_mean=("R2", "mean"),
        R2_std=("R2", "std"),
        MAE_mean=("MAE", "mean"),
        RMSE_mean=("RMSE", "mean"),
        Pearson_mean=("Pearson", "mean"),
        Spearman_mean=("Spearman", "mean"),
        n_repeats=("repeat", "count"),
    )
    .reset_index()
    .sort_values(["target", "model", "R2_mean"], ascending=[True, True, False])
)

site_summary.to_csv(SITE_SUMMARY_PATH, index=False, encoding="utf-8-sig")

all_mut_rows = []

unique_mutations = np.array(sorted(merged["mutation_id"].dropna().unique()))
n_test_mutations = int(len(unique_mutations) * TEST_FRACTION)

for repeat in range(1, N_REPEATS + 1):
    seed = RANDOM_STATE + repeat
    rng = np.random.default_rng(seed)

    shuffled = unique_mutations.copy()
    rng.shuffle(shuffled)

    heldout_mutations = set(shuffled[:n_test_mutations])

    test_idx = merged[merged["mutation_id"].isin(heldout_mutations)].index
    train_idx = merged[~merged["mutation_id"].isin(heldout_mutations)].index

    print("\n" + "=" * 80)
    print(f"Repeated mutation-held-out {repeat}/{N_REPEATS}")
    print(f"Held-out mutations: {len(heldout_mutations)}")
    print(f"Train rows: {len(train_idx)}")
    print(f"Test rows: {len(test_idx)}")
    print("=" * 80)

    for fs in feature_sets:
        for target in targets:
            rows = fit_eval(merged, train_idx, test_idx, target, fs, seed)
            for r in rows:
                r["repeat"] = repeat
                r["seed"] = seed
                r["validation"] = "mutation_heldout"
            all_mut_rows.extend(rows)

mut_detail = pd.DataFrame(all_mut_rows)
mut_detail.to_csv(MUT_DETAIL_PATH, index=False, encoding="utf-8-sig")

mut_summary = (
    mut_detail
    .groupby(["feature_set", "target", "model"])
    .agg(
        R2_mean=("R2", "mean"),
        R2_std=("R2", "std"),
        MAE_mean=("MAE", "mean"),
        RMSE_mean=("RMSE", "mean"),
        Pearson_mean=("Pearson", "mean"),
        Spearman_mean=("Spearman", "mean"),
        n_repeats=("repeat", "count"),
    )
    .reset_index()
    .sort_values(["target", "model", "R2_mean"], ascending=[True, True, False])
)

mut_summary.to_csv(MUT_SUMMARY_PATH, index=False, encoding="utf-8-sig")

print("\nSaved delta ESM site-held-out summary:")
print(SITE_SUMMARY_PATH)

print("\nDelta ESM repeated site-held-out summary:")
print(site_summary.to_string(index=False))

print("\nSaved delta ESM mutation-held-out summary:")
print(MUT_SUMMARY_PATH)

print("\nDelta ESM mutation-held-out summary:")
print(mut_summary.to_string(index=False))
