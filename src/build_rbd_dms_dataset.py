from pathlib import Path
import pandas as pd
import numpy as np
import re
import urllib.request

RAW_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")
OUT_DIR = Path("artifacts/tables")

RAW_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

URLS = {
    "wuhan_2020_single_mut_effects": {
        "url": "https://media.githubusercontent.com/media/jbloomlab/SARS-CoV-2-RBD_DMS/master/results/single_mut_effects/single_mut_effects.csv",
        "filename": "wuhan_2020_single_mut_effects.csv",
        "default_background": "Wuhan-Hu-1",
    },
    "variant_backgrounds_2022_final_scores": {
        "url": "https://media.githubusercontent.com/media/jbloomlab/SARS-CoV-2-RBD_DMS_variants/main/results/final_variant_scores/final_variant_scores.csv",
        "filename": "variant_backgrounds_2022_final_variant_scores.csv",
        "default_background": None,
    },
}

out_dataset_path = PROCESSED_DIR / "rbd_dms_dataset_v0.csv"
out_qc_path = OUT_DIR / "rbd_dms_dataset_qc_summary.csv"


def download_file(url, path):
    if path.exists() and path.stat().st_size > 100:
        print(f"Already exists: {path}")
        return

    print(f"Downloading: {url}")
    urllib.request.urlretrieve(url, path)

    head = path.read_bytes()[:120].decode(errors="ignore")
    if "git-lfs.github.com" in head:
        raise RuntimeError(
            f"{path} looks like a Git LFS pointer, not the real CSV. "
            "Use the media.githubusercontent.com URL."
        )


def find_col(df, candidates, contains_any=None):
    cols = list(df.columns)
    lower_map = {c.lower(): c for c in cols}

    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    if contains_any:
        for col in cols:
            col_l = col.lower()
            if any(key.lower() in col_l for key in contains_any):
                return col

    return None


def zscore(s):
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std()
    if pd.isna(sd) or sd == 0:
        return s * np.nan
    return (s - s.mean()) / sd


def standardize_rbd_dms(df, source_name, default_background=None):
    df = df.copy()

    print(f"\n[{source_name}] columns:")
    print(list(df.columns))

    site_col = find_col(df, ["site", "site_SARS2", "site_RBD", "position", "pos"])
    wt_col = find_col(df, ["wildtype", "wt", "wildtype_aa", "aa"])
    mut_col = find_col(df, ["mutant", "mutation", "mutant_aa", "aa_sub"])

    bg_col = find_col(
        df,
        ["background", "target", "library", "variant", "strain"],
        contains_any=["background", "target", "library", "variant", "strain"],
    )

    binding_col = find_col(
        df,
        ["bind_avg", "binding", "ace2_binding", "bind", "delta_bind"],
        contains_any=["bind"],
    )

    expression_col = find_col(
        df,
        ["expr_avg", "expression", "expr", "delta_expr"],
        contains_any=["expr"],
    )

    required = {
        "site": site_col,
        "wildtype": wt_col,
        "mutant": mut_col,
        "ace2_binding": binding_col,
        "expression": expression_col,
    }

    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise ValueError(
            f"[{source_name}] Could not find required columns: {missing}\n"
            f"Available columns are:\n{list(df.columns)}"
        )

    out = pd.DataFrame(index=df.index)
    out["source"] = source_name

    if bg_col is not None:
        out["background"] = df[bg_col].astype(str)
    else:
        out["background"] = default_background

    out["site"] = pd.to_numeric(df[site_col], errors="coerce").astype("Int64")
    out["wildtype"] = df[wt_col].astype(str)
    out["mutant"] = df[mut_col].astype(str)

    out["mutation"] = (
        out["wildtype"].astype(str)
        + out["site"].astype(str)
        + out["mutant"].astype(str)
    )

    out["ace2_binding"] = pd.to_numeric(df[binding_col], errors="coerce")
    out["expression"] = pd.to_numeric(df[expression_col], errors="coerce")

    out["ace2_binding_z"] = zscore(out["ace2_binding"])
    out["expression_z"] = zscore(out["expression"])
    out["functional_score_temp"] = out[["ace2_binding_z", "expression_z"]].mean(axis=1)

    out = out[
        [
            "source",
            "background",
            "site",
            "wildtype",
            "mutant",
            "mutation",
            "ace2_binding",
            "expression",
            "ace2_binding_z",
            "expression_z",
            "functional_score_temp",
        ]
    ]

    return out


all_tables = []

for source_name, info in URLS.items():
    raw_path = RAW_DIR / info["filename"]
    download_file(info["url"], raw_path)

    print(f"\nReading: {raw_path}")
    df = pd.read_csv(raw_path)

    standardized = standardize_rbd_dms(
        df,
        source_name=source_name,
        default_background=info["default_background"],
    )

    all_tables.append(standardized)


rbd = pd.concat(all_tables, ignore_index=True)

rbd = rbd.dropna(subset=["site", "wildtype", "mutant"])
rbd = rbd[rbd["wildtype"].str.len() == 1]
rbd = rbd[rbd["mutant"].str.len() == 1]

rbd.to_csv(out_dataset_path, index=False, encoding="utf-8-sig")

qc = (
    rbd.groupby(["source", "background"], dropna=False)
    .agg(
        n_rows=("mutation", "count"),
        n_sites=("site", "nunique"),
        n_mutations=("mutation", "nunique"),
        ace2_missing=("ace2_binding", lambda x: x.isna().sum()),
        expression_missing=("expression", lambda x: x.isna().sum()),
        ace2_mean=("ace2_binding", "mean"),
        expression_mean=("expression", "mean"),
    )
    .reset_index()
)

qc.to_csv(out_qc_path, index=False, encoding="utf-8-sig")

print("\nSaved dataset:")
print(out_dataset_path)

print("\nSaved QC summary:")
print(out_qc_path)

print("\nDataset shape:")
print(rbd.shape)

print("\nBackground counts:")
print(rbd["background"].value_counts(dropna=False).head(20))

print("\nPreview:")
print(rbd.head(20).to_string(index=False))
