from pathlib import Path
import pandas as pd

PROCESSED_DIR = Path("data/processed")
OUT_DIR = Path("artifacts/tables")

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

in_path = PROCESSED_DIR / "rbd_dms_dataset_v0.csv"
out_path = PROCESSED_DIR / "rbd_dms_substitutions_v0.csv"
qc_path = OUT_DIR / "rbd_dms_substitutions_qc_summary.csv"

print("Reading dataset...")
df = pd.read_csv(in_path)

print("\nOriginal shape:")
print(df.shape)

sub = df[df["wildtype"] != df["mutant"]].copy()

required_cols = [
    "source",
    "background",
    "site",
    "wildtype",
    "mutant",
    "mutation",
    "ace2_binding",
    "expression",
]

sub = sub.dropna(subset=required_cols).copy()

sub["site"] = sub["site"].astype(int)
sub["wildtype"] = sub["wildtype"].astype(str)
sub["mutant"] = sub["mutant"].astype(str)
sub["mutation"] = sub["mutation"].astype(str)
sub["background"] = sub["background"].astype(str)
sub["source"] = sub["source"].astype(str)

sub["mutation_key"] = (
    sub["background"]
    + "_"
    + sub["wildtype"]
    + sub["site"].astype(str)
    + sub["mutant"]
)

sub.to_csv(out_path, index=False, encoding="utf-8-sig")

qc = (
    sub.groupby(["source", "background"])
    .agg(
        n_rows=("mutation", "count"),
        n_sites=("site", "nunique"),
        n_mutations=("mutation", "nunique"),
        ace2_mean=("ace2_binding", "mean"),
        ace2_min=("ace2_binding", "min"),
        ace2_max=("ace2_binding", "max"),
        expression_mean=("expression", "mean"),
        expression_min=("expression", "min"),
        expression_max=("expression", "max"),
    )
    .reset_index()
)

qc.to_csv(qc_path, index=False, encoding="utf-8-sig")

print("\nFiltered substitution dataset saved:")
print(out_path)

print("\nQC summary saved:")
print(qc_path)

print("\nFiltered shape:")
print(sub.shape)

print("\nBackground counts:")
print(sub["background"].value_counts())

print("\nSource counts:")
print(sub["source"].value_counts())

print("\nPreview:")
print(sub.head(20).to_string(index=False))
