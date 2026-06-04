# %% [markdown]
# # RBD DMS Interpolation vs Extrapolation
#
# This notebook summarizes RBD DMS effect prediction results and compares interpolation-style validation with unseen-site extrapolation.
#
# **GitHub repository:** [sars2-rbd-dms-effect-prediction](https://github.com/binivin/sars2-rbd-dms-effect-prediction)
#
# **Kaggle dataset:** [RBD DMS Effect Prediction](https://www.kaggle.com/datasets/binivin/rbd-dms-effect-prediction)
#
# The full GitHub repository contains the reproducible analysis scripts, checkpoint summaries, and project documentation. This Kaggle notebook provides a concise visualization-centered summary of the final result tables.

# %%
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

INPUT_ROOT = Path("/kaggle/input")

print("Available Kaggle input folders:")
for p in sorted(INPUT_ROOT.iterdir()):
    if p.is_dir():
        print("-", p.name)

preferred = INPUT_ROOT / "rbd-dms-effect-prediction"

if preferred.exists():
    BASE_DIR = preferred
else:
    candidate_dirs = [
        p for p in INPUT_ROOT.iterdir()
        if p.is_dir() and list(p.rglob("*.csv"))
    ]
    if not candidate_dirs:
        raise FileNotFoundError(
            "No Kaggle input dataset folder containing CSV files was found. "
            "Add the Kaggle Dataset to this notebook using Add Data."
        )
    BASE_DIR = candidate_dirs[0]

print("\nUsing BASE_DIR:", BASE_DIR)
print("\nFiles found:")
for p in sorted(BASE_DIR.rglob("*")):
    if p.is_file():
        print("-", p.relative_to(BASE_DIR))

# %%
def find_file(filename):
    matches = list(BASE_DIR.rglob(filename))
    if not matches:
        available_csvs = [str(p.relative_to(BASE_DIR)) for p in BASE_DIR.rglob("*.csv")]
        raise FileNotFoundError(
            f"Could not find {filename} under {BASE_DIR}.\n"
            f"Available CSV files:\n" + "\n".join(available_csvs)
        )
    return matches[0]


def read_csv_file(filename):
    path = find_file(filename)
    print(f"Reading {path.relative_to(BASE_DIR)}")
    return pd.read_csv(path)

# %% [markdown]
# ## 1. Load summary tables

# %%
background_summary = read_csv_file("background_heldout_summary_metrics.csv")
mutation_summary = read_csv_file("mutation_heldout_summary_metrics.csv")
site_summary = read_csv_file("delta_esm_repeated_site_heldout_summary_metrics.csv")
feature_importance = read_csv_file("feature_importance_grouped_permutation.csv")
top_sites = read_csv_file("top_difficult_sites.csv")

# %% [markdown]
# ## 2. Background-held-out validation
#
# Background-held-out validation tests whether mutation-effect patterns transfer across measured backgrounds.

# %%
background_summary.sort_values(
    ["target", "model", "R2_mean"],
    ascending=[True, True, False]
).head(20)

# %%
rf_bg = background_summary[background_summary["model"] == "RandomForest"].copy()
rf_bg = rf_bg.sort_values("R2_mean", ascending=False)

plt.figure(figsize=(10, 5))
plt.barh(rf_bg["target"] + " | " + rf_bg["feature_set"], rf_bg["R2_mean"])
plt.xlabel("Mean R2")
plt.ylabel("Target | Feature set")
plt.title("Background-held-out validation")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 3. Mutation-held-out validation
#
# Mutation-held-out validation tests unseen substitution identities. The same residue position may still be represented by other substitutions.

# %%
mutation_summary.sort_values(
    ["target", "model", "R2_mean"],
    ascending=[True, True, False]
).head(20)

# %%
rf_mut = mutation_summary[mutation_summary["model"] == "RandomForest"].copy()
rf_mut = rf_mut.sort_values("R2_mean", ascending=False)

plt.figure(figsize=(10, 5))
plt.barh(rf_mut["target"] + " | " + rf_mut["feature_set"], rf_mut["R2_mean"])
plt.xlabel("Mean R2")
plt.ylabel("Target | Feature set")
plt.title("Mutation-held-out validation")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Site-held-out validation
#
# Site-held-out validation is the hardest setting because entire residue positions are removed from training.

# %%
site_summary.sort_values(
    ["target", "model", "R2_mean"],
    ascending=[True, True, False]
).head(20)

# %%
rf_site = site_summary[site_summary["model"] == "RandomForest"].copy()
rf_site = rf_site.sort_values("R2_mean", ascending=False)

plt.figure(figsize=(10, 5))
plt.barh(rf_site["target"] + " | " + rf_site["feature_set"], rf_site["R2_mean"])
plt.xlabel("Mean R2")
plt.ylabel("Target | Feature set")
plt.title("Repeated site-held-out validation")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 5. Interpolation vs extrapolation comparison
#
# The key modeling issue is that high performance in easier splits does not guarantee generalization to unseen residue positions.

# %%
comparison_rows = []

for name, df in [
    ("Background-held-out", background_summary),
    ("Mutation-held-out", mutation_summary),
    ("Site-held-out", site_summary),
]:
    rf = df[df["model"] == "RandomForest"].copy()
    best = rf.sort_values("R2_mean", ascending=False).groupby("target").head(1)
    for _, row in best.iterrows():
        comparison_rows.append({
            "validation": name,
            "target": row["target"],
            "best_feature_set": row["feature_set"],
            "R2_mean": row["R2_mean"],
            "R2_std": row.get("R2_std", None),
        })

comparison = pd.DataFrame(comparison_rows)
comparison

# %%
plt.figure(figsize=(8, 5))
for target in comparison["target"].unique():
    sub = comparison[comparison["target"] == target]
    plt.plot(sub["validation"], sub["R2_mean"], marker="o", label=target)

plt.axhline(0, linestyle="--")
plt.ylabel("Best RandomForest mean R2")
plt.xlabel("Validation setting")
plt.title("Performance drops from interpolation to unseen-site extrapolation")
plt.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Feature-importance analysis
#
# Grouped permutation importance estimates how much model performance drops when a feature group is shuffled.

# %%
feature_importance.sort_values(
    ["validation", "target", "permutation_R2_drop_mean"],
    ascending=[True, True, False]
).head(30)

# %%
for validation in feature_importance["validation"].unique():
    for target in feature_importance["target"].unique():
        sub = feature_importance[
            (feature_importance["validation"] == validation) &
            (feature_importance["target"] == target)
        ].copy()

        if sub.empty:
            continue

        sub = sub.sort_values("permutation_R2_drop_mean", ascending=True)

        plt.figure(figsize=(8, 4))
        plt.barh(sub["group"], sub["permutation_R2_drop_mean"])
        plt.xlabel("Mean R2 drop after permutation")
        plt.ylabel("Feature group")
        plt.title(f"Grouped permutation importance: {validation} | {target}")
        plt.tight_layout()
        plt.show()

# %% [markdown]
# ## 7. Difficult-site analysis
#
# Site-held-out errors were not uniformly distributed across the RBD.

# %%
top_sites.head(30)

# %%
plot_sites = top_sites.sort_values("difficulty_rank").head(20).copy()

plt.figure(figsize=(10, 5))
plt.bar(plot_sites["site"].astype(str), plot_sites["combined_mae"])
plt.xlabel("Residue site")
plt.ylabel("Relative difficulty score")
plt.title("Top difficult sites under site-held-out validation")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 8. Final interpretation
#
# Background-held-out and mutation-held-out validation showed strong performance, while site-held-out validation remained difficult.
#
# This means the model can reproduce DMS patterns and interpolate mutation effects at already represented residue positions, but it does not reliably extrapolate to completely unseen residue positions.
#
# Therefore, this project should be interpreted as a DMS effect reproduction, interpolation, and model-limitation analysis rather than a fully general predictor of unseen residue effects.
