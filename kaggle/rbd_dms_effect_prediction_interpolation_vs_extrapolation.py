# %% [markdown]
# # RBD DMS Effect Prediction: Interpolation vs Extrapolation
#
# This Kaggle notebook summarizes a computational biology project on RBD deep mutational scanning (DMS) effect prediction.
#
# The goal is to evaluate whether mutation-level features, structure-derived annotations, reference ESM-2 embeddings, and mutant-sequence delta ESM features can reproduce experimentally measured receptor-binding and expression/folding scores.
#
# The main conclusion is that the model performs well for interpolation at already represented residue positions, but extrapolation to completely unseen residue positions remains difficult.

# %% [markdown]
# ## 1. Load summary tables
#
# This notebook is designed to run on a Kaggle Dataset containing generated summary tables from the GitHub pipeline.
#
# Expected files include:
#
# - `background_heldout_summary_metrics.csv`
# - `mutation_heldout_summary_metrics.csv`
# - `delta_esm_repeated_site_heldout_summary_metrics.csv`
# - `delta_esm_mutation_heldout_summary_metrics.csv`
# - `feature_importance_grouped_permutation.csv`
# - `top_difficult_sites.csv`

# %%
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path('/kaggle/input/rbd-dms-effect-prediction-generalization-analysis')

print('Input directory exists:', BASE_DIR.exists())
if BASE_DIR.exists():
    print('Available files:')
    for p in sorted(BASE_DIR.rglob('*')):
        if p.is_file():
            print('-', p.relative_to(BASE_DIR))

# %%
def find_file(filename):
    matches = list(BASE_DIR.rglob(filename))
    if not matches:
        raise FileNotFoundError(f'Could not find {filename} under {BASE_DIR}')
    return matches[0]


def read_csv_if_exists(filename):
    path = find_file(filename)
    print(f'Reading {path.relative_to(BASE_DIR)}')
    return pd.read_csv(path)

# %%
background_summary = read_csv_if_exists('background_heldout_summary_metrics.csv')
mutation_summary = read_csv_if_exists('mutation_heldout_summary_metrics.csv')
site_summary = read_csv_if_exists('delta_esm_repeated_site_heldout_summary_metrics.csv')
delta_mutation_summary = read_csv_if_exists('delta_esm_mutation_heldout_summary_metrics.csv')
feature_importance = read_csv_if_exists('feature_importance_grouped_permutation.csv')
top_sites = read_csv_if_exists('top_difficult_sites.csv')

# %% [markdown]
# ## 2. Background-held-out validation
#
# Background-held-out validation tests whether mutation-effect patterns transfer across measured backgrounds.

# %%
background_summary.sort_values(['target', 'model', 'R2_mean'], ascending=[True, True, False]).head(20)

# %%
rf_bg = background_summary[background_summary['model'] == 'RandomForest'].copy()
rf_bg = rf_bg.sort_values('R2_mean', ascending=False)

plt.figure(figsize=(10, 5))
plt.barh(rf_bg['target'] + ' | ' + rf_bg['feature_set'], rf_bg['R2_mean'])
plt.xlabel('Mean R2')
plt.ylabel('Target | Feature set')
plt.title('Background-held-out validation')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 3. Mutation-held-out validation
#
# Mutation-held-out validation tests unseen substitution identities. The same residue position may still be represented by other substitutions.

# %%
mutation_summary.sort_values(['target', 'model', 'R2_mean'], ascending=[True, True, False]).head(20)

# %%
rf_mut = mutation_summary[mutation_summary['model'] == 'RandomForest'].copy()
rf_mut = rf_mut.sort_values('R2_mean', ascending=False)

plt.figure(figsize=(10, 5))
plt.barh(rf_mut['target'] + ' | ' + rf_mut['feature_set'], rf_mut['R2_mean'])
plt.xlabel('Mean R2')
plt.ylabel('Target | Feature set')
plt.title('Mutation-held-out validation')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Site-held-out validation
#
# Site-held-out validation is the hardest setting because entire residue positions are removed from training.

# %%
site_summary.sort_values(['target', 'model', 'R2_mean'], ascending=[True, True, False]).head(20)

# %%
rf_site = site_summary[site_summary['model'] == 'RandomForest'].copy()
rf_site = rf_site.sort_values('R2_mean', ascending=False)

plt.figure(figsize=(10, 5))
plt.barh(rf_site['target'] + ' | ' + rf_site['feature_set'], rf_site['R2_mean'])
plt.xlabel('Mean R2')
plt.ylabel('Target | Feature set')
plt.title('Repeated site-held-out validation')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 5. Interpolation vs extrapolation comparison
#
# The key biological modeling issue is that high performance in easier splits does not guarantee generalization to unseen residue positions.

# %%
comparison_rows = []

for name, df in [
    ('background-held-out', background_summary),
    ('mutation-held-out', mutation_summary),
    ('site-held-out', site_summary),
]:
    rf = df[df['model'] == 'RandomForest'].copy()
    best = rf.sort_values('R2_mean', ascending=False).groupby('target').head(1)
    for _, row in best.iterrows():
        comparison_rows.append({
            'validation': name,
            'target': row['target'],
            'best_feature_set': row['feature_set'],
            'R2_mean': row['R2_mean'],
            'R2_std': row.get('R2_std', None),
        })

comparison = pd.DataFrame(comparison_rows)
comparison

# %%
plt.figure(figsize=(8, 5))
for target in comparison['target'].unique():
    sub = comparison[comparison['target'] == target]
    plt.plot(sub['validation'], sub['R2_mean'], marker='o', label=target)
plt.axhline(0, linestyle='--')
plt.ylabel('Best RandomForest mean R2')
plt.xlabel('Validation setting')
plt.title('Performance drops from interpolation to unseen-site extrapolation')
plt.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Feature-importance analysis
#
# Grouped permutation importance estimates how much model performance drops when a feature group is shuffled.

# %%
feature_importance.sort_values(['validation', 'target', 'permutation_R2_drop_mean'], ascending=[True, True, False]).head(30)

# %%
for validation in feature_importance['validation'].unique():
    for target in feature_importance['target'].unique():
        sub = feature_importance[
            (feature_importance['validation'] == validation) &
            (feature_importance['target'] == target)
        ].copy()
        if sub.empty:
            continue
        sub = sub.sort_values('permutation_R2_drop_mean', ascending=True)
        plt.figure(figsize=(8, 4))
        plt.barh(sub['group'], sub['permutation_R2_drop_mean'])
        plt.xlabel('Mean R2 drop after permutation')
        plt.ylabel('Feature group')
        plt.title(f'Grouped permutation importance: {validation} | {target}')
        plt.tight_layout()
        plt.show()

# %% [markdown]
# ## 7. Difficult-site analysis
#
# Site-held-out errors were not uniformly distributed across the RBD. This table lists positions with the largest combined binding and expression/folding error.

# %%
top_sites.head(30)

# %%
plot_sites = top_sites.head(20).copy()

plt.figure(figsize=(10, 5))
plt.bar(plot_sites['site'].astype(str), plot_sites['combined_mae'])
plt.xlabel('Residue site')
plt.ylabel('Combined MAE')
plt.title('Top difficult sites under site-held-out validation')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 8. Final interpretation
#
# The validation results support the following interpretation:
#
# - Background-held-out and mutation-held-out validation show strong performance.
# - Site-held-out validation remains difficult.
# - Reference ESM features dominate known-site interpolation.
# - Delta ESM features provide modest help for unseen-site prediction.
# - Difficult sites are enriched around receptor-binding-sensitive and folding-sensitive regions.
#
# Therefore, this project should be interpreted as a DMS effect reproduction, interpolation, and model-limitation analysis rather than a fully general predictor of unseen residue effects.
