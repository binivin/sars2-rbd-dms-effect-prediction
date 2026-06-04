# Kaggle Materials

This folder contains the Kaggle companion materials for this project.

## Published Kaggle Links

| Resource | Link |
|---|---|
| Kaggle Dataset | https://www.kaggle.com/datasets/binivin/rbd-dms-effect-prediction |
| Kaggle Notebook | https://www.kaggle.com/code/binivin/rbd-dms-interpolation-vs-extrapolation |

## Kaggle Dataset

Recommended dataset metadata:

```json
{
  "title": "RBD DMS Effect Prediction",
  "id": "binivin/rbd-dms-effect-prediction",
  "licenses": [
    {
      "name": "CC0-1.0"
    }
  ]
}
```

The Kaggle Dataset contains concise result tables and summary documents rather than large intermediate model files.

Recommended uploaded files:

```text
background_heldout_summary_metrics.csv
mutation_heldout_summary_metrics.csv
delta_esm_repeated_site_heldout_summary_metrics.csv
delta_esm_mutation_heldout_summary_metrics.csv
feature_importance_grouped_permutation.csv
top_difficult_sites.csv
final_project_summary.md
project_summary_ko.md
README.md
dataset-metadata.json
```

## Kaggle Notebook

Notebook title:

```text
RBD DMS Interpolation vs Extrapolation
```

Notebook link:

```text
https://www.kaggle.com/code/binivin/rbd-dms-interpolation-vs-extrapolation
```

The notebook summarizes:

1. Background-held-out validation
2. Mutation-held-out validation
3. Site-held-out validation
4. Interpolation vs extrapolation comparison
5. Grouped feature-importance analysis
6. Difficult-site analysis
7. Final interpretation

## Important Interpretation

This project should be interpreted as a DMS effect reproduction, interpolation, and model-limitation analysis.

The main conclusion is:

```text
The model can predict new substitutions at already represented residue positions fairly well,
but it does not reliably generalize to completely unseen residue positions.
```

Therefore, high performance under easier validation settings should not be presented as evidence of fully general mutation-effect prediction.
