# Kaggle Upload Guide

This folder contains Kaggle-ready text and metadata for publishing this project as a Kaggle Dataset and Kaggle Notebook.

## Recommended Kaggle Dataset Title

```text
RBD DMS Effect Prediction and Generalization Analysis
```

## Recommended Dataset Subtitle

```text
Processed result tables and summaries for RBD DMS score prediction using mutation, structure, reference ESM-2, and delta ESM features.
```

## Recommended Kaggle Dataset Description

This dataset contains generated summary tables and documentation from a computational biology project analyzing RBD deep mutational scanning scores.

The project evaluates whether mutation-level physicochemical features, structure-derived annotations, reference ESM-2 residue embeddings, and mutant-sequence delta ESM features can reproduce experimentally measured receptor-binding and expression/folding scores.

The key focus is not only model performance, but also validation design. The analysis compares random split, background-held-out, mutation-held-out, and site-held-out validation to distinguish easy interpolation from difficult extrapolation.

## Main Finding

The model performs well when predicting new substitutions at already represented residue positions, but generalization to completely unseen residue positions remains difficult.

## Recommended Files to Upload to Kaggle Dataset

Upload generated result files from your local project folder after running the pipeline.

Recommended folders/files:

```text
artifacts/tables/rbd_baseline_normalized_model_metrics.csv
artifacts/tables/rbd_baseline_structural_model_metrics.csv
artifacts/tables/rbd_baseline_esm_model_metrics.csv
artifacts/tables/repeated_site_heldout_summary_metrics.csv
artifacts/tables/background_heldout_summary_metrics.csv
artifacts/tables/mutation_heldout_summary_metrics.csv
artifacts/tables/delta_esm_repeated_site_heldout_summary_metrics.csv
artifacts/tables/delta_esm_mutation_heldout_summary_metrics.csv
artifacts/tables/feature_importance_grouped_permutation.csv
artifacts/tables/top_difficult_sites.csv
docs/final_project_summary.md
docs/project_summary_ko.md
```

Do not upload very large intermediate ESM embedding tables unless needed, because Kaggle pages are easier to review when the dataset contains concise result tables and summaries.

## Kaggle CLI Workflow

Install Kaggle CLI:

```bash
pip install kaggle
```

Create a Kaggle dataset folder locally, copy the generated tables and docs into it, then place `dataset-metadata.json` in the same folder.

Create the dataset:

```bash
kaggle datasets create -p kaggle_dataset_folder
```

Update the dataset later:

```bash
kaggle datasets version -p kaggle_dataset_folder -m "Update RBD DMS validation summaries"
```

## Recommended Kaggle Notebook Title

```text
RBD DMS Effect Prediction: Interpolation vs Extrapolation
```

## Recommended Notebook Outline

1. Project overview
2. Dataset summary
3. Validation strategy
4. Background-held-out results
5. Mutation-held-out results
6. Site-held-out results
7. Feature-importance analysis
8. Difficult-site analysis
9. Final interpretation

## Important Interpretation

This project should be interpreted as a DMS effect reproduction, interpolation, and model-limitation analysis. It should not be presented as a fully general predictor of unseen residue effects.
