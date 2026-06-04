# SARS-CoV-2 RBD DMS Effect Prediction

This repository analyzes public deep mutational scanning (DMS) data for the SARS-CoV-2 Spike receptor-binding domain (RBD).

The current goal is **not** to design, recommend, or optimize new viral variants. Instead, this project evaluates whether sequence-, structure-, and protein-language-model-based features can reproduce experimentally measured DMS scores for ACE2 binding and RBD expression/folding.

## Project Status

Current checkpoint: **mutant-sequence delta ESM analysis**

Completed steps:

1. Downloaded public RBD DMS datasets.
2. Built a unified RBD DMS table.
3. Removed non-substitution rows where the wild-type and mutant amino acids were identical.
4. Normalized ACE2 binding and expression scores within each source/background group.
5. Trained baseline models using mutation-level physicochemical features.
6. Added structure-derived residue features from the 6M0J RBD-ACE2 complex.
7. Added ESM-2 reference residue-level embeddings using `facebook/esm2_t6_8M_UR50D`.
8. Evaluated multiple validation settings: random split, site-held-out, repeated site-held-out, background-held-out, and mutation-held-out.
9. Added mutant-sequence delta ESM features and compared them with reference ESM features.

## Dataset Summary

The processed substitution dataset contains:

- 22,572 RBD amino-acid substitutions
- 5 major backgrounds: Wuhan-Hu-1, Delta, Beta, E484K, and N501Y
- 2 prediction targets:
  - ACE2 binding score
  - RBD expression/folding score

Generated files are excluded from GitHub by `.gitignore` and are recreated by running the pipeline.

## Main Findings

The validation results show a clear difference between interpolation and generalization difficulty.

### Easier validation settings

Background-held-out and mutation-held-out validation showed strong performance, especially when ESM-2 features were used.

- Background-held-out, RandomForest, ESM features:
  - ACE2 binding: mean R2 about 0.91
  - Expression: mean R2 about 0.94
- Mutation-held-out, RandomForest, reference ESM + structure + background features:
  - ACE2 binding: mean R2 about 0.81
  - Expression: mean R2 about 0.81

These results suggest that the model can predict new substitutions within already represented residue positions and transfer reasonably well across measured variant backgrounds.

### Hardest validation setting

Repeated site-held-out validation remained difficult.

- ACE2 binding, delta ESM + structure + background, RandomForest: mean R2 about 0.03
- Expression, delta ESM + structure + background, RandomForest: mean R2 about 0.12

This suggests that mutant-sequence delta ESM features provide a modest improvement for unseen-site prediction, but full generalization to completely unseen RBD residue positions remains limited.

## Current Interpretation

The model is strongest when the residue position has already been represented in the training data. It can learn mutation-level patterns and transfer them across backgrounds or unseen substitutions at known sites.

Reference ESM features are strongest for known-site interpolation, especially mutation-held-out validation. Delta ESM features are more useful for the hardest site-held-out setting, but the improvement is still small.

The current checkpoint supports the following conclusion:

> ESM-2 and structural features improve mutation-level interpolation, while mutant-sequence delta ESM features slightly improve unseen-site generalization. However, unseen-site prediction remains the main limitation.

## Main Scripts

```text
src/build_rbd_dms_dataset.py
src/make_rbd_substitutions.py
src/train_rbd_baseline.py
src/train_rbd_baseline_normalized.py
src/add_rbd_structural_features.py
src/train_rbd_baseline_with_structure.py
src/run_esm_embedding_baseline.py
src/repeat_site_heldout_validation.py
src/run_background_heldout_validation.py
src/run_mutation_heldout_validation.py
src/run_delta_esm_validation.py
```

## How to Run

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the main pipeline:

```bash
python src/build_rbd_dms_dataset.py
python src/make_rbd_substitutions.py
python src/train_rbd_baseline_normalized.py
python src/add_rbd_structural_features.py
python src/train_rbd_baseline_with_structure.py
python src/run_esm_embedding_baseline.py
python src/repeat_site_heldout_validation.py
python src/run_background_heldout_validation.py
python src/run_mutation_heldout_validation.py
python src/run_delta_esm_validation.py
```

## Documentation

Result summaries are stored in `docs/`.

```text
docs/baseline_result_summary.md
docs/structural_feature_result_summary.md
docs/checkpoint2_validation_summary.md
docs/checkpoint3_delta_esm_summary.md
```

## Next Step

The next planned checkpoint is feature-importance and model-limitation analysis.

Recommended next tasks:

1. Compare feature importance across mutation, structure, reference ESM, and delta ESM features.
2. Identify difficult residue positions in repeated site-held-out validation.
3. Summarize which validation settings represent interpolation and which represent extrapolation.
4. Prepare figures and tables for a final report or presentation.
