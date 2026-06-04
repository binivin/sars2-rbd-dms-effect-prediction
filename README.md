# RBD DMS Effect Prediction

This repository analyzes public deep mutational scanning data for the Spike receptor-binding domain. The goal is to evaluate whether sequence, structure, and protein-language-model features can reproduce experimentally measured DMS scores for receptor binding and protein expression or folding.

## Project Status

Current checkpoint: **final project summary**

Completed steps:

1. Downloaded public RBD DMS datasets.
2. Built a unified DMS table.
3. Removed non-substitution rows where the wild-type and mutant amino acids were identical.
4. Normalized binding and expression scores within each source/background group.
5. Trained baseline models using mutation-level physicochemical features.
6. Added structure-derived residue features from the 6M0J receptor-binding complex.
7. Added ESM-2 reference residue-level embeddings using `facebook/esm2_t6_8M_UR50D`.
8. Evaluated random split, site-held-out, repeated site-held-out, background-held-out, and mutation-held-out validation settings.
9. Added mutant-sequence delta ESM features and compared them with reference ESM features.
10. Performed grouped feature-importance and difficult-site analyses.
11. Wrote the final project summary.

## Dataset Summary

The processed substitution dataset contains:

- 22,572 RBD amino-acid substitutions
- 5 measured backgrounds: Wuhan-Hu-1, Delta, Beta, E484K, and N501Y
- 2 prediction targets:
  - receptor-binding score
  - protein expression/folding score

Generated files are excluded from GitHub by `.gitignore` and are recreated by running the pipeline.

## Main Findings

The validation results show a clear difference between interpolation and generalization difficulty.

### Easier validation settings

Background-held-out and mutation-held-out validation showed strong performance, especially when ESM-2 features were used.

- Background-held-out, RandomForest, ESM features:
  - binding: mean R2 about 0.91
  - expression: mean R2 about 0.94
- Mutation-held-out, RandomForest, reference ESM + structure + background features:
  - binding: mean R2 about 0.81
  - expression: mean R2 about 0.81

These results suggest that the model can predict new substitutions within already represented residue positions and transfer reasonably well across measured backgrounds.

### Hardest validation setting

Repeated site-held-out validation remained difficult.

- Binding, delta ESM + structure + background, RandomForest: mean R2 about 0.03
- Expression, delta ESM + structure + background, RandomForest: mean R2 about 0.12

This suggests that mutant-sequence delta ESM features provide a modest improvement for unseen-site prediction, but full generalization to completely unseen residue positions remains limited.

### Feature importance and difficult sites

Grouped permutation importance showed that reference ESM features dominate mutation-held-out prediction, while delta ESM provides modest information for site-held-out binding prediction.

Difficult-site analysis showed that site-held-out errors were not uniformly distributed across the RBD. Binding errors were enriched around receptor-binding-motif-associated positions, while expression/folding errors appeared across both receptor-binding-motif and non-receptor-binding-motif sites.

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
src/analyze_difficult_sites.py
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
python src/analyze_difficult_sites.py
```

## Documentation

Result summaries are stored in `docs/`.

```text
docs/baseline_result_summary.md
docs/structural_feature_result_summary.md
docs/checkpoint2_validation_summary.md
docs/checkpoint3_delta_esm_summary.md
docs/checkpoint4_interpretability_summary.md
docs/final_project_summary.md
```

## Final Takeaway

This project demonstrates that careful validation design is essential in biological effect prediction. A model can appear highly accurate under random or mutation-held-out validation but still fail when asked to predict entirely unseen residue positions.
