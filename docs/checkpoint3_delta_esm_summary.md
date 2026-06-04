# Checkpoint 3 Delta ESM Summary

## Purpose

This checkpoint tested whether mutant-sequence ESM delta embeddings improve prediction of normalized RBD DMS scores.

Checkpoint 2 used reference residue-level ESM-2 embeddings. Those features describe the sequence context of each residue in the reference/background sequence, but they do not directly encode how a specific amino-acid substitution changes the protein-language-model representation.

Checkpoint 3 therefore tested a delta representation:

```text
delta_embedding = mutant_sequence_embedding - background_wildtype_sequence_embedding
```

The goal was to evaluate whether mutation-induced ESM representation changes help difficult validation settings, especially repeated site-held-out validation.

## Features Compared

The main feature sets were:

- `mutation_only`
- `mutation_structure`
- `mutation_structure_refesm_background`
- `mutation_structure_deltaesm_background`
- `mutation_structure_refesm_deltaesm_background`

The prediction targets were:

- `ace2_binding_bg_z`
- `expression_bg_z`

## Repeated Site-held-out Results

Site-held-out validation remains the hardest setting because entire residue positions are removed from training.

For ACE2 binding with RandomForest:

```text
mutation_only:                          R2_mean = -0.053
mutation_structure:                     R2_mean = -0.129
mutation_structure_refesm_background:   R2_mean = -0.029
mutation_structure_deltaesm_background: R2_mean =  0.033
```

For expression with RandomForest:

```text
mutation_only:                          R2_mean = 0.061
mutation_structure:                     R2_mean = 0.080
mutation_structure_deltaesm_background: R2_mean = 0.117
```

Interpretation: delta ESM features produced the best repeated site-held-out performance among the tested feature sets. The improvement was modest, but it suggests that mutation-induced ESM representation changes contain useful information for unseen-site generalization.

## Mutation-held-out Results

Mutation-held-out validation tests unseen substitution identities while residue positions may still be represented by other substitutions.

For ACE2 binding with RandomForest:

```text
mutation_structure_refesm_background:   R2_mean = 0.805
mutation_structure_deltaesm_background: R2_mean = 0.598
```

For expression with RandomForest:

```text
mutation_structure_refesm_background:   R2_mean = 0.812
mutation_structure_deltaesm_background: R2_mean = 0.664
```

Interpretation: reference ESM features remained stronger for known-site interpolation. Delta ESM features did not improve mutation-held-out performance and, in this setting, reduced performance compared with reference ESM.

## Main Conclusion

Delta ESM and reference ESM behave differently depending on the validation setting.

```text
Reference ESM:
Strong for known-site mutation interpolation.
Best in mutation-held-out validation.

Delta ESM:
Modestly helpful for unseen-site generalization.
Best in repeated site-held-out validation among tested feature sets.
```

Therefore, the current conclusion is:

> Delta ESM improves the hardest site-held-out setting slightly, but it does not solve unseen-site prediction. Reference ESM remains better for predicting unseen substitutions at already represented residue positions.

## Updated Research Interpretation

The full validation pattern now suggests that the model has strong interpolation ability but limited extrapolation ability.

- Background-held-out: strong transfer across measured backgrounds
- Mutation-held-out: strong prediction for unseen substitutions at known sites
- Site-held-out: weak prediction for entirely unseen residue positions
- Delta ESM: small improvement for site-held-out, but still limited

This supports the broader interpretation that RBD DMS effects are strongly residue-position dependent.

## Next Step

The next checkpoint should focus on interpretability rather than adding more feature types immediately.

Recommended next analysis:

```text
Feature importance and model limitation analysis
```

Possible tasks:

1. Compare feature importance for RandomForest models.
2. Identify which validation settings are driven by site, amino-acid property changes, structure, reference ESM, or delta ESM.
3. Analyze difficult residue positions in site-held-out validation.
4. Produce figures and tables for a final report or presentation.
