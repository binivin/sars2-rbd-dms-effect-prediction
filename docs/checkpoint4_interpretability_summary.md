# Checkpoint 4 Interpretability and Difficult-site Summary

## Purpose

This checkpoint interpreted the model results after the reference ESM and delta ESM analyses.

The goals were:

1. Identify which feature groups the model relies on.
2. Determine why mutation-held-out validation performs well.
3. Determine why site-held-out validation remains difficult.
4. Identify RBD residue positions with large site-held-out prediction errors.

## Feature-importance Analysis

Grouped permutation importance was used to estimate how much each feature group contributed to RandomForest prediction performance.

The feature groups were:

- `site_position`
- `aa_identity`
- `aa_property_change`
- `structure`
- `reference_esm`
- `delta_esm`
- `source_background`

## Main Feature-importance Findings

### Mutation-held-out validation

In mutation-held-out validation, reference ESM was the dominant feature group.

For ACE2 binding:

```text
reference_esm R2 drop = 1.065
delta_esm     R2 drop = 0.097
aa_property   R2 drop = 0.042
site_position R2 drop = 0.035
```

For expression:

```text
reference_esm R2 drop = 1.278
delta_esm     R2 drop = 0.076
aa_property   R2 drop = 0.044
```

Interpretation: the strong mutation-held-out performance mainly comes from reference ESM residue-context information. The model uses sequence-context embeddings more strongly than simple amino-acid property or structural features when predicting unseen substitutions at represented residue positions.

### Site-held-out validation

Site-held-out validation remained weak, but feature importance still showed a useful pattern.

For ACE2 binding:

```text
reference_esm R2 drop = 0.121
delta_esm     R2 drop = 0.054
structure     R2 drop = 0.015
site_position R2 drop = 0.006
```

Interpretation: reference ESM and delta ESM provide some information even for unseen sites, but the overall predictive performance remains low. This means that ESM-derived features are useful but not sufficient for fully general site-level extrapolation.

For expression, reference ESM and amino-acid property features were the most informative groups, while delta ESM and structure were less useful in the tested model.

## Difficult-site Analysis

The difficult-site analysis summarized site-held-out prediction errors by RBD residue position.

The most difficult sites by combined ACE2 binding and expression MAE included:

```text
454, 442, 355, 398, 379, 467, 490, 350, 423, 461
```

## ACE2-binding Difficult Sites

For ACE2 binding, many high-error sites were in or near the receptor-binding motif (RBM):

```text
490, 504, 442, 454, 487, 483, 499
```

This suggests that ACE2-binding prediction errors are enriched around receptor-binding-related regions.

However, several high-error RBM sites were not labeled as 5 Å ACE2-interface sites by the simple 6M0J distance feature. This indicates that a single static structure and a simple contact-distance threshold do not fully capture functional sensitivity in the DMS data.

## Expression/Folding Difficult Sites

For expression, high-error sites included both RBM and non-RBM positions:

```text
379, 355, 461, 398, 467, 454, 380, 430, 492, 417
```

This suggests that expression/folding effects are not limited to the ACE2-binding surface and may depend on broader RBD stability or folding-sensitive regions.

## Error Direction

The difficult-site analysis also showed that prediction bias differs by site.

Examples:

- Sites 490 and 504 had negative ACE2-binding bias, meaning the model tended to underpredict their normalized ACE2-binding score.
- Sites 442 and 454 had positive ACE2-binding bias, meaning the model tended to overpredict their normalized ACE2-binding score.

This indicates that the model is not simply globally overpredicting or underpredicting. Instead, prediction errors are site-specific.

## Main Conclusion

The interpretability analyses support the following conclusions:

1. Reference ESM is the most important feature group for known-site mutation interpolation.
2. Delta ESM provides modest additional information for unseen-site ACE2-binding prediction.
3. Simple structural features have limited influence in the current model.
4. Source/background labels are not the main driver of prediction performance after normalization.
5. Site-held-out failure is concentrated in biologically sensitive positions, especially RBM-associated ACE2-binding sites and some broader expression/folding-sensitive regions.
6. Current features improve interpolation, but they do not fully solve extrapolation to entirely unseen RBD positions.

## Updated Research Interpretation

The model should be interpreted as a DMS effect reproduction and interpolation model rather than a fully general effect predictor for unseen residue positions.

The strongest result is not that the model can predict every unseen RBD site, but that the validation framework reveals where the model generalizes and where it fails:

```text
Background-held-out: strong transfer across measured backgrounds
Mutation-held-out: strong prediction for unseen substitutions at known sites
Site-held-out: weak prediction for entirely unseen positions
Difficult-site analysis: errors enriched in RBM and folding-sensitive regions
```

## Next Step

The next step should be to prepare final report materials:

1. A clean method summary.
2. A validation-result table.
3. A feature-importance table.
4. A difficult-site table.
5. A final interpretation section explaining interpolation vs extrapolation.

Additional modeling could be explored later, but at this stage the project already has a coherent computational-biology research narrative.
