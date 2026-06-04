# Final Project Summary

## Project Title

ESM-2 and Structure-based RBD DMS Effect Prediction: Interpolation and Generalization Limits

## Core Research Question

This project asked whether amino-acid substitution effects measured by RBD deep mutational scanning can be reproduced using mutation-level, structure-level, and protein-language-model features.

A second and more important question was whether high prediction performance reflects true generalization or only interpolation within already represented residue positions.

## Data and Targets

The project used public RBD deep mutational scanning datasets and built a unified substitution-level table.

The processed dataset contained:

```text
22,572 amino-acid substitutions
5 measured backgrounds: Wuhan-Hu-1, Delta, Beta, E484K, and N501Y
2 prediction targets: receptor-binding score and expression/folding score
```

Both targets were normalized within each source/background group so that model evaluation focused on relative mutation effects rather than raw score-scale differences across datasets.

## Feature Sets

The project progressively added the following feature groups:

1. Mutation-level physicochemical features
   - residue position
   - wild-type and mutant amino-acid identity
   - hydropathy, volume, and charge changes
   - proline/cysteine indicators

2. Structure-derived features from 6M0J
   - receptor-binding-motif annotation
   - minimum distance to receptor chain
   - 5 angstrom interface indicator
   - number of receptor-contact residues
   - structure-site availability

3. Reference ESM-2 residue embeddings
   - residue-level sequence-context embeddings from `facebook/esm2_t6_8M_UR50D`

4. Mutant-sequence delta ESM features
   - mutation-induced representation change:

```text
delta_embedding = mutant_sequence_embedding - background_wildtype_sequence_embedding
```

## Validation Strategy

The most important part of the project was the validation design.

### Random split

This is the easiest setting. Similar sites and substitutions can appear in both training and test sets, so high performance mainly indicates within-distribution interpolation.

### Background-held-out split

One measured background is held out at a time. This tests whether mutation-effect patterns transfer across backgrounds.

### Mutation-held-out split

A subset of mutation identities is held out across all backgrounds. This tests whether the model can predict unseen substitutions at positions that may still be represented by other substitutions.

### Site-held-out split

Entire residue positions are held out from training. This is the hardest setting and tests whether the model can generalize to completely unseen residue positions.

## Main Results

### Background-held-out validation

Background-held-out validation showed strong performance.

```text
RandomForest with ESM features:
Binding R2_mean    about 0.91
Expression R2_mean about 0.94
```

This indicates that the model can transfer learned mutation-effect patterns across measured backgrounds.

### Mutation-held-out validation

Mutation-held-out validation also showed strong performance.

```text
RandomForest with reference ESM + structure + background:
Binding R2_mean    about 0.81
Expression R2_mean about 0.81
```

This indicates that the model can predict unseen substitutions reasonably well when the residue position is already represented in the training data.

### Site-held-out validation

Site-held-out validation remained difficult.

```text
RandomForest with delta ESM + structure + background:
Binding R2_mean    about 0.03
Expression R2_mean about 0.12
```

Delta ESM produced a modest improvement in the hardest setting, but the overall predictive performance remained weak.

## Feature-importance Findings

Grouped permutation importance showed that reference ESM was the dominant feature group in mutation-held-out validation.

For mutation-held-out binding prediction:

```text
reference_esm R2 drop = 1.065
delta_esm     R2 drop = 0.097
aa_property   R2 drop = 0.042
site_position R2 drop = 0.035
```

For mutation-held-out expression prediction:

```text
reference_esm R2 drop = 1.278
delta_esm     R2 drop = 0.076
aa_property   R2 drop = 0.044
```

This suggests that strong mutation-held-out performance mainly comes from sequence-context information captured by reference ESM-2 embeddings.

In site-held-out validation, reference ESM and delta ESM still provided some information, but they were not sufficient for robust unseen-site prediction.

## Difficult-site Findings

Difficult-site analysis showed that site-held-out errors were not uniformly distributed across the RBD.

The most difficult sites by combined binding and expression error included:

```text
454, 442, 355, 398, 379, 467, 490, 350, 423, 461
```

For binding prediction, high-error sites were enriched around receptor-binding-motif-associated positions such as:

```text
490, 504, 442, 454, 487, 483, 499
```

For expression/folding prediction, difficult sites included both receptor-binding-motif and non-receptor-binding-motif positions:

```text
379, 355, 461, 398, 467, 454, 380, 430, 492, 417
```

This suggests that binding errors are concentrated near receptor-binding-sensitive regions, while expression/folding errors reflect broader stability or folding-sensitive regions.

## Final Interpretation

The most important conclusion is:

```text
The model predicts new substitutions at already represented residue positions fairly well, but it does not reliably predict mutation effects at completely unseen residue positions.
```

Therefore, this project should be interpreted as a DMS effect reproduction and interpolation study, not as a fully general predictor of unseen residue effects.

The strength of the project is not simply high model performance. Its main value is that it separates easy interpolation settings from difficult extrapolation settings and identifies where the model fails.

## Limitations

1. The structure features were based on a single static complex structure and simple distance thresholds.
2. The ESM experiment used a small ESM-2 model; larger models may change performance.
3. Reference ESM embeddings were strong for represented positions but did not solve unseen-site prediction.
4. Delta ESM improved site-held-out prediction only modestly.
5. The model is not suitable for claiming reliable prediction at entirely unseen residue positions.

## Future Work

Future analyses could improve unseen-site generalization by adding:

- larger protein language model embeddings
- multiple structural conformations
- solvent accessibility
- residue conservation scores
- structural dynamics features
- more detailed interface annotations
- focused analysis of high-error residue positions

## Final Takeaway

This project demonstrates that careful validation design is essential in biological effect prediction. A model can appear highly accurate under random or mutation-held-out validation but still fail when asked to predict entirely unseen residue positions.
