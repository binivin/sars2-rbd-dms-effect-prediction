# Checkpoint 2 Validation Summary

## Purpose

This checkpoint evaluated whether structure-derived features and ESM-2 residue-level embeddings improve prediction of normalized SARS-CoV-2 RBD DMS scores.

The prediction targets were:

- `ace2_binding_bg_z`: source/background-normalized ACE2 binding score
- `expression_bg_z`: source/background-normalized RBD expression/folding score

The major question was whether the model can generalize beyond simple random split performance.

## Features Compared

The experiments compared combinations of the following feature groups:

- Mutation-level features: site, wild-type amino acid, mutant amino acid, physicochemical changes
- Structure features: RBM annotation, ACE2 interface distance, 5 Å contact annotation from 6M0J
- ESM-2 reference residue embeddings: 320-dimensional embeddings from `facebook/esm2_t6_8M_UR50D`
- Background metadata: source and variant background labels

## Validation Settings

Four validation types were considered.

### 1. Random split

This is the easiest setting because similar sites and mutations can appear in both training and test sets. It mainly checks whether the model can reproduce patterns within the same data distribution.

### 2. Site-held-out split

Entire RBD residue positions are held out from training. This tests whether the model can predict mutation effects at residue positions it has never seen.

This was the hardest validation setting.

### 3. Background-held-out split

One variant background is held out at a time. This tests whether mutation-effect patterns transfer across measured backgrounds.

### 4. Mutation-held-out split

A subset of mutation identities, such as `N501Y`, is held out across all backgrounds. This tests whether the model can predict unseen substitutions at residue positions that may still be represented by other substitutions.

## Key Results

### Repeated site-held-out validation

Repeated site-held-out validation showed weak generalization to unseen residue positions.

For ACE2 binding, RandomForest results were approximately:

```text
mutation_only:          R2_mean = -0.047
mutation_structure:     R2_mean = -0.127
mutation_esm:           R2_mean = -0.043
mutation_structure_esm: R2_mean = -0.017
```

For expression, RandomForest results were approximately:

```text
mutation_only:          R2_mean = 0.066
mutation_structure:     R2_mean = 0.077
mutation_esm:           R2_mean = -0.130
mutation_structure_esm: R2_mean = -0.128
```

Interpretation: ESM-2 features slightly improved ACE2 binding in the hardest split, but the overall site-held-out performance remained weak.

### Background-held-out validation

Background-held-out validation showed strong performance.

Best RandomForest results:

```text
ACE2 binding, mutation_esm:
R2_mean = 0.907

Expression, mutation_esm:
R2_mean = 0.941
```

Interpretation: the model can transfer mutation-effect patterns across measured variant backgrounds when the mutation positions and substitutions are broadly represented in training.

### Mutation-held-out validation

Mutation-held-out validation also showed strong performance.

Best RandomForest results:

```text
ACE2 binding, mutation_structure_esm_background:
R2_mean = 0.807

Expression, mutation_structure_esm_background:
R2_mean = 0.813
```

Interpretation: the model can predict unseen amino-acid substitutions at known or represented residue positions reasonably well.

## Main Conclusion

The model shows strong interpolation ability but limited unseen-site generalization.

In practical terms:

```text
Random split: easy, high performance
Background-held-out: strong transfer across measured backgrounds
Mutation-held-out: strong prediction for unseen substitutions at known sites
Site-held-out: weak prediction for completely unseen residue positions
```

Therefore, the current model should be interpreted as a DMS score reproduction and interpolation model, not as a fully general residue-effect predictor.

## Scientific Interpretation

The results suggest that RBD DMS effects are highly position-dependent. Once a residue position is represented in training, the model can learn useful mutation-level patterns. However, when the entire residue position is unseen, simple physicochemical features, structural annotations, and reference residue-level ESM-2 embeddings are not sufficient for robust prediction.

## Next Experiment

The next experiment should test mutant-sequence ESM delta embeddings.

The current ESM feature represents the reference residue context only. It does not directly encode how the sequence representation changes after each mutation. A stronger approach is to compute:

```text
delta_embedding = mutant_sequence_embedding - wildtype_sequence_embedding
```

This can test whether mutation-induced changes in protein language model representation improve difficult validation settings, especially site-held-out prediction.
