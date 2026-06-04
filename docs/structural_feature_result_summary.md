# Structural Feature Baseline Result Summary

## Purpose

This checkpoint tested whether simple structure-derived residue features improve generalization of RBD DMS score prediction, especially under site-held-out validation.

The structural features were generated from the RBD-ACE2 complex structure 6M0J and included:

- Rough RBM annotation
- Minimum distance to ACE2
- Whether the residue is within 5 Å of ACE2
- Number of ACE2 contact residues within 5 Å
- Whether the residue is available in the structure

## Feature Sets Compared

Four feature sets were tested:

1. `mutation_only`
2. `mutation_plus_background`
3. `mutation_plus_structure`
4. `mutation_plus_structure_background`

The prediction targets were source/background-normalized DMS scores:

- `ace2_binding_bg_z`
- `expression_bg_z`

## Key Results

### Random split

Random split performance remained high for RandomForest models across feature sets.

For example:

```text
ACE2 binding, RandomForest, mutation_only:
R2 = 0.900

ACE2 binding, RandomForest, mutation_plus_structure_background:
R2 = 0.920

Expression, RandomForest, mutation_only:
R2 = 0.929

Expression, RandomForest, mutation_plus_structure_background:
R2 = 0.938
```

This indicates that the model can interpolate DMS effects within already represented residue positions.

### Site-held-out split

Site-held-out performance remained much weaker, but simple structure features improved some models slightly.

For ACE2 binding:

```text
RandomForest, mutation_only:
R2 = -0.174

RandomForest, mutation_plus_structure:
R2 = -0.127

Ridge, mutation_only:
R2 = -0.008

Ridge, mutation_plus_structure:
R2 = 0.025
```

For expression:

```text
RandomForest, mutation_only:
R2 = 0.100

RandomForest, mutation_plus_structure:
R2 = 0.134

Ridge, mutation_only:
R2 = 0.117

Ridge, mutation_plus_structure:
R2 = 0.197
```

## Interpretation

Simple structure-derived features provide a modest improvement, especially for the Ridge model and for expression prediction under site-held-out validation.

However, RandomForest site-held-out performance remains weak, and ACE2 binding prediction still shows negative R2 even after adding structural features. This means that the current structural features are not sufficient for robust generalization to unseen RBD positions.

The key conclusion is:

> Simple structural context helps slightly, but unseen-site generalization remains limited.

## Next Step

The next step should add stronger sequence-context features, especially protein language model embeddings, and compare them against the current mutation-only and structure-enhanced baselines.

A reasonable next checkpoint is:

```text
Add ESM-2 residue-level embeddings and test whether they improve site-held-out performance.
```
