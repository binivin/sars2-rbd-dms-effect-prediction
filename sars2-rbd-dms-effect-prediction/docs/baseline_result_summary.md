# Baseline Result Summary

## Purpose

This checkpoint tested whether simple mutation-level features can reproduce experimentally measured RBD DMS scores.

The prediction targets were:

- ACE2 binding
- RBD expression/folding

The analysis was performed twice:

1. Raw-score baseline
2. Source/background-normalized baseline

## Processed Dataset

After removing non-substitution rows, the processed dataset contained:

```text
Filtered shape: (22572, 12)
```

Background counts:

```text
Wuhan-Hu-1    7495
Delta         3818
Beta          3801
E484K         3741
N501Y         3717
```

## Raw-score Baseline

Raw-score models showed very high performance.

However, raw ACE2 binding and expression scores had different scales across datasets and backgrounds. Therefore, this result may be inflated by source/background-specific score-scale differences.

## Normalized Baseline

After source/background-level z-score normalization, random split performance remained high, but site-held-out performance dropped.

### Normalized Random Split

```text
ACE2 binding, RandomForest:
R2 = 0.900

Expression, RandomForest:
R2 = 0.929
```

### Normalized Site-held-out Split

```text
ACE2 binding, RandomForest:
R2 = -0.164

Expression, RandomForest:
R2 = 0.143
```

## Interpretation

The baseline model can reproduce DMS scores for mutations at familiar RBD positions, but it does not generalize well to unseen residue positions.

This suggests that simple features such as position, amino-acid property changes, and background labels are not sufficient for robust unseen-site generalization.

## Next Step

The next analysis should add residue-level context, such as:

- RBD/RBM annotation
- ACE2 interface proximity
- structural accessibility
- protein language model embeddings

The key question for the next experiment is:

> Do structural or sequence-embedding features improve site-held-out generalization?
