# SARS-CoV-2 RBD DMS Effect Prediction

This repository analyzes public deep mutational scanning (DMS) data for the SARS-CoV-2 Spike receptor-binding domain (RBD).

The current goal is **not** to design, recommend, or optimize new viral variants. Instead, this project evaluates whether simple mutation-level features can reproduce experimentally measured DMS scores for ACE2 binding and RBD expression/folding.

## Project Status

Current checkpoint: **baseline DMS score reproduction analysis**

Completed steps:

1. Downloaded public RBD DMS datasets.
2. Built a unified RBD DMS table.
3. Removed non-substitution rows where the wild-type and mutant amino acids were identical.
4. Trained baseline models on raw ACE2 binding and expression scores.
5. Normalized scores within each source/background group.
6. Compared random split and site-held-out split performance.

## Dataset Summary

The processed substitution dataset contains:

- 22,572 RBD amino-acid substitutions
- 5 major backgrounds: Wuhan-Hu-1, Delta, Beta, E484K, and N501Y
- 2 prediction targets:
  - ACE2 binding score
  - RBD expression/folding score

Generated file:

```text
data/processed/rbd_dms_substitutions_v0.csv
```

Generated files are excluded from GitHub by `.gitignore`.

## Baseline Modeling

Two baseline experiments were performed.

### 1. Raw-score baseline

Raw DMS scores were used directly as prediction targets. Random forest models produced very high scores, but this was likely inflated by score-scale differences between datasets and backgrounds.

### 2. Normalized-score baseline

ACE2 binding and expression scores were normalized within each source/background group using z-scores.

After normalization:

- Random split performance remained high.
- Site-held-out performance dropped substantially.

This suggests that the model can interpolate mutation effects at familiar RBD positions, but simple features are not sufficient for strong generalization to unseen RBD sites.

## Key Finding

The current baseline model mainly learns patterns within already observed RBD positions. When entire residue positions are held out, prediction performance decreases sharply.

Therefore, the next step is to add residue-level biological context, such as structural features or protein language model embeddings, and test whether these features improve site-held-out generalization.

## Main Scripts

```text
src/build_rbd_dms_dataset.py
src/make_rbd_substitutions.py
src/train_rbd_baseline.py
src/train_rbd_baseline_normalized.py
```

## How to Run

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the pipeline:

```bash
python src/build_rbd_dms_dataset.py
python src/make_rbd_substitutions.py
python src/train_rbd_baseline.py
python src/train_rbd_baseline_normalized.py
```

## Current Interpretation

The first baseline results should be interpreted as a preliminary reproducibility and validation checkpoint, not as a final biological prediction model.

The next planned step is to add structural context and evaluate whether it improves generalization under site-held-out validation.
