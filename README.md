# RBD DMS Effect Prediction and Generalization Analysis

[![Python](https://img.shields.io/badge/Python-3.x-blue?logo=python&logoColor=white)](https://www.python.org/)
[![Field](https://img.shields.io/badge/Field-Bioinformatics-green)](#)
[![ML](https://img.shields.io/badge/Method-ESM--2%20%7C%20RandomForest%20%7C%20Ridge-purple)](#)
[![Kaggle Dataset](https://img.shields.io/badge/Kaggle-Dataset-20BEFF?logo=kaggle&logoColor=white)](https://www.kaggle.com/datasets/binivin/rbd-dms-effect-prediction)
[![Kaggle Notebook](https://img.shields.io/badge/Kaggle-Notebook-20BEFF?logo=kaggle&logoColor=white)](https://www.kaggle.com/code/binivin/rbd-dms-interpolation-vs-extrapolation)
[![Status](https://img.shields.io/badge/Status-Final%20summary%20complete-brightgreen)](#)

> Deep mutational scanning score prediction using mutation features, structure-derived annotations, reference ESM-2 embeddings, delta ESM embeddings, and validation-split analysis.

This repository contains a computational biology pipeline for analyzing RBD deep mutational scanning data. The project evaluates whether amino-acid substitution effects on receptor binding and protein expression/folding can be reproduced using sequence, structure, and protein-language-model features.

The main finding is that the model performs well when predicting new substitutions at already represented residue positions, but generalization to completely unseen residue positions remains difficult.

---

## Kaggle

| Resource | Link |
|---|---|
| Kaggle Dataset | https://www.kaggle.com/datasets/binivin/rbd-dms-effect-prediction |
| Kaggle Notebook | https://www.kaggle.com/code/binivin/rbd-dms-interpolation-vs-extrapolation |

The Kaggle Dataset contains the final result summary tables and project summaries. The Kaggle Notebook provides a concise visualization-oriented version of the analysis, focusing on validation splits, feature importance, and difficult-site interpretation.

---

## Highlights

- Built a unified RBD DMS substitution dataset with **22,572 amino-acid substitutions**.
- Normalized binding and expression/folding scores within each source/background group.
- Engineered mutation-level physicochemical features.
- Added structure-derived residue annotations from the **6M0J** receptor-binding complex.
- Added **reference ESM-2 residue embeddings** using `facebook/esm2_t6_8M_UR50D`.
- Added **mutant-sequence delta ESM features** to represent mutation-induced embedding changes.
- Compared random split, background-held-out, mutation-held-out, and site-held-out validation.
- Performed grouped feature-importance analysis.
- Identified difficult residue positions under site-held-out validation.

---

## Project workflow

```text
Public RBD DMS data
  -> unified substitution dataset construction
  -> source/background-level score normalization
  -> mutation-level feature engineering
  -> structure feature annotation from 6M0J
  -> reference ESM-2 residue embedding
  -> mutant-sequence delta ESM embedding
  -> baseline model training
  -> validation split comparison
  -> grouped feature-importance analysis
  -> difficult-site analysis
  -> final interpretation
```

---

## Data sources

The analysis uses public RBD deep mutational scanning datasets and a public protein complex structure.

| Source | Role in this project |
|---|---|
| RBD DMS datasets | Experimental binding and expression/folding scores |
| 6M0J structure | Structure-derived receptor-interface and distance features |
| ESM-2 | Protein-language-model residue embeddings |
| Kaggle Dataset | Final result tables and summary documents |

Generated intermediate files and model output tables are excluded from GitHub by `.gitignore`. The concise final result tables are available through Kaggle.

---

## Main results

### 1. Random split performance can overestimate generalization

Random split validation produced high performance because similar residue positions and substitutions can appear in both training and test sets. Therefore, random split performance was treated as a basic reproduction check rather than the main evidence of generalization.

---

### 2. Background-held-out validation showed strong transfer across measured backgrounds

When one measured background was held out at a time, models with ESM features retained strong performance.

| Target | Best model setting | Mean R2 |
|---|---|---:|
| Binding | RandomForest with ESM features | ~0.91 |
| Expression/folding | RandomForest with ESM features | ~0.94 |

This suggests that mutation-effect patterns can transfer across the measured backgrounds in this dataset.

---

### 3. Mutation-held-out validation showed strong known-site interpolation

Mutation-held-out validation removed specific substitution identities from training. The best models still performed well when the residue position was represented by other substitutions.

| Target | Best model setting | Mean R2 |
|---|---|---:|
| Binding | RandomForest with reference ESM + structure + background | ~0.81 |
| Expression/folding | RandomForest with reference ESM + structure + background | ~0.81 |

This indicates that the model can predict unseen substitutions at already represented residue positions.

---

### 4. Site-held-out validation remained difficult

Site-held-out validation removed entire residue positions from training. This was the hardest and most biologically important validation setting.

| Target | Best model setting | Mean R2 |
|---|---|---:|
| Binding | RandomForest with delta ESM + structure + background | ~0.03 |
| Expression/folding | RandomForest with delta ESM + structure + background | ~0.12 |

Delta ESM features modestly improved unseen-site prediction, but the overall performance remained limited.

---

### 5. Reference ESM was the dominant feature group for mutation-held-out prediction

Grouped permutation importance showed that reference ESM features were the most important feature group in mutation-held-out validation.

| Validation | Target | Most important feature group | Interpretation |
|---|---|---|---|
| Mutation-held-out | Binding | Reference ESM | Sequence-context information drives known-site interpolation |
| Mutation-held-out | Expression/folding | Reference ESM | Residue context strongly supports prediction |
| Site-held-out | Binding | Reference ESM + delta ESM | ESM features help, but do not fully solve unseen-site extrapolation |

---

### 6. Difficult-site analysis revealed position-specific model failures

The most difficult sites by combined binding and expression/folding error included:

```text
454, 442, 355, 398, 379, 467, 490, 350, 423, 461
```

Binding prediction errors were enriched around receptor-binding-motif-associated positions such as:

```text
490, 504, 442, 454, 487, 483, 499
```

Expression/folding errors appeared across both receptor-binding-motif and non-receptor-binding-motif positions, suggesting that expression/folding effects may reflect broader stability-sensitive regions.

---

## Interpretation

The validation results support the following interpretation:

> The model can reproduce DMS patterns and interpolate mutation effects at already represented residue positions, but it does not reliably extrapolate to completely unseen residue positions.

This means the project should be interpreted as a DMS effect reproduction, interpolation, and model-limitation analysis rather than a fully general predictor of unseen residue effects.

The main value of this project is not only the model performance itself, but the validation design that separates easy interpolation from difficult extrapolation.

---

## Repository structure

```text
sars2-rbd-dms-effect-prediction/
├─ README.md
├─ requirements.txt
├─ src/
│  ├─ build_rbd_dms_dataset.py
│  ├─ make_rbd_substitutions.py
│  ├─ train_rbd_baseline.py
│  ├─ train_rbd_baseline_normalized.py
│  ├─ add_rbd_structural_features.py
│  ├─ train_rbd_baseline_with_structure.py
│  ├─ run_esm_embedding_baseline.py
│  ├─ repeat_site_heldout_validation.py
│  ├─ run_background_heldout_validation.py
│  ├─ run_mutation_heldout_validation.py
│  ├─ run_delta_esm_validation.py
│  └─ analyze_difficult_sites.py
├─ docs/
│  ├─ baseline_result_summary.md
│  ├─ structural_feature_result_summary.md
│  ├─ checkpoint2_validation_summary.md
│  ├─ checkpoint3_delta_esm_summary.md
│  ├─ checkpoint4_interpretability_summary.md
│  ├─ final_project_summary.md
│  └─ project_summary_ko.md
├─ kaggle/
│  ├─ README.md
│  ├─ dataset-metadata.json
│  ├─ kernel-metadata.json
│  └─ rbd_dms_effect_prediction_interpolation_vs_extrapolation.py
├─ data/                  # local generated data, excluded from GitHub
└─ artifacts/             # generated intermediate outputs, mostly excluded from GitHub
```

---

## How to reproduce

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Run the main analysis pipeline

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

The ESM-related scripts may take longer on CPU because they run protein-language-model inference.

### 3. View the Kaggle summary notebook

A visualization-oriented summary is available here:

```text
https://www.kaggle.com/code/binivin/rbd-dms-interpolation-vs-extrapolation
```

---

## Documentation

| Document | Description |
|---|---|
| `docs/baseline_result_summary.md` | Baseline and normalization summary |
| `docs/structural_feature_result_summary.md` | Structure feature experiment summary |
| `docs/checkpoint2_validation_summary.md` | Reference ESM and validation split analysis |
| `docs/checkpoint3_delta_esm_summary.md` | Delta ESM analysis |
| `docs/checkpoint4_interpretability_summary.md` | Feature-importance and difficult-site analysis |
| `docs/final_project_summary.md` | Final English project summary |
| `docs/project_summary_ko.md` | Korean project summary |
| `kaggle/README.md` | Kaggle dataset and notebook guide |

---

## Limitations

- The structure features rely on a single static complex structure and simple distance thresholds.
- The ESM experiment used a small ESM-2 model; larger models may change performance.
- Reference ESM embeddings are strong for represented positions but do not solve unseen-site prediction.
- Delta ESM improves site-held-out prediction only modestly.
- The model should not be interpreted as a reliable predictor for entirely unseen residue positions.

---

## Next steps

- Prepare final report figures and validation-result tables.
- Add visualization scripts for model-comparison plots.
- Analyze high-error residue positions in more biological detail.
- Test larger protein language models if computational resources allow.
- Add conservation, solvent accessibility, or multi-structure features.
- Convert the pipeline into a concise research report and presentation.
