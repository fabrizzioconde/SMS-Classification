# SMS-Classification

**SNPs Marker Selector (SMS) for Classification**

Implementation and experimental evaluation of the SMS method applied to binary classification with imbalanced genomic data, developed as part of an undergraduate thesis (TCC) at the Universidade Federal de Juiz de Fora (UFJF).

---

## Overview

This repository contains the R code, datasets, and results for a systematic evaluation of the SNPs Marker Selector (SMS) method under different experimental conditions, including:

- Balanced and imbalanced training/validation scenarios
- SMOTE-NC oversampling within stratified cross-validation folds
- F1-Score as the fitness metric for the Genetic Algorithm (GA)
- Reproduction of Simulation 6 from Oliveira (2015) with state-of-the-art techniques

## Repository Structure

```
SMS-Classification/
│
├── simulations/                       # All simulation experiments
│   │
│   ├── SMS_Oliveira2015_Sim6_Desbal_SMOTE_F1/          # Sim 6 — single execution
│   ├── SMS_Oliveira2015_Sim6_Desbal_SMOTE_F1_10exec/   # Sim 6 — 10 independent datasets
│   ├── SMS_Oliveira2015_Sim6_Desbal_SMOTE_F1_1base10exec/  # Sim 6 — 1 dataset × 10 runs
│   │
│   ├── SMS_Completo_Acuracia_1_Base_Balanceada_Validacao_Balanceada/
│   ├── SMS_Completo_Acuracia_1_Base_Balanceada_Validacao_Desbalanceada/
│   ├── SMS_Completo_Acuracia_1_Base_Desbalanceada_Validacao_Balanceada/
│   ├── SMS_Completo_Acuracia_1_Base_Desbalanceada_Validacao_Desbalanceada/
│   ├── SMS_Completo_Acuracia_1_Desbalanceada_CV_estrat_SMOTE/
│   ├── SMS_Completo_Acuracia_2_Balanceada_Validacao_Balanceada/
│   ├── SMS_Completo_Acuracia_2_Base_Balanceada_Validacao_Desbalanceada/
│   ├── SMS_Completo_Acuracia_2_Base_Desbalanceada_Validacao_Desbalanceada/
│   ├── SMS_Completo_Acuracia_2_Desbalanceada_CV_estrat_SMOTE/
│   └── SMS_Completo_Acuracia_2_Desbalanceada_Validacao_Balanceada/
│
├── scripts/                           # Utility R and Python scripts
│   ├── gvi_pvi_todas_bases.R          # GVI/PVI analysis across all scenarios
│   ├── gerar_figura1_gwas_real.R      # GWAS figure generator
│   ├── buscar_gwas_publicacoes.R      # GWAS publication search
│   └── add_citations.py              # Bibliography utility (Python)
│
├── report/                            # LaTeX source of the TCC report
│   ├── Relatorio_Completo_SMS_TCC.tex # Main document
│   ├── APENDICES.tex                  # Appendices
│   ├── capa_izabela.tex               # Cover page
│   ├── bibliography/                  # BibTeX files
│   │   ├── Referencias.bib
│   │   └── Bibliografia.bib
│   └── figures/                       # Figures used in the report
│
├── data/                              # Auxiliary data files
│   └── gwas_publicacoes_por_ano.csv   # GWAS publications per year
│
└── docs/                              # Documentation and reference materials
    ├── articles/                      # Reference articles
    ├── report-history/                # Previous versions of the report
    ├── TCC_1_Izabela.pdf              # TCC version 1 (draft)
    └── Tabelas___TCC_II.pdf           # Summary tables (TCC II)
```

## Simulations Description

### Simulation 6 — Oliveira (2015) Reproduction

Reproduces the binary imbalanced classification scenario from Section 7.1.6 of Oliveira (2015), with the following improvements:

| Feature | Oliveira (2015) | This work |
|---|---|---|
| Oversampling | — | SMOTE-NC (train folds only) |
| Cross-validation | Random k-fold | Stratified k-fold (k=10) |
| GA fitness metric | AUC-ROC | F1-Score |
| SVM type | SVR adapted | C-classification (binary) |

**Simulated model:** n = 1000, 100 SNPs, MAF ~ U[0.1, 0.4], class ratio 862:138.

### Complete Scenarios (Accuracy)

Eleven scenarios combining:

- **Dataset version:** 1 (original SMS) or 2 (with parameter tuning)
- **Training data:** Balanced or Imbalanced
- **Validation data:** Balanced, Imbalanced, or Stratified CV with SMOTE

Each scenario folder contains:
- `*.R` — Main simulation script
- `dados*.csv` — Simulated dataset(s)
- `metricas_causais_CLASSIF.csv` — Precision/Recall/F1 of causal SNP detection
- `Grafico_*.pdf` — F1 and GA evolution curves per SVM kernel
- `rank_global_classif/` — Global classification ranking results

## Requirements

- **R** ≥ 4.0
- **R packages:** `e1071`, `GA`, `smotefamily`, `caret`, `dplyr`, `ggplot2`

Install dependencies:

```r
install.packages(c("e1071", "GA", "smotefamily", "caret", "dplyr", "ggplot2"))
```

## How to Run

Each simulation is self-contained. Open the `.R` script inside any `simulations/` subfolder in RStudio and run it with `Source` (`Ctrl+Shift+Enter`). Scripts auto-detect their working directory.

**Example — Sim 6 single execution:**

```r
# Open and source:
simulations/SMS_Oliveira2015_Sim6_Desbal_SMOTE_F1/SMS_Oliveira2015_Sim6_Desbal_SMOTE_F1.R
```

## Reference

> OLIVEIRA, Fabrízzio Condé de. *Um método para seleção de atributos em dados genômicos*. Tese (Doutorado em Modelagem Computacional) — Universidade Federal de Juiz de Fora, 2015.

## Author

**Izabela** — TCC II, UFJF  
Advisor: Prof. Fabrízzio Condé de Oliveira

---

*This project is part of ongoing research on genomic SNP selection for classification at UFJF.*
