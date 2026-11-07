
# Integrative Analysis of NEET Genes in Breast Cancer

## Overview

This project investigates the expression of three NEET protein-coding genes — **CISD1, CISD2, and CISD3** — in breast cancer using publicly available TCGA-BRCA data.

The objective is to explore differences in gene expression between tumor and normal breast tissue and to investigate expression patterns across PAM50 molecular subtypes.

## Data

The analysis uses gene expression and clinical annotation data from the TCGA Breast Invasive Carcinoma (TCGA-BRCA) cohort, retrieved through UCSC Xena.

Gene expression values are provided as log2-transformed normalized counts.

## Analysis workflow

1. Data retrieval and preparation.
2. Tumor versus normal expression comparison.
3. Expression analysis across PAM50 molecular subtypes.
4. Statistical analysis using Wilcoxon, Kruskal-Wallis and Dunn tests.
5. Data visualization using boxplots.

## Genes of interest

- CISD1
- CISD2
- CISD3

## Main observations

- All three genes show higher median expression in tumor samples than in normal samples.
- Their expression distributions differ across PAM50 molecular subtypes.
- CISD3 shows the largest descriptive tumor–normal difference and the strongest subtype-associated effect size among the three genes in this analysis.

These observations are exploratory and do not establish causality or therapeutic relevance.

## Repository structure

```text
Projet_NEET/
├── data/
├── figures/
├── reports/
├── results/
├── scripts/
│   ├── 01_data_preparation.R
│   ├── 02_tumor_normal_comparison.R
│   ├── 03_PAM50_analysis.R
│   └── 04_visualization.R
└── README.md
```

## Limitations

- The analysis is based on bulk RNA-seq data.
- Gene expression does not directly measure protein abundance.
- The study is observational and exploratory.
- The PAM50 Normal-like category is a molecular subtype and should not be confused with normal breast tissue samples.

## Tools

- R
- UCSCXenaTools
- dplyr
- tidyr
- ggplot2
- FSA

## Reproducibility

Run the scripts in numerical order from the project root directory.

The first script retrieves and prepares the data. The subsequent scripts perform the statistical analyses and generate figures.

## Author

Bioinformatics and Genomics
