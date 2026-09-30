# bioinfo-portfolio
Transcriptomic Data Analysis with R & Python
# Breast Cancer Transcriptomics: ER+ vs. ER- Differential Expression Analysis

## Overview
This project presents an end-to-end transcriptomic analysis identifying Differentially Expressed Genes (DEGs) between **Estrogen Receptor-positive (ER+)** and **Estrogen Receptor-negative (ER-)** breast cancer subtypes. Using RNA-Seq count data, the goal is to elucidate the transcriptional downstream targets driven by estrogen signaling and characterize the molecular footprint of ER status.

---

## Key Findings & Biological Insights
* **Target Gene Downregulation:** Primary DEGs identified (including `GENE_17`, `GENE_24`, and `GENE_25`) demonstrate significant downregulation in ER- condition ($\log_2\text{FC} \approx -2.1$, $p_{\text{adj}} < 10^{-50}$).
* **Coordinated Transcriptional Module:** The heatmaps and volcano plots confirm a tight, co-regulated gene block whose suppression characterizes the ER- phenotype, reflecting the loss of ESR1-mediated transcriptional activation.

---

## Technical & Statistical Pipeline
The analysis is implemented in **R** using Bioconductor packages following standard differential expression guidelines:

1. **Preprocessing & Filtering:** Low-count genes with total counts $< 10$ across samples were removed to reduce multiple testing burden.
2. **Model Fitting (`DESeq2`):** Size-factor normalization, empirical Bayes dispersion estimation, and Generalized Linear Model (GLM) fitting using the Negative Binomial distribution.
3. **Hypothesis Testing:** Wald test for differential expression; $p$-values adjusted for Multiple Testing using the Benjamini-Hochberg (FDR) procedure.
4. **Data Visualization:** High-resolution generation of Volcano plots (`EnhancedVolcano`) and hierarchical clustering heatmaps (`pheatmap`).

---

## Repository Structure
```text
breast-cancer-rnaseq/
├── data/               # Raw count matrix and metadata
├── scripts/            # R script (01_breast_cancer_deseq2.R)
├── results/            # Output CSV (DEGs) and PNG visualizations
└── README.md           # Project documentation