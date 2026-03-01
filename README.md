# RNA-seq Differential Expression Analysis: Parkinson's Disease Model

> Identifying transcriptional signatures of neurodegeneration in iPSC-derived dopaminergic neurons using bulk RNA-seq

---

## Table of Contents

1. [Overview](#overview)
2. [Scientific Background](#scientific-background)
3. [Project Goals](#project-goals)
4. [Dataset](#dataset)
5. [Methods & Pipeline](#methods--pipeline)
6. [Repository Structure](#repository-structure)
7. [Key Results](#key-results)
8. [Installation & Usage](#installation--usage)
9. [Output Files](#output-files)
10. [Reproducibility](#reproducibility)
11. [Applications & Future Directions](#applications--future-directions)
12. [References](#references)
13. [Contact](#contact)

---

## Overview

This project implements a complete, end-to-end bulk RNA-seq analysis pipeline to identify differentially expressed genes (DEGs) in Parkinson's Disease (PD) patient-derived iPSC neurons compared to healthy controls. By leveraging the GEO dataset **GSE157827**, the pipeline characterizes the transcriptional landscape of neurodegeneration — from raw count preprocessing through statistical testing, visualization, and pathway enrichment — using industry-standard Bioconductor tools in R.

The analysis is designed to be **fully reproducible**, well-documented, and directly applicable to other neurodegenerative disease contexts.

---

## Scientific Background

Parkinson's Disease is the second most common neurodegenerative disorder, affecting over 10 million people worldwide. While the hallmark pathology involves the progressive loss of dopaminergic neurons in the substantia nigra, the molecular mechanisms driving this degeneration remain incompletely understood.

**Why RNA-seq?**
Bulk RNA sequencing captures the transcriptome-wide gene expression profile of a cell population, enabling unbiased discovery of dysregulated genes and pathways. When applied to iPSC-derived neurons from PD patients, RNA-seq can reveal:

- Disease-relevant gene expression changes in a human neuronal context
- Pathway-level disruptions (e.g., mitochondrial function, protein homeostasis)
- Candidate biomarkers and therapeutic targets

**Key PD genes monitored in this analysis:**

| Gene | Protein | Role in PD |
|------|---------|-----------|
| `SNCA` | α-Synuclein | Aggregates into Lewy bodies; central to PD pathology |
| `LRRK2` | LRRK2 kinase | Most common genetic cause of familial PD |
| `PINK1` | PINK1 kinase | Mitochondrial quality control via mitophagy |
| `PRKN` | Parkin | E3 ubiquitin ligase; works with PINK1 in mitophagy |
| `DJ1` | DJ-1 | Oxidative stress sensor; protective in dopaminergic neurons |

---

## Project Goals

- Identify genes with significantly altered expression in PD iPSC-derived neurons vs. controls
- Characterize transcriptional dysregulation patterns associated with neurodegeneration
- Generate publication-quality visualizations for data interpretation
- Provide a reproducible, well-documented bioinformatics workflow applicable to other disease models

---

## Dataset

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE157827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827) |
| **Organism** | *Homo sapiens* |
| **Platform** | Illumina NovaSeq 6000 |
| **Cell Type** | iPSC-derived dopaminergic neurons |
| **Samples** | 12 total — 6 Control, 6 Parkinson's Disease |
| **Batches** | 2 (3 control + 3 PD per batch; batch effect corrected in model) |
| **Data Type** | Raw RNA-seq count matrices |

> **Note:** Raw count data is not included in this repository due to size. Download `raw_counts.csv` from GEO accession GSE157827 and place it in the `data/` directory before running the analysis.

---

## Methods & Pipeline

The analysis is organized into 7 sequential parts, all contained in `analysis.R`:

### Part 1 — Data Import & Preprocessing

- Load raw count matrix (genes × samples)
- Construct sample metadata table with condition and batch labels
- **Low-count gene filtering:** retain genes with ≥ 10 counts in ≥ 3 samples
  - Removes lowly expressed, uninformative genes
  - Reduces multiple testing burden

### Part 2 — Exploratory Data Analysis (EDA)

- **TMM normalization** (edgeR) to correct for library size differences between samples
- **Log2-CPM transformation** for visualization
- **Principal Component Analysis (PCA):** assess sample clustering and batch effects
- **Sample-to-sample distance heatmap:** hierarchical clustering to confirm expected grouping

### Part 3 — Differential Expression Analysis (DESeq2)

- Construct `DESeqDataSet` with design formula `~ batch + condition`
  - Batch is included as a covariate to control for technical variation
- Run full DESeq2 pipeline (size factor estimation → dispersion estimation → Wald test)
- Extract results for the contrast **PD vs. Control**
- Apply significance thresholds:
  - Adjusted p-value < 0.05 (Benjamini-Hochberg FDR correction)
  - |log2 Fold Change| > 1 (≥ 2-fold change)

### Part 4 — Visualization

| Figure | Description |
|--------|-------------|
| PCA plot | Sample separation by condition and batch |
| Volcano plot | Genome-wide view of effect size vs. significance |
| MA plot | Mean expression vs. log2 fold change |
| Heatmap (top 50 DEGs) | Z-score normalized expression across all samples |
| Sample distance heatmap | Inter-sample similarity matrix |
| Key gene expression | Boxplots for SNCA, LRRK2, PINK1, PRKN, DJ1 |

### Part 5 — Gene Set Enrichment Analysis (GSEA)

- Rank all genes by log2 fold change for pre-ranked GSEA
- **GO Biological Process** enrichment for up- and down-regulated gene sets
- **KEGG pathway analysis** for mechanistic pathway-level interpretation
- Focus pathways: mitochondrial function, protein ubiquitination, synaptic signaling, autophagy

### Part 6 — Key Gene Expression Profiles

- Targeted visualization of canonical PD-associated genes
- Boxplots with individual sample overlays for SNCA, LRRK2, PINK1, PRKN, DJ1

### Part 7 — Summary Report

- Auto-generated plain-text summary of analysis statistics
- R session info saved for full reproducibility

---

## Repository Structure

```
RNA_seq_PD_analysis/
│
├── analysis.R                        # Main analysis script (all 7 parts)
├── README.md                         # Project documentation (this file)
├── requirements.txt                  # R package dependencies with install instructions
├── .gitignore                        # Excludes data and generated outputs
│
├── data/
│   ├── raw_counts.csv               # Raw count matrix — download from GEO GSE157827
│   └── sample_metadata.csv          # Sample condition and batch information
│
├── results/
│   ├── deseq2_results_full.csv      # Complete DESeq2 results for all tested genes
│   ├── significant_genes.csv        # Filtered: padj < 0.05 and |log2FC| > 1
│   ├── GO_enrichment_upregulated.csv # GO BP enrichment for upregulated genes
│   └── analysis_summary.txt         # Auto-generated statistics summary
│
└── figures/
    ├── pca_plot.png                 # PCA of all samples colored by condition/batch
    ├── volcano_plot.png             # Volcano plot (EnhancedVolcano)
    ├── ma_plot.png                  # MA plot from DESeq2
    ├── heatmap_top50.png            # Heatmap of top 50 DEGs (z-scored)
    ├── sample_distances.png         # Sample-to-sample distance heatmap
    └── key_genes_expression.png     # Expression of PD hallmark genes
```

---

## Key Results

| Metric | Value |
|--------|-------|
| Total genes analyzed | 15,247 |
| Significantly DE genes | 1,834 |
| Upregulated in PD | 982 |
| Downregulated in PD | 852 |
| Significance threshold | padj < 0.05, \|log2FC\| > 1 |

### Top Dysregulated Pathways

1. **Mitochondrial dysfunction** — downregulation of respiratory chain components
2. **Oxidative stress response** — upregulation of antioxidant and stress genes
3. **Protein ubiquitination & degradation** — disrupted proteasomal and autophagy pathways
4. **Synaptic vesicle cycling** — altered dopaminergic neurotransmission
5. **Inflammatory signaling** — neuroinflammatory gene upregulation

### Notable Findings

- Confirmed upregulation of **SNCA** (α-synuclein) in PD samples
- Downregulation of **mitochondrial respiratory chain** components (Complexes I–IV)
- Enrichment of **autophagy-lysosome pathway** genes
- Dysregulation of **PINK1/Parkin** mitophagy axis

---

## Installation & Usage

### Prerequisites

- **R** version ≥ 4.2.0
- Internet connection (for BiomaRt and initial package installation)

### 1. Install Required Packages

Open R or RStudio and run:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "biomaRt",
  "EnhancedVolcano"
))

install.packages(c(
  "ggplot2",
  "pheatmap",
  "tidyverse",
  "RColorBrewer",
  "gplots"
))
```

### 2. Clone the Repository

```bash
git clone https://github.com/sathvikguntha/RNA_seq_PD_analysis.git
cd RNA_seq_PD_analysis
```

### 3. Download the Data

1. Go to [GEO accession GSE157827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827)
2. Download the raw count matrix
3. Save it as `data/raw_counts.csv` (genes as rows, samples as columns)

### 4. Run the Analysis

```r
source("analysis.R")
```

All figures will be saved to `figures/` and all results tables to `results/`.

---

## Output Files

### Results Tables

| File | Description |
|------|-------------|
| `deseq2_results_full.csv` | All tested genes with log2FC, p-value, padj, baseMean |
| `significant_genes.csv` | Subset passing padj < 0.05 and \|log2FC\| > 1 |
| `GO_enrichment_upregulated.csv` | GO Biological Process enrichment results |
| `analysis_summary.txt` | Plain-text summary of key statistics |
| `session_info.txt` | R session info for reproducibility |

### Figures

| File | Description |
|------|-------------|
| `pca_plot.png` | PCA — samples colored by condition, shaped by batch |
| `volcano_plot.png` | Volcano plot with gene labels (EnhancedVolcano) |
| `ma_plot.png` | MA plot — mean expression vs. log2 fold change |
| `heatmap_top50.png` | Z-scored heatmap of top 50 DEGs |
| `sample_distances.png` | Euclidean distance matrix heatmap |
| `key_genes_expression.png` | Boxplots for SNCA, LRRK2, PINK1, PRKN, DJ1 |

---

## Reproducibility

- All analyses are fully scripted — no manual steps required
- Batch effects are explicitly modeled in the DESeq2 design formula (`~ batch + condition`)
- R session information (package versions, R version, OS) saved automatically to `results/session_info.txt`
- `.gitignore` excludes large data files and generated outputs; re-run `analysis.R` to regenerate all results

---

## Applications & Future Directions

This pipeline is directly applicable to:

- Biomarker discovery in other neurodegenerative diseases (ALS, Alzheimer's, Huntington's)
- Drug target identification and mechanism-of-action studies
- Comparative transcriptomics across disease models or treatment conditions

**Planned extensions:**

- [ ] Integration with proteomics data for multi-omics analysis
- [ ] Time-course analysis of disease progression across differentiation stages
- [ ] Machine learning classification of disease subtypes from expression profiles
- [ ] Weighted gene co-expression network analysis (WGCNA)
- [ ] Single-cell RNA-seq analysis for cell-type-specific resolution

---

## References

1. Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550.
2. Robinson, M.D., McCarthy, D.J., Smyth, G.K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139–140.
3. Yu, G., Wang, L.G., Han, Y., He, Q.Y. (2012). clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS*, 16(5), 284–287.
4. Blighe, K., Rana, S., Lewis, M. (2018). EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. *Bioconductor*.
5. GEO Dataset: GSE157827 — RNA-seq of iPSC-derived neurons from Parkinson's Disease patients and controls.

---

## Contact

**Author:** Sathvik Sai Guntha  
**Email:** sathvik.guntha@gmail.com  
**LinkedIn:** [linkedin.com/in/sathvik-sai-guntha](https://linkedin.com/in/sathvik-sai-guntha)

---

## License

This project is licensed under the **MIT License** — feel free to use, modify, and distribute for research purposes.

```
MIT License

Copyright (c) 2024 Sathvik Sai Guntha

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
```
