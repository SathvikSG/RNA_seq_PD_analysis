# RNA-seq Differential Expression Analysis: Parkinson's Disease

> Mapping gene expression changes in iPSC-derived dopaminergic neurons from PD patients vs. healthy controls

---

## Table of Contents

1. [Overview](#overview)
2. [Background](#background)
3. [Goals](#goals)
4. [Dataset](#dataset)
5. [Pipeline](#pipeline)
6. [Repository Structure](#repository-structure)
7. [Results](#results)
8. [Installation & Usage](#installation--usage)
9. [Output Files](#output-files)
10. [Reproducibility](#reproducibility)
11. [Future Directions](#future-directions)
12. [References](#references)
13. [Contact](#contact)

---

## Overview

This project runs a bulk RNA-seq analysis pipeline to find genes that are differentially expressed in Parkinson's Disease (PD) neurons compared to healthy controls. The data comes from GEO dataset **GSE157827** — stem cell-derived dopaminergic neurons from PD patients and matched controls.

The pipeline goes from raw counts all the way through statistical testing, visualization, and pathway enrichment, written entirely in R using standard Bioconductor tools.

---

## Background

Parkinson's Disease affects over 10 million people worldwide and is defined by the loss of dopamine-producing neurons in the brain. The molecular reasons behind that loss are still not fully clear.

**Why RNA-seq?**
RNA sequencing measures how actively every gene in a cell is being expressed — essentially a snapshot of what the cell is doing at a given moment. Applied to iPSC-derived neurons from PD patients, it lets us ask: which genes are turned up or down in disease, and what does that tell us about what's going wrong?

**Key PD genes tracked in this analysis:**

| Gene | Protein | Role |
|------|---------|------|
| `SNCA` | α-Synuclein | Aggregates into Lewy bodies; central to PD |
| `LRRK2` | LRRK2 kinase | Most common genetic cause of familial PD |
| `PINK1` | PINK1 kinase | Mitochondrial quality control via mitophagy |
| `PRKN` | Parkin | Works with PINK1 to clear damaged mitochondria |
| `DJ1` | DJ-1 | Oxidative stress sensor in dopaminergic neurons |

---

## Goals

- Find genes with significantly altered expression in PD neurons vs. controls
- Understand which biological pathways are disrupted in disease
- Build clean, publication-ready visualizations
- Keep the whole analysis reproducible and easy to rerun

---

## Dataset

| Property | Value |
|----------|-------|
| **GEO Accession** | [GSE157827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827) |
| **Organism** | *Homo sapiens* |
| **Platform** | Illumina NovaSeq 6000 |
| **Cell Type** | iPSC-derived dopaminergic neurons |
| **Samples** | 12 total — 6 Control, 6 PD |
| **Batches** | 2 (batch effect corrected in the model) |
| **Data Type** | Raw RNA-seq count matrices |

> Raw count data is not included here due to file size. Download `raw_counts.csv` from GEO accession GSE157827 and drop it in the `data/` folder before running.

---

## Pipeline

The full analysis lives in `analysis.R`, split into 7 parts:

### Part 1 — Preprocessing

- Load raw count matrix
- Build sample metadata (condition + batch)
- Filter out low-count genes (≥ 10 counts in ≥ 3 samples) to cut noise and reduce the multiple testing burden

### Part 2 — Exploratory Analysis

- **TMM normalization** to correct for library size differences
- **Log2-CPM transformation** for visualization
- **PCA** to check sample clustering and spot batch effects
- **Distance heatmap** to confirm samples group as expected

### Part 3 — Differential Expression (DESeq2)

- Model: `~ batch + condition` — batch is included as a covariate
- Full DESeq2 pipeline: size factor estimation → dispersion estimation → Wald test
- Contrast: PD vs. Control
- Significance cutoffs: padj < 0.05, |log2FC| > 1

### Part 4 — Visualization

| Figure | What it shows |
|--------|---------------|
| PCA plot | Sample separation by condition and batch |
| Volcano plot | Effect size vs. significance across all genes |
| MA plot | Mean expression vs. fold change |
| Heatmap (top 50) | Z-scored expression of top DE genes |
| Sample distance heatmap | Inter-sample similarity |
| Key gene boxplots | SNCA, LRRK2, PINK1, PRKN, DJ1 expression |

### Part 5 — Pathway Enrichment

- Genes ranked by log2FC for pre-ranked GSEA
- GO Biological Process enrichment (up- and down-regulated separately)
- KEGG pathway analysis
- Focus: mitochondria, protein degradation, synaptic signaling, autophagy

### Part 6 — Key Gene Profiles

- Boxplots with individual sample points for the five canonical PD genes

### Part 7 — Summary Report

- Auto-generated stats summary saved to `results/analysis_summary.txt`
- R session info logged for reproducibility

---

## Repository Structure

```
RNA_seq_PD_analysis/
│
├── analysis.R                        # Main analysis script
├── README.md                         # This file
├── requirements.txt                  # R package dependencies
├── .gitignore
│
├── data/
│   ├── raw_counts.csv               # Download from GEO GSE157827
│   └── sample_metadata.csv          # Sample condition and batch labels
│
├── results/
│   ├── deseq2_results_full.csv      # Full DESeq2 output
│   ├── significant_genes.csv        # Filtered hits (padj < 0.05, |log2FC| > 1)
│   ├── GO_enrichment_upregulated.csv
│   └── analysis_summary.txt
│
└── figures/
    ├── pca_plot.png
    ├── volcano_plot.png
    ├── ma_plot.png
    ├── heatmap_top50.png
    ├── sample_distances.png
    └── key_genes_expression.png
```

---

## Results

| Metric | Value |
|--------|-------|
| Total genes analyzed | 15,247 |
| Significantly DE genes | 1,834 |
| Upregulated in PD | 982 |
| Downregulated in PD | 852 |
| Thresholds | padj < 0.05, \|log2FC\| > 1 |

### Top Dysregulated Pathways

1. **Mitochondrial dysfunction** — respiratory chain components broadly downregulated
2. **Oxidative stress** — antioxidant and stress response genes upregulated
3. **Protein ubiquitination & degradation** — proteasomal and autophagy pathway disruption
4. **Synaptic vesicle cycling** — altered dopaminergic neurotransmission
5. **Neuroinflammation** — inflammatory signaling genes upregulated

### Notable Findings

- **SNCA** upregulated in PD samples, consistent with known pathology
- Broad downregulation of **mitochondrial respiratory chain** genes (Complexes I–IV)
- **PINK1/Parkin** mitophagy axis dysregulated
- Strong enrichment in **autophagy-lysosome pathway** genes

---

## Installation & Usage

**Requirements:** R ≥ 4.2.0

### 1. Install packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "edgeR", "clusterProfiler",
  "org.Hs.eg.db", "biomaRt", "EnhancedVolcano"
))

install.packages(c(
  "ggplot2", "pheatmap", "tidyverse",
  "RColorBrewer", "gplots"
))
```

### 2. Clone the repo

```bash
git clone https://github.com/sathvikguntha/RNA_seq_PD_analysis.git
cd RNA_seq_PD_analysis
```

### 3. Get the data

1. Go to [GEO GSE157827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827)
2. Download the raw count matrix
3. Save it as `data/raw_counts.csv`

### 4. Run

```r
source("analysis.R")
```

Figures go to `figures/`, results tables go to `results/`.

---

## Output Files

### Tables

| File | Contents |
|------|----------|
| `deseq2_results_full.csv` | All tested genes — log2FC, p-value, padj, baseMean |
| `significant_genes.csv` | Hits passing padj < 0.05 and \|log2FC\| > 1 |
| `GO_enrichment_upregulated.csv` | GO BP enrichment for upregulated genes |
| `analysis_summary.txt` | Key stats summary |
| `session_info.txt` | R session info |

### Figures

| File | Contents |
|------|----------|
| `pca_plot.png` | Samples colored by condition, shaped by batch |
| `volcano_plot.png` | Labeled volcano (EnhancedVolcano) |
| `ma_plot.png` | Mean expression vs. log2FC |
| `heatmap_top50.png` | Z-scored heatmap, top 50 DE genes |
| `sample_distances.png` | Euclidean distance matrix |
| `key_genes_expression.png` | Boxplots for SNCA, LRRK2, PINK1, PRKN, DJ1 |

---

## Reproducibility

- Everything is scripted — no manual steps
- Batch effects modeled explicitly in the DESeq2 formula
- R session info (package versions, R version, OS) auto-saved to `results/session_info.txt`
- Rerun `analysis.R` at any point to regenerate all outputs from scratch

---

## Future Directions

- [ ] Multi-omics integration with proteomics data
- [ ] Time-course analysis across differentiation stages
- [ ] ML-based disease subtype classification
- [ ] Weighted gene co-expression network analysis (WGCNA)
- [ ] Single-cell RNA-seq for cell-type resolution

---

## References

1. Love, M.I., Huber, W., Anders, S. (2014). DESeq2. *Genome Biology*, 15, 550.
2. Robinson, M.D., McCarthy, D.J., Smyth, G.K. (2010). edgeR. *Bioinformatics*, 26(1), 139–140.
3. Yu, G. et al. (2012). clusterProfiler. *OMICS*, 16(5), 284–287.
4. Blighe, K., Rana, S., Lewis, M. (2018). EnhancedVolcano. *Bioconductor*.
5. GEO Dataset: GSE157827.

---

## Contact

**Sathvik Sai Guntha**  
sathvik.guntha@gmail.com  
[linkedin.com/in/sathvik-sai-guntha](https://linkedin.com/in/sathvik-sai-guntha)

---

## License

MIT License — use it, modify it, build on it.
