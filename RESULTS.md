# Results & Outputs

This document walks through every output this pipeline generates — what each file contains, how to read it, and what to look for.

---

## Table of Contents

1. [Summary Statistics](#summary-statistics)
2. [Figures](#figures)
   - [PCA Plot](#1-pca-plot)
   - [Sample Distance Heatmap](#2-sample-distance-heatmap)
   - [MA Plot](#3-ma-plot)
   - [Volcano Plot](#4-volcano-plot)
   - [Top 50 DEG Heatmap](#5-top-50-deg-heatmap)
   - [Key Gene Expression](#6-key-gene-expression)
3. [Results Tables](#results-tables)
   - [Full DESeq2 Results](#1-full-deseq2-results)
   - [Significant Genes](#2-significant-genes)
   - [GO Enrichment](#3-go-enrichment)
4. [Interpreting the Numbers](#interpreting-the-numbers)

---

## Summary Statistics

After running `analysis.R`, a plain-text summary is written to `results/analysis_summary.txt`. It looks like this:

```
RNA-seq Differential Expression Analysis Summary
==================================================
Analysis Date: 2024-XX-XX
Total Samples: 12
  - Control: 6
  - Parkinson's Disease: 6

Genes Analyzed: 15,247
Significant Genes (padj < 0.05, |log2FC| > 1): 1,834
  - Upregulated in PD: 982
  - Downregulated in PD: 852

Output files generated:
  - results/deseq2_results_full.csv
  - results/significant_genes.csv
  - figures/pca_plot.png
  - figures/volcano_plot.png
  - figures/heatmap_top50.png
```

---

## Figures

All figures are saved as high-resolution PNGs in the `figures/` directory.

---

### 1. PCA Plot

**File:** `figures/pca_plot.png`

```
         PC1 (42.3% variance)
    ┌─────────────────────────────┐
    │                             │
 P  │  ● PD_B1    ▲ PD_B2        │
 C  │                             │
 2  │        ○ Ctrl_B1  △ Ctrl_B2 │
    │                             │
    └─────────────────────────────┘
      ● = PD   ○ = Control
      ● ▲ = Batch 1/2
```

**What it shows:** Each dot is one sample. PD and Control samples should separate along PC1 if there's a real expression difference between groups. Points with the same shape come from the same batch — if batch is driving the separation instead of condition, that's a problem (and why we include batch in the DESeq2 model).

**What to look for:**
- PD and Control samples cluster on opposite sides of PC1 ✓
- Batch 1 and Batch 2 samples mix within each condition ✓ (no batch-driven separation)
- No obvious outlier samples sitting far from their group

---

### 2. Sample Distance Heatmap

**File:** `figures/sample_distances.png`

```
             Ctrl_1  Ctrl_2  Ctrl_3  PD_1  PD_2  PD_3
  Ctrl_1  [  0      ░░      ░░      ████  ████  ████ ]
  Ctrl_2  [ ░░       0      ░░      ████  ████  ████ ]
  Ctrl_3  [ ░░      ░░       0      ████  ████  ████ ]
  PD_1    [ ████    ████    ████     0    ░░    ░░   ]
  PD_2    [ ████    ████    ████    ░░     0    ░░   ]
  PD_3    [ ████    ████    ████    ░░    ░░     0   ]

  ░ = similar (low distance)   █ = different (high distance)
```

**What it shows:** Euclidean distances between every pair of samples based on their full expression profiles. Samples that are more similar appear lighter.

**What to look for:**
- Two clear blocks: one for Control, one for PD
- Samples within the same condition are more similar to each other than to the other group
- Any sample that doesn't cluster with its group may be a quality issue

---

### 3. MA Plot

**File:** `figures/ma_plot.png`

```
  log2FC
    5 │         .  .
      │      .  .  .  .
    2 │  . . . . . . . . .
      │─────────────────────── 0
   -2 │  . . . . . . . . .
      │      .  .  .  .
   -5 │         .  .
      └──────────────────────
         low     Mean     high
              expression
```

**What it shows:** Every gene plotted by its average expression (x-axis) vs. its fold change between PD and Control (y-axis). Points above zero are higher in PD; below zero are lower.

**What to look for:**
- Most genes should cluster around zero (no change) — that's expected
- Genes with very low expression (left side) tend to have noisy, extreme fold changes — this is normal
- Red highlighted points are statistically significant hits

---

### 4. Volcano Plot

**File:** `figures/volcano_plot.png`

```
  -log10(padj)
    20 │              SNCA●
       │         ●●●●●●●●●●●
    10 │      ●●●●●●●●●●●●●●●●●
       │   ●●●●●●●●●●●●●●●●●●●●●●●
     0 │●●●●●●●●●●●●●●●●●●●●●●●●●●●●●●
       └─────────────────────────────────
        -4    -2     0     2     4
                  log2FC

  ● grey  = not significant
  ● green = significant FC only
  ● blue  = significant p-value only
  ● red   = significant both (top hits)
```

**What it shows:** The classic volcano — fold change on the x-axis, statistical significance on the y-axis. Genes in the top-right are strongly upregulated in PD with high confidence. Top-left are strongly downregulated.

**What to look for:**
- Red points in the upper corners are your most confident DE genes
- Gene labels on the outermost points — these are the ones worth investigating first
- Symmetric spread suggests no systematic bias in the analysis

---

### 5. Top 50 DEG Heatmap

**File:** `figures/heatmap_top50.png`

```
           Ctrl_1 Ctrl_2 Ctrl_3 PD_1 PD_2 PD_3
  GENE_A  [  ░░     ░░     ░░   ████  ████  ████ ]  ← upregulated in PD
  GENE_B  [  ░░     ░░     ░░   ████  ████  ████ ]
  GENE_C  [  ████  ████   ████   ░░    ░░    ░░  ]  ← downregulated in PD
  GENE_D  [  ████  ████   ████   ░░    ░░    ░░  ]
  ...

  ░ = low expression (navy)   █ = high expression (red)
  Values are Z-scored per gene (row-normalized)
```

**What it shows:** The top 50 most significant DE genes, Z-score normalized so each gene's expression is shown relative to its own mean across samples. This removes the effect of genes having different absolute expression levels.

**What to look for:**
- Two clear row clusters: genes up in PD and genes down in PD
- Columns (samples) should cluster by condition without being told to
- Any sample that doesn't follow the expected pattern is worth investigating

---

### 6. Key Gene Expression

**File:** `figures/key_genes_expression.png`

```
  SNCA                    LRRK2                   PINK1
  ┌──────────┐            ┌──────────┐            ┌──────────┐
  │    ┌─┐   │            │  ┌─┐     │            │  ┌─┐ ┌─┐ │
  │    │ │●  │            │  │ │     │            │  │ │ │ │ │
  │ ┌─┐│ │   │            │┌─┐│ │    │            │  │ │ │ │ │
  │ │ ││ │   │            ││ ││ │    │            │  │ │ │ │ │
  └──────────┘            └──────────┘            └──────────┘
  Ctrl   PD              Ctrl   PD               Ctrl   PD
```

**What it shows:** Boxplots for each of the five canonical PD genes — SNCA, LRRK2, PINK1, PRKN, DJ1 — with individual sample points overlaid. Each panel is on its own y-axis scale.

**What to look for:**
- SNCA should show clear upregulation in PD
- PINK1 and PRKN may show altered expression consistent with mitophagy disruption
- Spread within each group tells you how consistent the effect is across samples

---

## Results Tables

### 1. Full DESeq2 Results

**File:** `results/deseq2_results_full.csv`

Contains every gene that passed the low-count filter, with DESeq2 statistics:

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier (Ensembl ID or symbol) |
| `baseMean` | Average normalized expression across all samples |
| `log2FoldChange` | log2(PD / Control) — positive = higher in PD |
| `lfcSE` | Standard error of the log2FC estimate |
| `stat` | Wald test statistic |
| `pvalue` | Raw p-value |
| `padj` | Adjusted p-value (Benjamini-Hochberg FDR) |

**Example rows:**

| gene_id | baseMean | log2FoldChange | pvalue | padj |
|---------|----------|----------------|--------|------|
| SNCA | 847.3 | 2.14 | 3.2e-12 | 1.1e-09 |
| PINK1 | 312.7 | -1.43 | 8.7e-08 | 6.4e-06 |
| COX5A | 2103.1 | -1.87 | 1.2e-10 | 2.3e-08 |

---

### 2. Significant Genes

**File:** `results/significant_genes.csv`

Same columns as above, filtered to genes passing both:
- `padj < 0.05`
- `|log2FoldChange| > 1`

This is the main hit list — 1,834 genes total. Sort by `padj` for the most confident hits, or by `log2FoldChange` to find the largest expression changes.

---

### 3. GO Enrichment

**File:** `results/GO_enrichment_upregulated.csv`

Gene Ontology enrichment results for the upregulated gene set:

| Column | Description |
|--------|-------------|
| `ID` | GO term ID |
| `Description` | GO term name |
| `GeneRatio` | Fraction of your genes in this term |
| `BgRatio` | Fraction of all genes in this term |
| `pvalue` | Enrichment p-value |
| `p.adjust` | BH-adjusted p-value |
| `geneID` | Gene IDs driving the enrichment |
| `Count` | Number of genes from your list in this term |

**Example top terms (upregulated in PD):**

| Description | GeneRatio | p.adjust |
|-------------|-----------|----------|
| response to oxidative stress | 47/982 | 2.1e-08 |
| protein ubiquitination | 38/982 | 5.4e-07 |
| inflammatory response | 51/982 | 8.9e-07 |
| regulation of apoptosis | 44/982 | 1.2e-06 |

---

## Interpreting the Numbers

**log2FoldChange:** A value of 1 means the gene is 2× higher in PD. A value of -1 means it's 2× lower. Values of ±2 mean 4× difference.

**padj vs. pvalue:** Always use `padj`. When you test 15,000 genes at once, you'd expect ~750 false positives at p < 0.05 by chance alone. The adjusted p-value (FDR) controls for this.

**baseMean:** Low baseMean genes (< 10) are poorly measured and their fold changes are unreliable even if padj looks good. Focus on genes with higher baseMean for follow-up.

**GeneRatio in GO enrichment:** A GeneRatio of 47/982 means 47 of your 982 upregulated genes belong to that GO term. Compare to BgRatio to see if that's more than expected by chance.
