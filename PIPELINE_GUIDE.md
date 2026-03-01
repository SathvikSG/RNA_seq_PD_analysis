# Pipeline Guide

A step-by-step walkthrough of the RNA-seq differential expression analysis — from raw sequencing data to biological interpretation. Written to be readable whether you're new to bioinformatics or just want to understand what each step is actually doing.

---

## The Big Picture

```
┌─────────────────────────────────────────────────────────────────────┐
│                     RNA-seq Analysis Pipeline                       │
│                                                                     │
│   RAW DATA                                                          │
│      │                                                              │
│      ▼                                                              │
│  ┌─────────┐    Filter low-count genes                             │
│  │ PART 1  │    Build sample metadata                              │
│  │  PREP   │                                                        │
│  └────┬────┘                                                        │
│       │                                                             │
│       ▼                                                             │
│  ┌─────────┐    Normalize (TMM)                                    │
│  │ PART 2  │    PCA + distance heatmap                             │
│  │   EDA   │    → Are samples behaving as expected?               │
│  └────┬────┘                                                        │
│       │                                                             │
│       ▼                                                             │
│  ┌─────────┐    DESeq2 statistical model                           │
│  │ PART 3  │    Test every gene: PD vs. Control                   │
│  │   DE    │    → Which genes are significantly different?         │
│  └────┬────┘                                                        │
│       │                                                             │
│       ▼                                                             │
│  ┌─────────┐    Volcano, MA, heatmap, boxplots                    │
│  │ PART 4  │                                                        │
│  │  PLOTS  │    → Visualize the results                           │
│  └────┬────┘                                                        │
│       │                                                             │
│       ▼                                                             │
│  ┌─────────┐    GO enrichment + KEGG pathways                     │
│  │ PART 5  │                                                        │
│  │  GSEA   │    → What biology is being disrupted?                │
│  └────┬────┘                                                        │
│       │                                                             │
│       ▼                                                             │
│  ┌─────────┐    Targeted plots for SNCA, LRRK2, PINK1...          │
│  │ PART 6  │                                                        │
│  │KEY GENES│                                                        │
│  └────┬────┘                                                        │
│       │                                                             │
│       ▼                                                             │
│  ┌─────────┐    Summary stats + session info                      │
│  │ PART 7  │                                                        │
│  │ REPORT  │                                                        │
│  └─────────┘                                                        │
│                                                                     │
│   OUTPUTS: figures/ + results/                                      │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Part 1 — Data Import & Preprocessing

### What's happening

Before any analysis, we need to load the data and clean it up. RNA-seq produces a count matrix — a table where every row is a gene, every column is a sample, and each cell is the number of sequencing reads that mapped to that gene in that sample.

```
           Ctrl_1  Ctrl_2  Ctrl_3  PD_1  PD_2  PD_3
  SNCA        45      52      48    189   201   178
  LRRK2      312     298     321    287   301   294
  PINK1      104      98     112     61    58    67
  OBSCURE_1    0       1       0      0     0     1   ← filtered out
  OBSCURE_2    2       0       1      3     1     0   ← filtered out
```

### Low-count filtering

Genes with very few reads across samples are essentially unmeasured — they carry no useful signal and just add noise to the statistical tests. We remove any gene that doesn't have at least **10 counts in at least 3 samples**.

```
Before filtering:  ~20,000 genes
After filtering:    15,247 genes
Removed:            ~4,753 genes (mostly unexpressed or near-zero)
```

**Why this matters:** Every gene we test is a separate statistical test. Running 20,000 tests instead of 15,000 means more false positives to correct for. Removing uninformative genes makes the analysis cleaner and more powerful.

### Sample metadata

We also build a table that tells the pipeline which sample belongs to which condition and which batch:

```r
sample_info <- data.frame(
  condition = c("Control","Control","Control","Control","Control","Control",
                "PD","PD","PD","PD","PD","PD"),
  batch     = c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2)
)
```

Batch refers to when the samples were processed in the lab. Samples processed at different times can have technical differences unrelated to biology — we track this so we can correct for it later.

---

## Part 2 — Exploratory Data Analysis

### What's happening

Before running any statistics, we want to visually inspect the data. This is a sanity check — do the samples look like they should? Do PD and Control samples separate? Are there any obvious outliers?

### TMM Normalization

Raw counts can't be compared directly between samples because different samples have different total read counts (library sizes). A gene with 100 counts in a sample with 10M total reads is expressed very differently than a gene with 100 counts in a sample with 50M total reads.

**TMM (Trimmed Mean of M-values)** calculates a scaling factor for each sample that accounts for both library size and composition differences.

```
Before normalization:
  Sample A: 10M reads total  →  SNCA = 100 counts
  Sample B: 50M reads total  →  SNCA = 100 counts
  (looks the same, but B is actually 5× lower relative expression)

After TMM normalization:
  Sample A: SNCA = 100 / 0.2  = 500 normalized counts
  Sample B: SNCA = 100 / 1.0  = 100 normalized counts
  (now comparable)
```

### Principal Component Analysis (PCA)

PCA takes the expression of all 15,247 genes across 12 samples and compresses that information into a small number of dimensions (principal components) that capture the most variation.

Think of it like this: if you had 15,247 measurements per person, PCA finds the 2-3 "summary scores" that best explain the differences between people.

```
Full data (15,247 dimensions):
  Sample 1: [45, 312, 104, 89, 201, ...]  ← 15,247 numbers
  Sample 2: [52, 298,  98, 91, 198, ...]
  ...

After PCA (2 dimensions):
  Sample 1: PC1 = -3.2,  PC2 = 0.8   ← Control cluster
  Sample 2: PC1 = -2.9,  PC2 = 1.1   ← Control cluster
  Sample 7: PC1 =  4.1,  PC2 = -0.3  ← PD cluster
  Sample 8: PC1 =  3.8,  PC2 =  0.5  ← PD cluster
```

If the biology is real, PD and Control samples should separate along PC1. If batch is the dominant effect, samples from the same batch will cluster together regardless of condition — that's a red flag.

---

## Part 3 — Differential Expression with DESeq2

### What's happening

This is the core of the analysis. For every one of the 15,247 genes, we ask: **is the expression level in PD samples significantly different from Control?**

DESeq2 uses a statistical model designed specifically for count data. It runs three steps:

```
Step 1: Size factor estimation
        Normalize each sample so they're comparable
        (similar to TMM, but DESeq2's own method)

Step 2: Dispersion estimation
        For each gene, estimate how variable it is across samples
        Genes with high variability need stronger evidence to call significant

Step 3: Wald test
        For each gene, fit a model and test:
        H₀: log2FC = 0  (no difference between PD and Control)
        H₁: log2FC ≠ 0  (there is a difference)
```

### The design formula

```r
design = ~ batch + condition
```

This tells DESeq2: "when testing for condition effects, account for batch first." Batch is a nuisance variable — we want to remove its influence so it doesn't confuse the PD vs. Control comparison.

```
Without batch correction:
  Gene X looks upregulated in PD...
  ...but actually Batch 2 just happened to have higher expression
  and PD samples were over-represented in Batch 2.

With batch correction:
  DESeq2 separates the batch effect from the condition effect
  and only reports the true PD vs. Control difference.
```

### Significance thresholds

```
padj < 0.05        →  less than 5% chance this is a false positive
|log2FC| > 1       →  at least a 2-fold expression difference

Both must be true to call a gene "significant"
```

**Why adjusted p-value?** We're running ~15,000 tests simultaneously. At p < 0.05, we'd expect ~750 false positives by chance. The Benjamini-Hochberg correction controls the false discovery rate — padj < 0.05 means at most 5% of our significant genes are expected to be false positives.

```
Example output for one gene:

  gene_id:         SNCA
  baseMean:        847.3       ← well-expressed gene, reliable estimate
  log2FoldChange:  2.14        ← ~4.4× higher in PD
  pvalue:          3.2e-12     ← extremely unlikely by chance
  padj:            1.1e-09     ← still significant after correction
```

---

## Part 4 — Visualization

### Volcano Plot — reading it

```
High significance
(small padj)
      ▲
      │         ●SNCA
   20 │      ●●●●●●●●●
      │   ●●●●●●●●●●●●●●
      │●●●●●●●●●●●●●●●●●●●●
    0 │─────────────────────── ▶ log2FC
      │●●●●●●●●●●●●●●●●●●●●
      │   ●●●●●●●●●●●●●●
      │      ●●●●●●●●●
      │
      ▼
Low significance

      ◄──────────┬──────────►
           lower in PD  higher in PD
```

- **Top right:** high confidence, upregulated in PD — strongest hits
- **Top left:** high confidence, downregulated in PD
- **Bottom center:** low fold change, not significant — most genes live here

### Heatmap — reading it

```
Z-score: how many standard deviations above/below a gene's own mean

Gene A mean expression = 500 counts, std = 100
  Control sample: 400 counts  →  Z = (400-500)/100 = -1.0  (below average)
  PD sample:      700 counts  →  Z = (700-500)/100 = +2.0  (above average)

This lets you compare genes with very different absolute expression levels
on the same color scale.

  navy = below average for that gene
  white = at the average
  red = above average
```

---

## Part 5 — Pathway Enrichment

### What's happening

A list of 1,834 significant genes is hard to interpret on its own. Pathway enrichment asks: **do these genes cluster into known biological processes more than you'd expect by chance?**

### Gene Ontology (GO) enrichment

GO is a database that groups genes into functional categories (terms) like "mitochondrial electron transport" or "response to oxidative stress."

```
Your upregulated genes: 982 genes

GO term: "response to oxidative stress"
  - Contains 300 genes in the full human genome
  - 47 of your 982 upregulated genes are in this term
  - Expected by chance: ~15 genes
  - That's a 3× enrichment → p-value = 2.1e-08 → significant
```

### KEGG pathway analysis

KEGG maps genes onto curated biological pathways with known interactions — things like "Parkinson's disease pathway" or "oxidative phosphorylation."

```
Top enriched KEGG pathways in downregulated genes:

  Pathway                          Genes hit    p.adjust
  ─────────────────────────────────────────────────────
  Oxidative phosphorylation           38         1.2e-11
  Parkinson's disease                 29         3.4e-09
  Thermogenesis                       31         8.7e-08
  Alzheimer's disease                 27         2.1e-07
```

The overlap with the "Parkinson's disease" KEGG pathway is a positive validation — it means the genes we found are the same ones already known to be involved in PD.

---

## Part 6 — Key Gene Profiles

### Why these five genes

| Gene | Why it matters |
|------|---------------|
| `SNCA` | Encodes α-synuclein, the protein that aggregates into Lewy bodies — the defining pathological feature of PD |
| `LRRK2` | Mutations in LRRK2 are the most common genetic cause of familial PD; its kinase activity is a major drug target |
| `PINK1` | Detects damaged mitochondria and triggers their removal (mitophagy); loss-of-function causes early-onset PD |
| `PRKN` | Works downstream of PINK1 to tag damaged mitochondria for degradation; mutations cause juvenile PD |
| `DJ1` | Protects neurons from oxidative stress; loss-of-function is another cause of early-onset PD |

These five genes represent the major genetic and molecular axes of PD — mitochondrial quality control (PINK1/PRKN), protein aggregation (SNCA), kinase signaling (LRRK2), and oxidative stress (DJ1).

### What the boxplots show

```
  Log2 CPM
     8 │         ┌───┐
       │         │   │  ●
     6 │  ┌───┐  │   │
       │  │   │  └───┘
     4 │  └───┘
       │
       └──────────────
          Control    PD

  The box spans the interquartile range (25th–75th percentile)
  The line inside is the median
  Dots are individual samples
  Wider separation = stronger expression difference
```

---

## Part 7 — Summary & Reproducibility

### Session info

At the end of the run, R saves a full record of the environment to `results/session_info.txt`:

```
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20
Running under: macOS Ventura 13.5

Attached packages:
  DESeq2_1.40.2
  edgeR_3.42.4
  ggplot2_3.4.3
  ...
```

This means anyone can reproduce the exact analysis by matching these package versions — important for peer review and long-term reproducibility.

---

## Glossary

| Term | Plain-English definition |
|------|--------------------------|
| **RNA-seq** | Technique that sequences all RNA in a cell to measure gene expression |
| **Count matrix** | Table of raw sequencing read counts per gene per sample |
| **CPM** | Counts Per Million — normalizes for library size |
| **TMM** | Normalization method that accounts for library size and composition |
| **DEG** | Differentially Expressed Gene — significantly different between conditions |
| **log2FC** | log2 Fold Change — log2(PD/Control); 1 = 2× higher, -1 = 2× lower |
| **padj** | Adjusted p-value — corrected for multiple testing across all genes |
| **FDR** | False Discovery Rate — expected fraction of false positives in significant hits |
| **PCA** | Principal Component Analysis — dimensionality reduction for visualization |
| **GO** | Gene Ontology — database of functional gene categories |
| **KEGG** | Database of curated biological pathways |
| **GSEA** | Gene Set Enrichment Analysis — tests whether gene sets are enriched |
| **iPSC** | Induced Pluripotent Stem Cell — adult cells reprogrammed to stem cell state |
| **Mitophagy** | Cellular process that removes damaged mitochondria |
| **Lewy body** | Protein aggregate (mainly α-synuclein) found in PD neurons |
| **Batch effect** | Technical variation between samples processed at different times |
