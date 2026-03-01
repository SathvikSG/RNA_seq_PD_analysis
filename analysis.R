# RNA-seq Differential Expression Analysis
# Analyzing gene expression changes in neurodegenerative disease models
# Author: Sathvik Sai Guntha
# Dataset: GSE157827 (Parkinson's Disease iPSC-derived neurons vs. controls)

# Load required libraries
library(DESeq2)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(tidyverse)
library(RColorBrewer)
library(gplots)

# Set working directory and create output folders
setwd("~/RNA_seq_analysis")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# ============================================================================
# PART 1: Data Import and Preprocessing
# ============================================================================

# Read in raw count matrix (rows = genes, columns = samples)
counts_raw <- read.csv("data/raw_counts.csv", row.names = 1)

# Create sample metadata
sample_info <- data.frame(
  sample_id = colnames(counts_raw),
  condition = factor(c(rep("Control", 6), rep("PD", 6)),
                     levels = c("Control", "PD")),
  batch = factor(c(rep(1, 3), rep(2, 3), rep(1, 3), rep(2, 3))),
  row.names = colnames(counts_raw)
)

# Quality control: Filter low-count genes
# Keep genes with at least 10 counts in at least 3 samples
keep <- rowSums(counts_raw >= 10) >= 3
counts_filtered <- counts_raw[keep, ]

cat("Initial genes:", nrow(counts_raw), "\n")
cat("After filtering:", nrow(counts_filtered), "\n")
cat("Genes removed:", nrow(counts_raw) - nrow(counts_filtered), "\n")

# ============================================================================
# PART 2: Exploratory Data Analysis
# ============================================================================

# Calculate TPM for visualization (not for DE analysis)
calculate_tpm <- function(counts, gene_lengths) {
  rate <- counts / gene_lengths
  tpm <- rate / sum(rate) * 1e6
  return(tpm)
}

# Normalize using TMM (trimmed mean of M-values)
dge <- DGEList(counts = counts_filtered, group = sample_info$condition)
dge <- calcNormFactors(dge, method = "TMM")

# Log-transform for visualization
logcpm <- cpm(dge, log = TRUE, prior.count = 2)

# PCA analysis
pca_result <- prcomp(t(logcpm), scale. = TRUE)
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  condition = sample_info$condition,
  batch = sample_info$batch
)

# Plot PCA
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of RNA-seq Samples",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")) +
  scale_color_manual(values = c("Control" = "#3498db", "PD" = "#e74c3c")) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/pca_plot.png", pca_plot, width = 8, height = 6, dpi = 300)

# Sample-to-sample distance heatmap
sample_dists <- dist(t(logcpm))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste(sample_info$condition,
                                       sample_info$batch, sep = "_")
colnames(sample_dist_matrix) <- rownames(sample_dist_matrix)

png("figures/sample_distances.png", width = 800, height = 800, res = 120)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-Sample Distances")
dev.off()

# ============================================================================
# PART 3: Differential Expression Analysis with DESeq2
# ============================================================================

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_info,
  design = ~ batch + condition
)

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "PD", "Control"))
res_ordered <- res[order(res$padj), ]

# Summary of results
summary(res)

# Save results table
res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj))

write.csv(res_df, "results/deseq2_results_full.csv", row.names = FALSE)

# Extract significant genes (padj < 0.05, |log2FC| > 1)
sig_genes <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

write.csv(sig_genes, "results/significant_genes.csv", row.names = FALSE)

cat("\nDifferential Expression Summary:\n")
cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant genes (padj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")
cat("Upregulated in PD:", sum(sig_genes$log2FoldChange > 1), "\n")
cat("Downregulated in PD:", sum(sig_genes$log2FoldChange < -1), "\n")

# ============================================================================
# PART 4: Visualization of Results
# ============================================================================

# MA plot
png("figures/ma_plot.png", width = 1000, height = 800, res = 120)
plotMA(res, ylim = c(-5, 5), main = "MA Plot: PD vs Control")
dev.off()

# Volcano plot
volcano_plot <- EnhancedVolcano(res_df,
                                lab = res_df$gene_id,
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Differential Expression: PD vs Control',
                                pCutoff = 0.05,
                                FCcutoff = 1,
                                pointSize = 2.0,
                                labSize = 4.0,
                                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                                colAlpha = 0.7,
                                legendPosition = 'right',
                                legendLabSize = 12,
                                legendIconSize = 4.0)

ggsave("figures/volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)

# Heatmap of top 50 differentially expressed genes
top_genes <- head(sig_genes$gene_id, 50)
mat <- logcpm[top_genes, ]

# Z-score normalization for heatmap
mat_scaled <- t(scale(t(mat)))

annotation_col <- data.frame(
  Condition = sample_info$condition,
  row.names = colnames(mat)
)

png("figures/heatmap_top50.png", width = 1000, height = 1200, res = 120)
pheatmap(mat_scaled,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Top 50 Differentially Expressed Genes",
         fontsize_row = 6)
dev.off()

# ============================================================================
# PART 5: Gene Set Enrichment Analysis
# ============================================================================

# Prepare gene list for GSEA
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$gene_id
gene_list <- sort(gene_list, decreasing = TRUE)

# Get gene symbols from Ensembl IDs (if needed)
# Assuming gene_id format is Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# For this example, we'll use the gene IDs directly
# In practice, you would convert Ensembl IDs to Entrez IDs

# GO enrichment analysis for upregulated genes
upreg_genes <- sig_genes %>%
  filter(log2FoldChange > 1) %>%
  pull(gene_id)

downreg_genes <- sig_genes %>%
  filter(log2FoldChange < -1) %>%
  pull(gene_id)

# Note: This requires gene symbols or Entrez IDs
# For demonstration, assuming conversion has been done
# ego_up <- enrichGO(gene = upreg_genes,
#                    OrgDb = org.Hs.eg.db,
#                    keyType = "SYMBOL",
#                    ont = "BP",
#                    pAdjustMethod = "BH",
#                    pvalueCutoff = 0.05,
#                    qvalueCutoff = 0.2)

# KEGG pathway analysis
# kegg_up <- enrichKEGG(gene = upreg_genes_entrez,
#                       organism = 'hsa',
#                       pvalueCutoff = 0.05)

# Save enrichment results
# write.csv(as.data.frame(ego_up), "results/GO_enrichment_upregulated.csv", row.names = FALSE)

# ============================================================================
# PART 6: Expression Profiles of Key Genes
# ============================================================================

# Plot expression of specific genes of interest
genes_of_interest <- c("SNCA", "LRRK2", "PINK1", "PRKN", "DJ1")  # PD-related genes

# Filter for genes present in dataset
goi_present <- genes_of_interest[genes_of_interest %in% rownames(logcpm)]

if (length(goi_present) > 0) {
  goi_expr <- logcpm[goi_present, , drop = FALSE]
  goi_df <- goi_expr %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    left_join(sample_info %>% rownames_to_column("sample"), by = "sample")

  goi_plot <- ggplot(goi_df, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    facet_wrap(~ gene, scales = "free_y", ncol = 3) +
    theme_minimal() +
    labs(title = "Expression of PD-Related Genes",
         y = "Log2 CPM",
         x = "Condition") +
    scale_fill_manual(values = c("Control" = "#3498db", "PD" = "#e74c3c")) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.text = element_text(face = "bold"))

  ggsave("figures/key_genes_expression.png", goi_plot, width = 12, height = 8, dpi = 300)
}

# ============================================================================
# PART 7: Generate Summary Report
# ============================================================================

# Create a summary statistics file
summary_stats <- list(
  total_samples = ncol(counts_filtered),
  control_samples = sum(sample_info$condition == "Control"),
  pd_samples = sum(sample_info$condition == "PD"),
  genes_analyzed = nrow(counts_filtered),
  significant_genes = nrow(sig_genes),
  upregulated = sum(sig_genes$log2FoldChange > 1),
  downregulated = sum(sig_genes$log2FoldChange < -1),
  threshold_padj = 0.05,
  threshold_log2fc = 1
)

writeLines(
  c("RNA-seq Differential Expression Analysis Summary",
    paste(rep("=", 50), collapse = ""),
    paste("Analysis Date:", Sys.Date()),
    paste("Total Samples:", summary_stats$total_samples),
    paste("  - Control:", summary_stats$control_samples),
    paste("  - Parkinson's Disease:", summary_stats$pd_samples),
    "",
    paste("Genes Analyzed:", summary_stats$genes_analyzed),
    paste("Significant Genes (padj < 0.05, |log2FC| > 1):", summary_stats$significant_genes),
    paste("  - Upregulated in PD:", summary_stats$upregulated),
    paste("  - Downregulated in PD:", summary_stats$downregulated),
    "",
    "Output files generated:",
    "  - results/deseq2_results_full.csv",
    "  - results/significant_genes.csv",
    "  - figures/pca_plot.png",
    "  - figures/volcano_plot.png",
    "  - figures/heatmap_top50.png"
  ),
  "results/analysis_summary.txt"
)

cat("\nAnalysis complete! Check the 'results' and 'figures' directories for outputs.\n")

# Save R session info for reproducibility
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
