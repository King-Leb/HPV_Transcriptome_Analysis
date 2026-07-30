############################################################
# Load required packages
############################################################
library(dplyr)
library(purrr)
library(DESeq2)
library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(apeglm)
library(VennDiagram)
library(EnhancedVolcano)
library(WGCNA)
library(matrixStats)
library(biomaRt)
library(clusterProfiler)
library(igraph)
library(ggraph)
library(ggplot2)
library(dynamicTreeCut)
library(knitr)
library(enrichplot)
library(tidyverse)
library(tidygraph)
library(RCy3)
library(CorLevelPlot)
library(reshape2)
library(STRINGdb)
library(RColorBrewer)
library(pheatmap)
library(vioplot)
library(writexl)
library(corrplot)
library(Cairo)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)
library(patchwork)
library(svglite)
library(grid)  # for unit()
library(gridGraphics)
library(ggplotify)

############################################################
# Import count matrices and metadata
############################################################
vscc_meta <- read.delim("vscc_metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
opscc_meta <- read.delim("opscc_metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
hnc_meta <- read.delim("hnc_metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
colnames(vscc_meta) <- colnames(opscc_meta) <- colnames(hnc_meta)

# Keep only count columns (drop Length if present)
# Read the count matrices
vscc_counts <- read.delim("vscc_counts_matrix.txt", stringsAsFactors = FALSE, check.names = FALSE)
opscc_counts <- read.delim("opscc_counts_matrix.txt", stringsAsFactors = FALSE, check.names = FALSE)
hnc_counts   <- read.delim("hnc_counts_matrix.txt", stringsAsFactors = FALSE, check.names = FALSE)

# Keep gene IDs and counts
vscc_counts <- vscc_counts[, c("Geneid", grep("^SRR", colnames(vscc_counts), value = TRUE))]
opscc_counts <- opscc_counts[, c("Geneid", grep("^SRR", colnames(opscc_counts), value = TRUE))]
hnc_counts   <- hnc_counts[,   c("Geneid", grep("^SRR", colnames(hnc_counts), value = TRUE))]

vscc_counts[is.na(vscc_counts)] <- 0
opscc_counts[is.na(opscc_counts)] <- 0
hnc_counts[is.na(hnc_counts)] <- 0

vscc_counts$Geneid <- sub("\\..*$", "", vscc_counts$Geneid)
opscc_counts$Geneid <- sub("\\..*$", "", opscc_counts$Geneid)
hnc_counts$Geneid <- sub("\\..*$", "", hnc_counts$Geneid)

# Make Geneid the row names
rownames(vscc_counts) <- vscc_counts$Geneid
rownames(opscc_counts) <- opscc_counts$Geneid
rownames(hnc_counts) <- hnc_counts$Geneid

# Drop the Geneid column so it doesn't interfere
vscc_counts <- vscc_counts[ , -1]
opscc_counts <- opscc_counts[ , -1]
hnc_counts <- hnc_counts[ , -1]

## should return TRUE
all(colnames(vscc_counts) == vscc_meta$SRA)
all(colnames(opscc_counts) == opscc_meta$SRA)
all(colnames(hnc_counts) == hnc_meta$SRA)

############################################################
# 3. Differential expression analysis (DESeq2: design = ~ HPV_Status [+ Gender])
############################################################
# --------------------------------------
# HNC (PCA + DESeq)
# --------------------------------------
dds_hnc <- DESeqDataSetFromMatrix(
  countData = hnc_counts,
  colData = hnc_meta,
  design = ~ HPV_Status + Gender
)
vsd_hnc <- vst(dds_hnc, blind = TRUE)
pca_data_hnc <- plotPCA(vsd_hnc, intgroup = c("Gender", "HPV_Status"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data_hnc, "percentVar"))
pca_plot_batch_hnc <- ggplot(pca_data_hnc, aes(PC1, PC2, color = Gender, shape = HPV_Status)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA of Variance-Stabilized Counts of HNC")
ggsave(filename = "PCA_plot_batch_hnc.png", plot = pca_plot_batch_hnc, dpi = 900, width = 7, height = 5, units = "in")

keep_hnc_deseq <- rowSums(counts(dds_hnc) >= 10) >= 10
dds_hnc <- dds_hnc[keep_hnc_deseq, ]
dds_hnc$HPV_Status <- relevel(dds_hnc$HPV_Status, ref = 'neg')
deseq_hnc <- DESeq(dds_hnc)
results_hnc <- results(deseq_hnc, contrast = c('HPV_Status', 'pos', 'neg'))
results_hnc <- as.data.frame(results_hnc)
results_hnc <- results_hnc[order(results_hnc$pvalue), ]
df_hnc <- as.data.frame(-log10(results_hnc$pvalue) * sign(results_hnc$log2FoldChange),
                        row.names = rownames(results_hnc))
colnames(df_hnc) <- "HNC"

# --------------------------------------
# OPSCC (PCA + DESeq)
# --------------------------------------
dds_opscc <- DESeqDataSetFromMatrix(
  countData = opscc_counts,
  colData = opscc_meta,
  design = ~ HPV_Status + Gender
)
vsd_opscc <- vst(dds_opscc, blind = TRUE)
pca_data_opscc <- plotPCA(vsd_opscc, intgroup = c("Gender", "HPV_Status"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data_opscc, "percentVar"))
pca_plot_batch_opscc <- ggplot(pca_data_opscc, aes(PC1, PC2, color = Gender, shape = HPV_Status)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA of Variance-Stabilized Counts of OPSCC")
ggsave(filename = "PCA_plot_batch_opscc.png", plot = pca_plot_batch_opscc, dpi = 900, width = 7, height = 5, units = "in")

keep_opscc_deseq <- rowSums(counts(dds_opscc) >= 10) >= 10
dds_opscc <- dds_opscc[keep_opscc_deseq, ]
dds_opscc$HPV_Status <- relevel(dds_opscc$HPV_Status, ref = 'neg')
deseq_opscc <- DESeq(dds_opscc)
results_opscc <- results(deseq_opscc, contrast = c('HPV_Status', 'pos', 'neg'))
results_opscc <- as.data.frame(results_opscc)
results_opscc <- results_opscc[order(results_opscc$pvalue), ]
df_opscc <- as.data.frame(-log10(results_opscc$pvalue) * sign(results_opscc$log2FoldChange),
                          row.names = rownames(results_opscc))
colnames(df_opscc) <- "OPSCC"

# --------------------------------------
# VSCC (PCA + DESeq)
# --------------------------------------
dds_vscc <- DESeqDataSetFromMatrix(
  countData = vscc_counts,
  colData = vscc_meta,
  design = ~ HPV_Status
)
vsd_vscc <- vst(dds_vscc, blind = TRUE)
pca_data_vscc <- plotPCA(vsd_vscc, intgroup = c("HPV_Status"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data_vscc, "percentVar"))
pca_plot_batch_vscc <- ggplot(pca_data_vscc, aes(PC1, PC2, shape = HPV_Status)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA of Variance-Stabilized Counts of VSCC")
ggsave(filename = "PCA_plot_batch_vscc.png", plot = pca_plot_batch_vscc, dpi = 900, width = 7, height = 5, units = "in")

keep_vscc_deseq <- rowSums(counts(dds_vscc) >= 10) >= 10
dds_vscc <- dds_vscc[keep_vscc_deseq, ]
dds_vscc$HPV_Status <- relevel(dds_vscc$HPV_Status, ref = 'neg')
deseq_vscc <- DESeq(dds_vscc)
results_vscc <- results(deseq_vscc, contrast = c('HPV_Status', 'pos', 'neg'))
results_vscc <- as.data.frame(results_vscc)
results_vscc <- results_vscc[order(results_vscc$pvalue), ]
df_vscc <- as.data.frame(-log10(results_vscc$pvalue) * sign(results_vscc$log2FoldChange),
                         row.names = rownames(results_vscc))
colnames(df_vscc) <- "VSCC"

# --------------------------------------------------
# Merge the dataframes into a signature matrix
# --------------------------------------------------
dflist <- list(df_hnc, df_opscc, df_vscc)
dflist_named <- lapply(dflist, rownames_to_column, var = "gene")
signatureMatrix <- purrr::reduce(
  dflist_named,
  dplyr::full_join,
  by = "gene"
) %>%
  tibble::column_to_rownames("gene")

# ----------------------------------------
# Add gene symbols
# ----------------------------------------
ensembl_ids <- rownames(signatureMatrix)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
# Add as a new column
signatureMatrix <- cbind(GeneSymbol = gene_symbols[rownames(signatureMatrix)],
                         signatureMatrix)
save(signatureMatrix, file = "signatureMatrix.rda")
save(df_hnc, df_opscc, df_vscc, signatureMatrix, file = "merged_data.rda")
write_xlsx(signatureMatrix, "signatureMatrix.xlsx")
# load(file = "signatureMatrix.rda")
signatureMatrix <- signatureMatrix[!is.na(signatureMatrix$GeneSymbol), ]
rownames(signatureMatrix) <- make.unique(as.character(signatureMatrix$GeneSymbol))
sig_matrix <- signatureMatrix[, -which(names(signatureMatrix) == "GeneSymbol")]

# ----------------------------------------------
# Apply LFC shrinkage using apeglm
# ----------------------------------------------
res_hnc_shrunk <- lfcShrink(deseq_hnc,
                            coef = "HPV_Status_pos_vs_neg",
                            type = "apeglm")
res_opscc_shrunk <- lfcShrink(deseq_opscc,
                              coef = "HPV_Status_pos_vs_neg",
                              type = "apeglm")
res_vscc_shrunk <- lfcShrink(deseq_vscc,
                             coef = "HPV_Status_pos_vs_neg",
                             type = "apeglm")

# ----------------------------------------------
# Convert to data frames and add gene symbols using biomaRt
# ----------------------------------------------
add_gene_symbols <- function(df, mart) {
  ensembl_ids <- rownames(df)
  annot <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = mart
  )
  df$ensembl_gene_id <- ensembl_ids
  df_annot <- merge(df, annot, by = "ensembl_gene_id", all.x = TRUE)
  return(df_annot)
}

# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

prepare_results <- function(res, mart) {
  df <- as.data.frame(res)
  ensembl_ids <- rownames(df)
  annot <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = mart
  )
  df$ensembl_gene_id <- ensembl_ids
  df_annot <- merge(df, annot, by = "ensembl_gene_id", all.x = TRUE)
  return(df_annot)
}

results_hnc_shrunk   <- prepare_results(res_hnc_shrunk, mart)
results_opscc_shrunk <- prepare_results(res_opscc_shrunk, mart)
results_vscc_shrunk  <- prepare_results(res_vscc_shrunk, mart)

# HNC
degs_hnc <- results_hnc_shrunk[!is.na(results_hnc_shrunk$padj) & results_hnc_shrunk$padj < 0.05, ]
up_hnc   <- nrow(degs_hnc[degs_hnc$log2FoldChange > 0, ])
down_hnc <- nrow(degs_hnc[degs_hnc$log2FoldChange < 0, ])
# OPSCC
degs_opscc <- results_opscc_shrunk[!is.na(results_opscc_shrunk$padj) & results_opscc_shrunk$padj < 0.05, ]
up_opscc   <- nrow(degs_opscc[degs_opscc$log2FoldChange > 0, ])
down_opscc <- nrow(degs_opscc[degs_opscc$log2FoldChange < 0, ])
# VSCC
degs_vscc <- results_vscc_shrunk[!is.na(results_vscc_shrunk$padj) & results_vscc_shrunk$padj < 0.05, ]
up_vscc   <- nrow(degs_vscc[degs_vscc$log2FoldChange > 0, ])
down_vscc <- nrow(degs_vscc[degs_vscc$log2FoldChange < 0, ])

cat("HNC  - Up:", up_hnc,  "| Down:", down_hnc,  "\n")
cat("OPSCC- Up:", up_opscc,"| Down:", down_opscc, "\n")
cat("VSCC - Up:", up_vscc, "| Down:", down_vscc,  "\n")
