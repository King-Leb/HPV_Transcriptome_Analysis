############################################################
#    Weighted Gene Co-expression Network Analysis (WGCNA)
#    NOTE: design = ~1 here (unsupervised network construction)
############################################################
options(stringsAsFactors = FALSE)
allowWGCNAThreads()  # parallelization

# ==============================================
# HNC
# ==============================================
dds_hnc_wg <- DESeqDataSetFromMatrix(
  countData = hnc_counts,
  colData = hnc_meta,
  design = ~1
)
dds_hnc_wg <- DESeq(dds_hnc_wg)
vsd_hnc_wg <- vst(dds_hnc_wg, blind = TRUE)
norm_counts_hnc <- assay(vsd_hnc_wg)
geneVars_hnc <- rowVars(norm_counts_hnc)
topGenes_hnc <- order(geneVars_hnc, decreasing = TRUE)[1:7000]  # keep top 7000 variable genes
exprMatrix_hnc <- norm_counts_hnc[topGenes_hnc, ]
datExpr_hnc <- t(exprMatrix_hnc)
gsg_hnc <- goodSamplesGenes(datExpr_hnc, verbose = 3)
if (!gsg_hnc$allOK){
  datExpr_hnc <- datExpr_hnc[gsg_hnc$goodSamples, gsg_hnc$goodGenes]
}
datExprMat_hnc <- as.matrix(datExpr_hnc)
cor_samples_hnc <- bicor(t(datExprMat_hnc), use = "pairwise.complete.obs")
adj_samples_hnc <- (1 + cor_samples_hnc) / 2
K_hnc <- rowSums(adj_samples_hnc) - diag(adj_samples_hnc)
zK_hnc <- scale(K_hnc)
zK_hnc <- as.numeric(zK_hnc)
names(zK_hnc) <- rownames(datExprMat_hnc)
outlier_threshold_hnc <- -2
outliers_conn_hnc <- names(zK_hnc)[which(zK_hnc < outlier_threshold_hnc)]
message("Outliers detected by connectivity in hnc (Z.K < ", outlier_threshold_hnc, "): ",
        paste(outliers_conn_hnc, collapse = ", "))
keep_samples_conn_hnc <- !(rownames(datExprMat_hnc) %in% outliers_conn_hnc)
datExpr_filtered_conn_hnc <- datExprMat_hnc[keep_samples_conn_hnc, , drop = FALSE]
metadata_filtered_conn_hnc <- hnc_meta[hnc_meta$SRA %in% rownames(datExpr_filtered_conn_hnc), , drop = FALSE]
sampleTree2_hnc <- hclust(dist(datExpr_filtered_conn_hnc), method = "average")
png("sample_clustering_filtered_hnc.png", width = 14, height = 5, units = "in", res = 900)
plot(sampleTree2_hnc, main = "Sample clustering after outlier removal", cex = 1.0)
dev.off()
traitData_hnc <- metadata_filtered_conn_hnc
Status_Colors_hnc <- setNames(c("grey", "red"), levels(factor(traitData_hnc$HPV_Status)))
traitData_hnc$status_col <- Status_Colors_hnc[traitData_hnc$HPV_Status]
traitColors_hnc <- data.frame(
  Status = traitData_hnc$status_col
)
png("sample_dendrogram_hnc.png", width = 12, height = 5, units = "in", res = 900)
par(mar = c(2, 4, 4, 2), xpd = TRUE)
plotDendroAndColors(sampleTree2_hnc,
                    colors = traitColors_hnc,
                    groupLabels = c("HPV Status"),
                    main = "Sample dendrogram with trait")
legend("bottomleft", legend = names(Status_Colors_hnc), fill = Status_Colors_hnc, title = "HPV Status", cex = 0.7, bty = "n")
dev.off()
softPower_hnc <- 12
netSingleBlock_hnc <- blockwiseModules(
  datExpr_filtered_conn_hnc,
  power = softPower_hnc,
  TOMType = "signed",
  minModuleSize = 40,
  reassignThreshold = 0,
  mergeCutHeight = 0.3,     # post-merge
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(datExpr_filtered_conn_hnc),# ensures all genes in one block
  deepSplit = 1,
  verbose = 3
)
mergedColors_hnc <- labels2colors(netSingleBlock_hnc$colors)
netSingleBlockUnmerged_hnc <- blockwiseModules(
  datExpr_filtered_conn_hnc,
  power = softPower_hnc,
  TOMType = "signed",
  minModuleSize = 40,
  reassignThreshold = 0,
  mergeCutHeight = 0,        # disables merging
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(datExpr_filtered_conn_hnc),
  verbose = 3
)
unmergedColors_hnc <- labels2colors(netSingleBlockUnmerged_hnc$colors)
png("Gene_dendrogram_all_genes_pre_vs_post_merge_hnc.png",
    width = 16, height = 8, units = "in", res = 900)
par(cex = 0.8, mar = c(5,5,5,5))
plotDendroAndColors(
  dendro = netSingleBlock_hnc$dendrograms[[1]],
  colors = cbind(unmergedColors_hnc, mergedColors_hnc),
  groupLabels = c("Dynamic Tree Cut", "Merged Dynamic"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "HNC Cluster Dendrogram"
)
dev.off()
nSamples_hnc <- nrow(datExpr_filtered_conn_hnc)
MEs0_hnc <- moduleEigengenes(datExpr_filtered_conn_hnc, mergedColors_hnc)$eigengenes
MEs_hnc <- orderMEs(MEs0_hnc)
traitData_num_hnc <- data.frame(
  HPV_Status = as.numeric(factor(traitData_hnc$HPV_Status)),
  Gender = as.numeric(factor(traitData_hnc$Gender))
)
moduleTraitCor_hnc <- cor(MEs_hnc, traitData_num_hnc, use = "p")
moduleTraitPvalue_hnc <- corPvalueStudent(moduleTraitCor_hnc, nSamples_hnc)
module_names_clean_hnc <- sub("^ME", "", names(MEs_hnc))
textMatrix_hnc <- paste0(signif(moduleTraitCor_hnc, 2), "\n(", signif(moduleTraitPvalue_hnc, 1), ")")
png("module_trait_relationships_hnc.png", width = 9, height = 7, units = "in", res = 900)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor_hnc,
  xLabels = names(traitData_num_hnc),
  yLabels = module_names_clean_hnc,
  ySymbols = module_names_clean_hnc,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_hnc,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "HNC Module - Trait Relationships"
)
dev.off()

# ==============================================
# OPSCC
# ==============================================
dds_opscc_wg <- DESeqDataSetFromMatrix(
  countData = opscc_counts,
  colData = opscc_meta,
  design = ~1
)
dds_opscc_wg <- DESeq(dds_opscc_wg)
vsd_opscc_wg <- vst(dds_opscc_wg, blind = TRUE)
norm_counts_opscc <- assay(vsd_opscc_wg)
geneVars_opscc <- rowVars(norm_counts_opscc)
topGenes_opscc <- order(geneVars_opscc, decreasing = TRUE)[1:7000]
exprMatrix_opscc <- norm_counts_opscc[topGenes_opscc, ]
datExpr_opscc <- t(exprMatrix_opscc)
gsg_opscc <- goodSamplesGenes(datExpr_opscc, verbose = 3)
if (!gsg_opscc$allOK){
  datExpr_opscc <- datExpr_opscc[gsg_opscc$goodSamples, gsg_opscc$goodGenes]
}
datExprMat_opscc <- as.matrix(datExpr_opscc)
cor_samples_opscc <- bicor(t(datExprMat_opscc), use = "pairwise.complete.obs")
adj_samples_opscc <- (1 + cor_samples_opscc) / 2
K_opscc <- rowSums(adj_samples_opscc) - diag(adj_samples_opscc)
zK_opscc <- scale(K_opscc)
zK_opscc <- as.numeric(zK_opscc)
names(zK_opscc) <- rownames(datExprMat_opscc)
outliers_conn_opscc <- names(zK_opscc)[which(zK_opscc < -2)]
keep_samples_conn_opscc <- !(rownames(datExprMat_opscc) %in% outliers_conn_opscc)
datExpr_filtered_conn_opscc <- datExprMat_opscc[keep_samples_conn_opscc, , drop = FALSE]
metadata_filtered_conn_opscc <- opscc_meta[opscc_meta$SRA %in% rownames(datExpr_filtered_conn_opscc), , drop = FALSE]
sampleTree2_opscc <- hclust(dist(datExpr_filtered_conn_opscc), method = "average")
png("sample_dendrogram_opscc.png", width = 12, height = 5, units = "in", res = 900)
par(mar = c(2, 4, 4, 2), xpd = TRUE)
plot(sampleTree2_opscc, main = "OPSCC Sample dendrogram")
dev.off()
traitData_opscc <- metadata_filtered_conn_opscc
Status_Colors_opscc <- setNames(c("grey", "red"), levels(factor(traitData_opscc$HPV_Status)))
traitData_opscc$status_col <- Status_Colors_opscc[traitData_opscc$HPV_Status]
traitColors_opscc <- data.frame(Status = traitData_opscc$status_col)
png("sample_dendrogram_trait_opscc.png", width = 12, height = 5, units = "in", res = 900)
plotDendroAndColors(sampleTree2_opscc, traitColors_opscc,
                    groupLabels = "HPV Status",
                    main = "OPSCC dendrogram with trait")
legend("bottomleft", legend = names(Status_Colors_opscc),
       fill = Status_Colors_opscc, cex = 0.7, bty = "n")
dev.off()
softPower_opscc <- 11
netSingleBlock_opscc <- blockwiseModules(
  datExpr_filtered_conn_opscc,
  power = softPower_opscc,
  TOMType = "signed",
  minModuleSize = 40,
  reassignThreshold = 0,
  mergeCutHeight = 0.3,     # post-merge
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(datExpr_filtered_conn_opscc),# ensures all genes in one block
  deepSplit = 2,
  verbose = 3
)
mergedColors_opscc <- labels2colors(netSingleBlock_opscc$colors)
netSingleBlockUnmerged_opscc <- blockwiseModules(
  datExpr_filtered_conn_opscc,
  power = softPower_opscc,
  TOMType = "signed",
  minModuleSize = 40,
  reassignThreshold = 0,
  mergeCutHeight = 0,        # disables merging
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(datExpr_filtered_conn_opscc),
  verbose = 3
)
unmergedColors_opscc <- labels2colors(netSingleBlockUnmerged_opscc$colors)
png("Gene_dendrogram_all_genes_pre_vs_post_merge_opscc.png",
    width = 16, height = 8, units = "in", res = 900)
par(cex = 0.8, mar = c(5,5,5,5))
plotDendroAndColors(
  dendro = netSingleBlock_opscc$dendrograms[[1]],
  colors = cbind(unmergedColors_opscc, mergedColors_opscc),
  groupLabels = c("Dynamic Tree Cut", "Merged Dynamic"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "OPSCC Cluster Dendrogram"
)
dev.off()
table(mergedColors_opscc)
nSamples_opscc <- nrow(datExpr_filtered_conn_opscc)
MEs0_opscc <- moduleEigengenes(datExpr_filtered_conn_opscc, mergedColors_opscc)$eigengenes
MEs_opscc <- orderMEs(MEs0_opscc)
traitData_num_opscc <- data.frame(
  HPV_Status = as.numeric(factor(traitData_opscc$HPV_Status)),
  Gender = as.numeric(factor(traitData_opscc$Gender))
)
moduleTraitCor_opscc <- cor(MEs_opscc, traitData_num_opscc, use = "p")
moduleTraitPvalue_opscc <- corPvalueStudent(moduleTraitCor_opscc, nSamples_opscc)
module_names_clean_opscc <- sub("^ME", "", names(MEs_opscc))
textMatrix_opscc <- paste0(signif(moduleTraitCor_opscc, 2), "\n(", signif(moduleTraitPvalue_opscc, 1), ")")
png("module_trait_relationships_opscc.png", width = 9, height = 7, units = "in", res = 900)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor_opscc,
  xLabels = names(traitData_num_opscc),
  yLabels = module_names_clean_opscc,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_opscc,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "OPSCC Module - Trait Relationships"
)
dev.off()

# ==============================================
# VSCC
# ==============================================
dds_vscc_wg <- DESeqDataSetFromMatrix(
  countData = vscc_counts,
  colData = vscc_meta,
  design = ~1
)
dds_vscc_wg <- DESeq(dds_vscc_wg)
vsd_vscc_wg <- vst(dds_vscc_wg, blind = TRUE)
norm_counts_vscc <- assay(vsd_vscc_wg)
geneVars_vscc <- rowVars(norm_counts_vscc)
topGenes_vscc <- order(geneVars_vscc, decreasing = TRUE)[1:7000]
exprMatrix_vscc <- norm_counts_vscc[topGenes_vscc, ]
datExpr_vscc <- t(exprMatrix_vscc)
gsg_vscc <- goodSamplesGenes(datExpr_vscc, verbose = 3)
if (!gsg_vscc$allOK){
  datExpr_vscc <- datExpr_vscc[gsg_vscc$goodSamples, gsg_vscc$goodGenes]
}
datExprMat_vscc <- as.matrix(datExpr_vscc)
cor_samples_vscc <- bicor(t(datExprMat_vscc), use = "pairwise.complete.obs")
adj_samples_vscc <- (1 + cor_samples_vscc) / 2
K_vscc <- rowSums(adj_samples_vscc) - diag(adj_samples_vscc)
zK_vscc <- as.numeric(scale(K_vscc))
names(zK_vscc) <- rownames(datExprMat_vscc)
outliers_conn_vscc <- names(zK_vscc)[which(zK_vscc < -2)]
keep_samples_conn_vscc <- !(rownames(datExprMat_vscc) %in% outliers_conn_vscc)
datExpr_filtered_conn_vscc <- datExprMat_vscc[keep_samples_conn_vscc, , drop = FALSE]
metadata_filtered_conn_vscc <- vscc_meta[vscc_meta$SRA %in% rownames(datExpr_filtered_conn_vscc), , drop = FALSE]
sampleTree2_vscc <- hclust(dist(datExpr_filtered_conn_vscc), method = "average")
png("sample_dendrogram_vscc.png", width = 12, height = 5, units = "in", res = 900)
plot(sampleTree2_vscc, main = "VSCC Sample dendrogram")
dev.off()
traitData_vscc <- metadata_filtered_conn_vscc
Status_Colors_vscc <- setNames(c("grey","red"), levels(factor(traitData_vscc$HPV_Status)))
traitData_vscc$status_col <- Status_Colors_vscc[traitData_vscc$HPV_Status]
traitColors_vscc <- data.frame(Status = traitData_vscc$status_col)
png("sample_dendrogram_trait_vscc.png", width = 12, height = 5, units = "in", res = 900)
plotDendroAndColors(sampleTree2_vscc, traitColors_vscc,
                    groupLabels = "HPV Status",
                    main = "VSCC dendrogram with trait")
legend("bottomleft", legend = names(Status_Colors_vscc), fill = Status_Colors_vscc, title = "HPV Status", cex = 0.7, bty = "n")
dev.off()
softPower_vscc <- 18
netSingleBlock_vscc <- blockwiseModules(
  datExpr_filtered_conn_vscc,
  power = softPower_vscc,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,     # post-merge
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(datExpr_filtered_conn_vscc),# ensures all genes in one block
  deepSplit = 1,
  verbose = 3
)
mergedColors_vscc <- labels2colors(netSingleBlock_vscc$colors)
netSingleBlockUnmerged_vscc <- blockwiseModules(
  datExpr_filtered_conn_vscc,
  power = softPower_vscc,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0,        # disables merging
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(datExpr_filtered_conn_vscc),
  verbose = 3
)
unmergedColors_vscc <- labels2colors(netSingleBlockUnmerged_vscc$colors)
png("Gene_dendrogram_all_genes_pre_vs_post_merge_vscc.png",
    width = 16, height = 8, units = "in", res = 900)
par(cex = 0.8, mar = c(5,5,5,5))
plotDendroAndColors(
  dendro = netSingleBlock_vscc$dendrograms[[1]],
  colors = cbind(unmergedColors_vscc, mergedColors_vscc),
  groupLabels = c("Dynamic Tree Cut", "Merged Dynamic"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "VSCC Cluster Dendrogram"
)
dev.off()
table(mergedColors_vscc)
nSamples_vscc <- nrow(datExpr_filtered_conn_vscc)
MEs0_vscc <- moduleEigengenes(datExpr_filtered_conn_vscc, mergedColors_vscc)$eigengenes
MEs_vscc <- orderMEs(MEs0_vscc)
traitData_num_vscc <- data.frame(
  HPV_Status = as.numeric(factor(traitData_vscc$HPV_Status))
)
moduleTraitCor_vscc <- cor(MEs_vscc, traitData_num_vscc, use = "p")
moduleTraitPvalue_vscc <- corPvalueStudent(moduleTraitCor_vscc, nSamples_vscc)
# Clean module names for display
module_names_clean_vscc <- sub("^ME", "", names(MEs_vscc))
# Combine correlation and p-values
textMatrix_vscc <- paste0(signif(moduleTraitCor_vscc, 2), "\n(", signif(moduleTraitPvalue_vscc, 1), ")")
# Plot correlation heatmap
png("module_trait_relationships_vscc.png", width = 9, height = 7, units = "in", res = 900)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor_vscc,
  xLabels = names(traitData_num_vscc),
  yLabels = module_names_clean_vscc,
  ySymbols = module_names_clean_vscc,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_vscc,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "VSCC Module - Trait Relationships"
)
dev.off()

