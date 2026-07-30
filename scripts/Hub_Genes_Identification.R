############################################################
# Hub Gene Identification
############################################################
# ==============================================
# HNC
# ==============================================
module_of_interest_1_hnc <- "blue"
module_of_interest_2_hnc <- "turquoise"

# Convert HPV status to numeric for correlation
trait_of_interest_hnc <- as.numeric(factor(traitData_hnc$HPV_Status))

# Gene Significance (GS) and p-value
geneTraitSignificance_hnc <- as.data.frame(cor(datExpr_filtered_conn_hnc, trait_of_interest_hnc, use = "p"))
names(geneTraitSignificance_hnc) <- "GS.HPV"
GSPvalue_hnc <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_hnc), nSamples_hnc))
names(GSPvalue_hnc) <- "GSPvalue"

# Module Membership (MM) and p-value
geneModuleMembership_hnc <- as.data.frame(cor(datExpr_filtered_conn_hnc, MEs_hnc, use = "p"))
MMPvalue_hnc <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_hnc), nSamples_hnc))
names(geneModuleMembership_hnc) <- paste0("MM.", sub("^ME", "", names(MEs_hnc)))
names(MMPvalue_hnc) <- paste0("MMP.", sub("^ME", "", names(MEs_hnc)))

# Genes in selected modules
module_genes_blue_hnc <- mergedColors_hnc == module_of_interest_1_hnc
module_genes_turquoise_hnc <- mergedColors_hnc == module_of_interest_2_hnc

# Thresholds
mm_thresh_hnc <- 0.8
gs_thresh_hnc <- 0.5

# Identify hub genes (high MM and GS)
hub_genes_blue_hnc <- colnames(datExpr_filtered_conn_hnc)[
  module_genes_blue_hnc &
    abs(geneModuleMembership_hnc[, paste0("MM.", module_of_interest_1_hnc)]) > mm_thresh_hnc &
    abs(geneTraitSignificance_hnc$GS.HPV) > gs_thresh_hnc
]
hub_genes_turquoise_hnc <- colnames(datExpr_filtered_conn_hnc)[
  module_genes_turquoise_hnc &
    abs(geneModuleMembership_hnc[, paste0("MM.", module_of_interest_2_hnc)]) > mm_thresh_hnc &
    abs(geneTraitSignificance_hnc$GS.HPV) > gs_thresh_hnc
]

# Connect to Ensembl
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Map Ensembl IDs to Gene Symbols
id_map_1_hnc <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = hub_genes_blue_hnc,
  mart = mart
)
id_map_2_hnc <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = hub_genes_turquoise_hnc,
  mart = mart
)

# Handle missing symbols
hub_symbols_blue_hnc <- id_map_1_hnc$hgnc_symbol[match(hub_genes_blue_hnc, id_map_1_hnc$ensembl_gene_id)]
hub_symbols_blue_hnc[is.na(hub_symbols_blue_hnc)] <- hub_genes_blue_hnc[is.na(hub_symbols_blue_hnc)]
hub_symbols_turquoise_hnc <- id_map_2_hnc$hgnc_symbol[match(hub_genes_turquoise_hnc, id_map_2_hnc$ensembl_gene_id)]
hub_symbols_turquoise_hnc[is.na(hub_symbols_turquoise_hnc)] <- hub_genes_turquoise_hnc[is.na(hub_symbols_turquoise_hnc)]

hub_table_blue_hnc <- data.frame(
  EnsemblID = hub_genes_blue_hnc,
  GeneSymbol = hub_symbols_blue_hnc,
  MM = geneModuleMembership_hnc[hub_genes_blue_hnc, paste0("MM.", module_of_interest_1_hnc)],
  MMPvalue = MMPvalue_hnc[hub_genes_blue_hnc, paste0("MMP.", module_of_interest_1_hnc)],
  GS = geneTraitSignificance_hnc[hub_genes_blue_hnc, "GS.HPV"],
  GSPvalue = GSPvalue_hnc[hub_genes_blue_hnc, "GSPvalue"]
)
hub_table_turquoise_hnc <- data.frame(
  EnsemblID = hub_genes_turquoise_hnc,
  GeneSymbol = hub_symbols_turquoise_hnc,
  MM = geneModuleMembership_hnc[hub_genes_turquoise_hnc, paste0("MM.", module_of_interest_2_hnc)],
  MMPvalue = MMPvalue_hnc[hub_genes_turquoise_hnc, paste0("MMP.", module_of_interest_2_hnc)],
  GS = geneTraitSignificance_hnc[hub_genes_turquoise_hnc, "GS.HPV"],
  GSPvalue = GSPvalue_hnc[hub_genes_turquoise_hnc, "GSPvalue"]
)

combined_hubs_hnc <- bind_rows(hub_table_blue_hnc, hub_table_turquoise_hnc)
combined_hubs_scored_hnc <- combined_hubs_hnc %>%
  mutate(combined_score = abs(MM * GS)) %>%
  arrange(desc(combined_score))

write.csv(combined_hubs_scored_hnc, "combined hubs of hnc.csv", row.names = FALSE)
write.csv(hub_table_blue_hnc, "blue hub genes of hnc.csv", row.names = FALSE)
write.csv(hub_table_turquoise_hnc, "turquoise hub genes of hnc.csv", row.names = FALSE)

# ==============================================
# OPSCC
# ==============================================
module_of_interest_1_opscc <- "black"
module_of_interest_2_opscc <- "brown"

# Convert HPV status to numeric for correlation
trait_of_interest_opscc <- as.numeric(factor(traitData_opscc$HPV_Status))

# Gene Significance (GS) and p-value
geneTraitSignificance_opscc <- as.data.frame(cor(datExpr_filtered_conn_opscc, trait_of_interest_opscc, use = "p"))
names(geneTraitSignificance_opscc) <- "GS.HPV"
GSPvalue_opscc <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_opscc), nSamples_opscc))
names(GSPvalue_opscc) <- "GSPvalue"

# Module Membership (MM) and p-value
geneModuleMembership_opscc <- as.data.frame(cor(datExpr_filtered_conn_opscc, MEs_opscc, use = "p"))
MMPvalue_opscc <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_opscc), nSamples_opscc))
names(geneModuleMembership_opscc) <- paste0("MM.", sub("^ME", "", names(MEs_opscc)))
names(MMPvalue_opscc) <- paste0("MMP.", sub("^ME", "", names(MEs_opscc)))

# Genes in selected modules
module_genes_black_opscc <- mergedColors_opscc == module_of_interest_1_opscc
module_genes_brown_opscc <- mergedColors_opscc == module_of_interest_2_opscc

# Thresholds
mm_thresh_opscc <- 0.8
gs_thresh_opscc <- 0.5

# Identify hub genes (high MM and GS)
hub_genes_black_opscc <- colnames(datExpr_filtered_conn_opscc)[
  module_genes_black_opscc &
    abs(geneModuleMembership_opscc[, paste0("MM.", module_of_interest_1_opscc)]) > mm_thresh_opscc &
    abs(geneTraitSignificance_opscc$GS.HPV) > gs_thresh_opscc
]
hub_genes_brown_opscc <- colnames(datExpr_filtered_conn_opscc)[
  module_genes_brown_opscc &
    abs(geneModuleMembership_opscc[, paste0("MM.", module_of_interest_2_opscc)]) > mm_thresh_opscc &
    abs(geneTraitSignificance_opscc$GS.HPV) > gs_thresh_opscc
]

# Connect to Ensembl
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Map Ensembl IDs to Gene Symbols
id_map_1_opscc <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = hub_genes_black_opscc,
  mart = mart
)
id_map_2_opscc <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = hub_genes_brown_opscc,
  mart = mart
)

# Handle missing symbols
hub_symbols_black_opscc <- id_map_1_opscc$hgnc_symbol[match(hub_genes_black_opscc, id_map_1_opscc$ensembl_gene_id)]
hub_symbols_black_opscc[is.na(hub_symbols_black_opscc)] <- hub_genes_black_opscc[is.na(hub_symbols_black_opscc)]
hub_symbols_brown_opscc <- id_map_2_opscc$hgnc_symbol[match(hub_genes_brown_opscc, id_map_2_opscc$ensembl_gene_id)]
hub_symbols_brown_opscc[is.na(hub_symbols_brown_opscc)] <- hub_genes_brown_opscc[is.na(hub_symbols_brown_opscc)]

hub_table_black_opscc <- data.frame(
  EnsemblID = hub_genes_black_opscc,
  GeneSymbol = hub_symbols_black_opscc,
  MM = geneModuleMembership_opscc[hub_genes_black_opscc, paste0("MM.", module_of_interest_1_opscc)],
  MMPvalue = MMPvalue_opscc[hub_genes_black_opscc, paste0("MMP.", module_of_interest_1_opscc)],
  GS = geneTraitSignificance_opscc[hub_genes_black_opscc, "GS.HPV"],
  GSPvalue = GSPvalue_opscc[hub_genes_black_opscc, "GSPvalue"]
)
hub_table_brown_opscc <- data.frame(
  EnsemblID = hub_genes_brown_opscc,
  GeneSymbol = hub_symbols_brown_opscc,
  MM = geneModuleMembership_opscc[hub_genes_brown_opscc, paste0("MM.", module_of_interest_2_opscc)],
  MMPvalue = MMPvalue_opscc[hub_genes_brown_opscc, paste0("MMP.", module_of_interest_2_opscc)],
  GS = geneTraitSignificance_opscc[hub_genes_brown_opscc, "GS.HPV"],
  GSPvalue = GSPvalue_opscc[hub_genes_brown_opscc, "GSPvalue"]
)

combined_hubs_opscc <- bind_rows(hub_table_black_opscc, hub_table_brown_opscc)
combined_hubs_scored_opscc <- combined_hubs_opscc %>%
  mutate(combined_score = abs(MM * GS)) %>%
  arrange(desc(combined_score))

write.csv(combined_hubs_scored_opscc, "combined hubs of opscc.csv", row.names = FALSE)
write.csv(hub_table_black_opscc, "black hub genes of opscc.csv", row.names = FALSE)
write.csv(hub_table_brown_opscc, "brown hub genes of opscc.csv", row.names = FALSE)

# ==============================================
# VSCC
# ==============================================
module_of_interest_1_vscc <- "green"
module_of_interest_2_vscc <- "black"

# Convert HPV status to numeric for correlation
trait_of_interest_vscc <- as.numeric(factor(traitData_vscc$HPV_Status))

# Gene Significance (GS) and p-value
geneTraitSignificance_vscc <- as.data.frame(cor(datExpr_filtered_conn_vscc, trait_of_interest_vscc, use = "p"))
names(geneTraitSignificance_vscc) <- "GS.HPV"
GSPvalue_vscc <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_vscc), nSamples_vscc))
names(GSPvalue_vscc) <- "GSPvalue"

# Module Membership (MM) and p-value
geneModuleMembership_vscc <- as.data.frame(cor(datExpr_filtered_conn_vscc, MEs_vscc, use = "p"))
MMPvalue_vscc <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_vscc), nSamples_vscc))
names(geneModuleMembership_vscc) <- paste0("MM.", sub("^ME", "", names(MEs_vscc)))
names(MMPvalue_vscc) <- paste0("MMP.", sub("^ME", "", names(MEs_vscc)))

# Genes in selected modules
module_genes_green_vscc <- mergedColors_vscc == module_of_interest_1_vscc
module_genes_black_vscc <- mergedColors_vscc == module_of_interest_2_vscc

# Thresholds
mm_thresh_vscc <- 0.8
gs_thresh_vscc <- 0.5

# Identify hub genes (high MM and GS)
hub_genes_green_vscc <- colnames(datExpr_filtered_conn_vscc)[
  module_genes_green_vscc &
    abs(geneModuleMembership_vscc[, paste0("MM.", module_of_interest_1_vscc)]) > mm_thresh_vscc &
    abs(geneTraitSignificance_vscc$GS.HPV) > gs_thresh_vscc
]
hub_genes_black_vscc <- colnames(datExpr_filtered_conn_vscc)[
  module_genes_black_vscc &
    abs(geneModuleMembership_vscc[, paste0("MM.", module_of_interest_2_vscc)]) > mm_thresh_vscc &
    abs(geneTraitSignificance_vscc$GS.HPV) > gs_thresh_vscc
]

# Connect to Ensembl
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Map Ensembl IDs to Gene Symbols
id_map_1_vscc <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = hub_genes_green_vscc,
  mart = mart
)
id_map_2_vscc <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = hub_genes_black_vscc,
  mart = mart
)

# Handle missing symbols
hub_symbols_green_vscc <- id_map_1_vscc$hgnc_symbol[match(hub_genes_green_vscc, id_map_1_vscc$ensembl_gene_id)]
hub_symbols_green_vscc[is.na(hub_symbols_green_vscc)] <- hub_genes_green_vscc[is.na(hub_symbols_green_vscc)]
hub_symbols_black_vscc <- id_map_2_vscc$hgnc_symbol[match(hub_genes_black_vscc, id_map_2_vscc$ensembl_gene_id)]
hub_symbols_black_vscc[is.na(hub_symbols_black_vscc)] <- hub_genes_black_vscc[is.na(hub_symbols_black_vscc)]

hub_table_green_vscc <- data.frame(
  EnsemblID = hub_genes_green_vscc,
  GeneSymbol = hub_symbols_green_vscc,
  MM = geneModuleMembership_vscc[hub_genes_green_vscc, paste0("MM.", module_of_interest_1_vscc)],
  MMPvalue = MMPvalue_vscc[hub_genes_green_vscc, paste0("MMP.", module_of_interest_1_vscc)],
  GS = geneTraitSignificance_vscc[hub_genes_green_vscc, "GS.HPV"],
  GSPvalue = GSPvalue_vscc[hub_genes_green_vscc, "GSPvalue"]
)
hub_table_black_vscc <- data.frame(
  EnsemblID = hub_genes_black_vscc,
  GeneSymbol = hub_symbols_black_vscc,
  MM = geneModuleMembership_vscc[hub_genes_black_vscc, paste0("MM.", module_of_interest_2_vscc)],
  MMPvalue = MMPvalue_vscc[hub_genes_black_vscc, paste0("MMP.", module_of_interest_2_vscc)],
  GS = geneTraitSignificance_vscc[hub_genes_black_vscc, "GS.HPV"],
  GSPvalue = GSPvalue_vscc[hub_genes_black_vscc, "GSPvalue"]
)

combined_hubs_vscc <- bind_rows(hub_table_green_vscc, hub_table_black_vscc)
combined_hubs_scored_vscc <- combined_hubs_vscc %>%
  mutate(combined_score = abs(MM * GS)) %>%
  arrange(desc(combined_score))

write.csv(combined_hubs_scored_vscc, "combined hubs of vscc.csv", row.names = FALSE)
write.csv(hub_table_green_vscc, "green hub genes of vscc.csv", row.names = FALSE)
write.csv(hub_table_black_vscc, "black hub genes of vscc.csv", row.names = FALSE)



nrow(combined_hubs_hnc)
nrow(combined_hubs_opscc)
nrow(combined_hubs_vscc)
