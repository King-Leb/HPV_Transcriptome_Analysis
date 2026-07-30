############################################################
# 4. Meta-Analysis (cross-cohort comparison, uses sig_matrix)
############################################################
# --------------------------------------
# SCATTERPLOT
# --------------------------------------
get_scatter_plot <- function(mat, col_x, col_y, top_n = 100,
                             x_label = NULL, y_label = NULL, main_title = NULL,
                             point_col = "gray85", up_col = "forestgreen", down_col = "red2",
                             text_size = 3) {
  
  x <- mat[[col_x]]
  y <- mat[[col_y]]
  df <- data.frame(x = x, y = y, gene = rownames(mat))
  df <- df[complete.cases(df), ]
  
  # --- HIGHLIGHTING LOGIC ---
  top_up_x <- df$gene[order(df$x, decreasing = TRUE)[1:min(top_n, nrow(df))]]
  top_down_x <- df$gene[order(df$x, decreasing = FALSE)[1:min(top_n, nrow(df))]]
  top_up_y <- df$gene[order(df$y, decreasing = TRUE)[1:min(top_n, nrow(df))]]
  top_down_y <- df$gene[order(df$y, decreasing = FALSE)[1:min(top_n, nrow(df))]]
  
  df$highlight <- "NS"
  df$highlight[df$gene %in% c(top_up_x, top_up_y)] <- "upregulated"
  df$highlight[df$gene %in% c(top_down_x, top_down_y)] <- "downregulated"
  df$highlight <- factor(df$highlight, levels = c("upregulated", "downregulated", "NS"))
  
  extreme_genes <- c(head(df$gene[order(df$x + df$y, decreasing = TRUE)], 15),
                     head(df$gene[order(df$x + df$y, decreasing = FALSE)], 15))
  

  cor_test <- cor.test(df$x, df$y, method = "pearson")
  r_val    <- unname(cor_test$estimate)
  p_val    <- cor_test$p.value
  n_genes  <- nrow(df)
  
  pval_text <- if (p_val < 2.2e-16) {
    "< 2.2e-16"
  } else {
    signif(p_val, 3)
  }
  

  spearman_test <- suppressWarnings(cor.test(df$x, df$y, method = "spearman"))
  

  effect_flag <- if (abs(r_val) < 0.1) " *** NEGLIGIBLE EFFECT SIZE (p reflects large n) ***" else ""
  
  message(sprintf(
    "[%s vs %s] Pearson r = %.3f (p = %.3g) | Spearman rho = %.3f (p = %.3g) | n = %d genes%s",
    col_x, col_y, r_val, p_val, unname(spearman_test$estimate), spearman_test$p.value, n_genes, effect_flag
  ))
  
  if(is.null(x_label)) x_label <- col_x
  if(is.null(y_label)) y_label <- col_y
  if(is.null(main_title)) main_title <- paste(col_x, "vs", col_y)
  
  # --- PLOT ---
  p <- ggplot(df, aes(x = x, y = y, color = highlight)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(values = c("upregulated" = up_col, "downregulated" = down_col, "NS" = point_col),
                       name = "Regulation") +
    geom_text_repel(data = subset(df, gene %in% extreme_genes),
                    aes(label = gene), size = text_size, max.overlaps = 30) +
    geom_smooth(method = "lm", color = "midnightblue", se = FALSE, size = 1) +
    annotate("label", x = min(df$x), y = max(df$y),
             label = paste0("r = ", round(r_val, 3), "\np = ", pval_text, "\nn = ", n_genes),
             hjust = 0, vjust = 1, size = 3.5, fill = "white", alpha = 0.8) +
    labs(x = x_label, y = y_label, title = main_title) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),   # remove grid lines
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}

# Assemble and save
s_hnc_opscc <- get_scatter_plot(sig_matrix, "HNC", "OPSCC")
s_hnc_vscc  <- get_scatter_plot(sig_matrix, "HNC", "VSCC")
s_opscc_vscc <- get_scatter_plot(sig_matrix, "OPSCC", "VSCC")

# Create 2x2 grid (A, B on top, C on bottom left)
combined_scatter <- (s_hnc_opscc + s_hnc_vscc) / (s_opscc_vscc + plot_spacer()) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'right')

# High DPI and vector exports
ggsave("Scatter_Grid_600dpi.png", combined_scatter, width = 12, height = 12, dpi = 600)
ggsave("Scatter_Grid.pdf", combined_scatter, width = 12, height = 12)
ggsave("Scatter_Grid.svg", combined_scatter, width = 12, height = 12)




# ----------------------------------------------
# INTEGRATION
# ----------------------------------------------

sig_matrix_complete <- sig_matrix[complete.cases(sig_matrix), ]

# Step 1: compute Stouffer coefficient across 3 tumors
stouffer <- function(x){
  z <- sum(x) / sqrt(length(x))
  return(z)
}
stouffer_coefficient <- apply(sig_matrix_complete, 1, stouffer)

# Step 2: compute SD across 3 tumors
sd_genes <- apply(sig_matrix_complete, 1, sd)

# Step 3: create integrated matrix
matIntegration <- data.frame(
  gene = rownames(sig_matrix_complete),
  stouffer = stouffer_coefficient,
  sd = sd_genes
)
write.csv(matIntegration, file = "Matrix Integration.csv")

# Step 4: classify genes by stouffer and sd
summary(matIntegration$stouffer)
summary(matIntegration$sd)
sd_sd <- sd(matIntegration$sd, na.rm = TRUE)
sd_sd
quant_filter <- quantile(matIntegration$sd, 0.9, na.rm = TRUE)
quant_filter

# filter at the 90th percentile, 2.131126
# Adaptive thresholds
low_sd_thresh <- quant_filter   # from quantile
stouffer_thresh <- 5

# Classification
matIntegration$expression <- 'NS'
matIntegration$expression[matIntegration$sd < low_sd_thresh & matIntegration$stouffer < -stouffer_thresh] <- 'Down & Low SD'
matIntegration$expression[matIntegration$sd < low_sd_thresh & matIntegration$stouffer >  stouffer_thresh] <- 'Up & Low SD'
matIntegration$expression[matIntegration$sd >= low_sd_thresh & matIntegration$stouffer < -stouffer_thresh] <- 'Down & High SD'
matIntegration$expression[matIntegration$sd >= low_sd_thresh & matIntegration$stouffer >  stouffer_thresh] <- 'Up & High SD'

# Step 5: labels
matIntegration$label <- NA
matIntegration$label[matIntegration$expression != 'NS'] <- matIntegration$gene[matIntegration$expression != 'NS']

# Define colors
expr_colors <- c(
  'Down & Low SD'='royalblue',
  'Up & Low SD'='firebrick1',
  'Down & High SD'='darkturquoise',
  'Up & High SD'='orange',
  'NS' = 'grey70'
)

# Base plot
integration_plot <- ggplot(matIntegration, aes(x = stouffer, y = sd, color = expression, label = label)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(
    size = 3,
    max.overlaps = 30,
    box.padding = 0.3,
    segment.color = 'grey50'
  ) +
  scale_color_manual(values = expr_colors) +
  labs(
    title = 'Genes deregulated across all tumors',
    x = 'Stouffer Coefficient',
    y = 'Standard Deviation',
    color = 'expression'  # removes default legend title
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save high-res PNG
ggsave("Integration_Final.png",
       plot = integration_plot,
       width = 10, height = 6, units = "in",
       dpi = 600)
# Save PDF
ggsave("Integration_Final.pdf",
       plot = integration_plot,
       width = 10, height = 6, units = "in")
# Save SVG
ggsave("Integration_Final.svg",
       plot = integration_plot,
       width = 10, height = 6, units = "in")

# Keep only genes classified as up/down & low SD
up_low  <- filter(matIntegration, expression == 'Up & Low SD')
dn_low  <- filter(matIntegration, expression == 'Down & Low SD')
up_high <- filter(matIntegration, expression == 'Up & High SD')
dn_high <- filter(matIntegration, expression == 'Down & High SD')

# Extract the corresponding rows from the signature matrix
genes_low  <- sig_matrix[rownames(sig_matrix) %in% c(up_low$gene, dn_low$gene), ]
genes_high <- sig_matrix[rownames(sig_matrix) %in% c(up_high$gene, dn_high$gene), ]

# Replace NA with 0
genes_low[is.na(genes_low)] <- 0
genes_high[is.na(genes_high)] <- 0

table(matIntegration$expression)
write.csv(matIntegration, file = "Integration.csv")

save(deseq_hnc, file = "deseq_hnc.rda")
save(deseq_opscc, file = "deseq_opscc.rda")
save(deseq_vscc, file = "deseq_vscc.rda")
write.csv(sig_matrix, file = "sig_matrix_data.csv")






############################################################
# Sensitivity Analysis & Multiple-Testing Correction
############################################################

matIntegration$meta_p    <- 2 * pnorm(-abs(matIntegration$stouffer))
matIntegration$meta_padj <- p.adjust(matIntegration$meta_p, method = "BH")


# ------------------------------------------------------------
primary_Z_thresh  <- 5
primary_SD_thresh <- quantile(matIntegration$sd, 0.90, na.rm = TRUE)

p_at_Z5 <- 2 * pnorm(-abs(primary_Z_thresh))

genes_at_primary <- matIntegration[abs(matIntegration$stouffer) > primary_Z_thresh, ]
# least-significant (largest) adjusted p-value among genes just clearing the |Z|>5 boundary
padj_at_primary_boundary <- if (nrow(genes_at_primary) > 0) {
  max(genes_at_primary$meta_padj)
} else {
  NA
}

primary_signature <- matIntegration$gene[
  abs(matIntegration$stouffer) > primary_Z_thresh &
    matIntegration$sd < primary_SD_thresh
]
n_up_primary   <- sum(matIntegration$stouffer[match(primary_signature, matIntegration$gene)] > 0)
n_down_primary <- length(primary_signature) - n_up_primary

cat("--- PRIMARY THRESHOLD SUMMARY (fill into Methods text) ---\n")
cat("[X] Z threshold                         :", primary_Z_thresh, "\n")
cat("    Nominal p at |Z| =", primary_Z_thresh, "             :", format(p_at_Z5, scientific = TRUE), "\n")
cat("[Z] Total genes in cross-cohort programme (the '413' number) :", nrow(genes_at_primary), "\n")
cat("[Y] Least-significant adjusted p among retained programme genes:",
    format(padj_at_primary_boundary, scientific = TRUE), "\n")
cat("[W] SD threshold (90th percentile)      :", round(as.numeric(primary_SD_thresh), 3), "\n")
cat("[A] Up-regulated in Pan-HPV signature    :", n_up_primary, "\n")
cat("[B] Down-regulated in Pan-HPV signature  :", n_down_primary, "\n")
cat("[A+B] Total Pan-HPV signature size       :", length(primary_signature), "\n\n")

# ------------------------------------------------------------
# Step 3: Sensitivity analysis -- vary |Z| and SD-percentile
# ------------------------------------------------------------
z_thresholds   <- c(4, 4.5, 5, 5.5, 6)
sd_percentiles <- c(0.85, 0.90, 0.95)

jaccard <- function(a, b) {
  if (length(union(a, b)) == 0) return(NA)
  length(intersect(a, b)) / length(union(a, b))
}

param_grid <- expand.grid(Z_threshold = z_thresholds, SD_percentile = sd_percentiles)
sensitivity_results <- data.frame()

for (i in seq_len(nrow(param_grid))) {
  z_thr  <- param_grid$Z_threshold[i]
  sd_pct <- param_grid$SD_percentile[i]
  sd_thr <- quantile(matIntegration$sd, sd_pct, na.rm = TRUE)
  
  programme_genes <- matIntegration$gene[abs(matIntegration$stouffer) > z_thr]
  signature_genes <- matIntegration$gene[
    abs(matIntegration$stouffer) > z_thr & matIntegration$sd < sd_thr
  ]
  n_up   <- sum(matIntegration$stouffer[match(signature_genes, matIntegration$gene)] > 0)
  n_down <- length(signature_genes) - n_up
  
  sensitivity_results <- rbind(sensitivity_results, data.frame(
    Z_threshold        = z_thr,
    SD_percentile       = sd_pct,
    SD_threshold        = round(as.numeric(sd_thr), 3),
    programme_n         = length(programme_genes),
    signature_n         = length(signature_genes),
    n_up                = n_up,
    n_down              = n_down,
    jaccard_vs_primary  = jaccard(signature_genes, primary_signature)
  ))
}

cat("--- SENSITIVITY ANALYSIS RESULTS ---\n")
print(sensitivity_results, row.names = FALSE)

write.csv(
  sensitivity_results,
  "Supplementary_Table_Sensitivity_Analysis_HPV_Signature.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Step 4: Export the full matIntegration table with p / padj
# ------------------------------------------------------------
write.csv(
  matIntegration,
  "Matrix_Integration_with_adjusted_pvalues.csv",
  row.names = FALSE
)

cat("\nDone. Two files written to the working directory:\n")
cat(" - Supplementary_Table_Sensitivity_Analysis_Pan_HPV_Signature.csv\n")
cat(" - Matrix_Integration_with_adjusted_pvalues.csv\n")



