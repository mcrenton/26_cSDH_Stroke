# =============================================================================
# Figure 1A - cSDH fluid protein heatmap (HORIZONTAL layout)
#   Patients = rows, analytes = columns. Analytes ordered by hierarchical
#   clustering (dendrogram on top); patients unclustered.
#   ONLY analytes quantified (>= LLoQ) in ALL 5 patients are shown, so every
#   cell is fully quantitative -- no ND or below-LLoQ cells remain.
# =============================================================================

## ---- 1. Packages ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")
if (!requireNamespace("viridisLite", quietly = TRUE)) install.packages("viridisLite")
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(grid)
## ---- 2. Load data -----------------------------------------------------------
# setwd("~/Documents/Claude/Projects/cSDH paper")   # uncomment / edit if needed
df <- read.csv("csdh_heatmap_data.csv", stringsAsFactors = FALSE, check.names = FALSE)
patient_levels  <- c("Pt 43", "Pt 44", "Pt 45", "Pt 74", "Pt 77")
category_levels <- unique(df[order(df$category_order), "category"])
## ---- 3. Reshape to matrices -------------------------------------------------
analytes <- unique(df$analyte[order(df$category_order)])
conc_mat   <- matrix(NA_real_, nrow = length(analytes), ncol = length(patient_levels),
                     dimnames = list(analytes, patient_levels))
status_mat <- matrix("ND",     nrow = length(analytes), ncol = length(patient_levels),
                     dimnames = list(analytes, patient_levels))
lod_vec    <- setNames(rep(NA_real_, length(analytes)), analytes)
for (i in seq_len(nrow(df))) {
  a <- df$analyte[i]; p <- df$patient[i]
  status_mat[a, p] <- df$status[i]
  if (!is.na(df$concentration[i]) && df$concentration[i] != "")
    conc_mat[a, p] <- as.numeric(df$concentration[i])
  if (!is.na(df$LoD[i]) && df$LoD[i] != "") lod_vec[a] <- as.numeric(df$LoD[i])
}
# --- Keep only analytes quantified (>= LLoQ) in ALL 5 patients ---------------
keep <- analytes[apply(status_mat[analytes, patient_levels, drop = FALSE], 1,
                       function(s) all(s == "QUANT"))]
analytes   <- keep
conc_mat   <- conc_mat[analytes, , drop = FALSE]
status_mat <- status_mat[analytes, , drop = FALSE]
lod_vec    <- lod_vec[analytes]
# log10 transform for the color scale (values span several orders of magnitude)
log_mat <- log10(conc_mat)
## ---- 4. Hierarchical clustering ---------------------------------------------
# Every cell is quantitative (no NAs), so cluster on log_mat directly.
# --- Clustering options (defaults: euclidean distance, complete linkage) ------
# To cluster analytes by co-variation PATTERN across patients rather than by
# absolute abundance, uncomment the scaling line below (row z-scores).
clust_input <- log_mat
# clust_input <- t(scale(t(log_mat)))   # <- optional: per-analyte z-score
analyte_dend <- as.dendrogram(hclust(dist(clust_input, method = "euclidean"), method = "complete"))
## ---- 5. Colors --------------------------------------------------------------
rng     <- range(log_mat, na.rm = TRUE)
col_fun <- colorRamp2(seq(rng[1], rng[2], length.out = 256), magma(256))
## ---- 7. Build heatmap (transposed = horizontal) -----------------------------
ht <- Heatmap(
  t(log_mat),                       # patients = rows, analytes = columns
  name               = "log10\n(pg/mL)",
  col                = col_fun,
  na_col             = "#dcdcdc",
  cluster_rows       = FALSE,        # patients: keep given order
  cluster_columns    = analyte_dend, # analytes: clustered, dendrogram on top
  column_dend_height = unit(15, "mm"),
  column_names_rot   = 90,           # analyte labels vertical along the bottom
  column_names_gp    = gpar(fontsize = 7.5),
  column_names_side  = "bottom",
  row_names_gp       = gpar(fontsize = 10),
  row_names_side     = "left",
  border             = TRUE,
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(3, "cm")
  ),
  column_title = "cSDH fluid protein profile (n = 5 patients; hierarchical clustering)\nAnalytes quantified >= LLoQ in all 5 patients",
  column_title_gp = gpar(fontsize = 10, fontface = "bold")
)
## ---- 8. Export (landscape) --------------------------------------------------
pdf("Figure1A_heatmap.pdf", width = 14, height = 5)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legends = TRUE, padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()
png("Figure1A_heatmap.png", width = 14, height = 5, units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legends = TRUE, padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()
message("Done: Figure1A_heatmap.pdf and Figure1A_heatmap.png written.")
