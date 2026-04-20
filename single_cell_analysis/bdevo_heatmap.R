## bdevo_heatmap.R
## Create a heatmap showing % cells expressing canonical B-cell developement
## genes tracking B-cell stages and accross different cytokine conditions

library(Seurat)
library(dplyr)
library(stringr)
library(pheatmap)


## -------------------------------------------------------------------------
## Marker heatmaps by group
## -------------------------------------------------------------------------

# B-cell development genes, ordered by stage progression
gene_set <- c(
  "EBF1", "PAX5", "CD19", "CD34", "ERG", "DNTT", "MKI67", "RAG1", "RAG2",
  "MME", "CD9", "IKZF3", "SPIB", "FCRLA", "TCL1A", "MS4A1", "CD83",
  "CD40", "FCER2", "CXCR5"
)

cluster_order <- c(
  "Early pro-B",
  "Late pro-B",
  "Large pre-B",
  "Small pre-B",
  "Immature B",
  "Maturing B"
)

plot_marker_heatmap <- function(seurat_obj, group_name, file_out,
                                gene_set, cluster_order) {
  
  expr_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[gene_set, , drop = FALSE]
  expr_bin <- expr_mat > 0
  
  marker_pct <- sapply(cluster_order, function(cl) {
    cells <- WhichCells(seurat_obj, idents = cl)
    rowMeans(expr_bin[, cells, drop = FALSE]) * 100
  })
  
  marker_pct <- as.matrix(marker_pct)
  rownames(marker_pct) <- gene_set
  colnames(marker_pct) <- cluster_order
  
  plot <- pheatmap(
    marker_pct,
    color = colorRampPalette(c("white", "#000040"))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    border_color = "grey80",
    main = paste0("% cells expressing marker per cluster\n", group_name)
  )
  
  ggsave(
    filename = file_out,
    plot     = plot,
    device   = "svg",
    bg       = "transparent",
    width    = 3,
    height   = 7,
    units    = "in"
  )
}

## -------------------------------------------------------------------------
## Generate/save plots for each batch
## -------------------------------------------------------------------------

plot_marker_heatmap(
  bm_total,
  "BM_total",
  "figures/BM_total_b_development.svg",
  gene_set,
  cluster_order
)

plot_marker_heatmap(
  il7_high,
  "il7_high",
  "figures/il7_high_b_development.svg",
  gene_set,
  cluster_order
)

plot_marker_heatmap(
  il7_low,
  "il7_low",
  "figures/il7_low_b_development.svg",
  gene_set,
  cluster_order
)