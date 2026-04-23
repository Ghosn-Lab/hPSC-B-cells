## vlnplot_by_batch.R
# Generate violin plots for a given stage by batch. 
# Identify the pre-BCR program genes as dysregulated

library(Seurat)
library(dplyr)


## -------------------------------------------------------------------------
## Pre-BCR violin plots by batch
## -------------------------------------------------------------------------

make_vlnplot_by_batch <- function(seurat_obj, cell_label, batch_names, plot_colors) {
  obj_sub <- subset(seurat_obj, subset = custom_labels == cell_label)
  obj_sub <- subset(obj_sub, subset = batch %in% batch_names)
  
  obj_sub$batch <- factor(obj_sub$batch, levels = batch_names)
  
  plot <- VlnPlot(
    obj_sub,
    features = c("CD79A", "CD9", "IKZF3", "SPIB", "FCRLA", "TCL1A"),
    split.by = "batch",
    split.plot = TRUE,
    cols = plot_colors,
    pt.size = 0
  )
  
  return(plot)
}

plot_il7_high <- make_vlnplot_by_batch(
  seurat_obj = combined,
  cell_label = "Large pre-B",
  batch_names = c("BM_total", "IL7_high_hPSC"),
  plot_colors = c("#D537B7", "steelblue")
)

ggsave(
  filename = "figures/vlnplot_bcr_genes_lpb_bm_vs_il7_high.svg",
  plot     = plot_il7_high,
  device   = "svg",
  bg       = "transparent",
  width    = 6,
  height   = 5,
  units    = "in"
)

plot_il7_low <- make_vlnplot_by_batch(
  seurat_obj = combined,
  cell_label = "Large pre-B",
  batch_names = c("BM_total", "IL7_low_hPSC"),
  plot_colors = c("#D537B7", "turquoise4")
)

ggsave(
  filename = "figures/vlnplot_bcr_genes_lpb_bm_vs_il7_low.svg",
  plot     = plot_il7_low,
  device   = "svg",
  bg       = "transparent",
  width    = 6,
  height   = 5,
  units    = "in"
)