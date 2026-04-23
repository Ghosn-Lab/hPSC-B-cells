## featureplot_by_batch.R
# Generate feature plots for a given batch. Over here we look at canonical
# B-cell developement genes, showing maturation trajectory

library(Seurat)


## -------------------------------------------------------------------------
## Maturation feature plots for IL-7 low batch
## -------------------------------------------------------------------------

make_featureplots_by_batch <- function(seurat_obj, batch_name, genes,
                                       reduction = "umap_flipped",
                                       out_dir = "figures",
                                       width = 4, height = 3.5) {
  obj_sub <- subset(seurat_obj, subset = batch == batch_name)
  
  for (gene in genes) {
    plot <- FeaturePlot(
      obj_sub,
      features = gene,
      reduction = reduction,
      order = TRUE
    )
    
    ggsave(
      filename = file.path(out_dir, paste0(batch_name, "_featureplot_", gene, ".svg")),
      plot = plot,
      device = "svg",
      bg = "transparent",
      width = width,
      height = height,
      units = "in"
    )
  }
}


development_genes <- c(
  "CD19", "CD34", "MKI67", "DNTT", "RAG1",
  "IGF2", "CD40", "MS4A1", "CXCR5", "MME"
)

make_featureplots_by_batch(
  seurat_obj = combined,
  batch_name = "il7_low",
  genes = development_genes
)