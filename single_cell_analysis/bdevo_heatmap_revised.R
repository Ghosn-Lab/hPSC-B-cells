## bdevo_heatmap.R
## Create a heatmap showing % cells expressing canonical B-cell developement
## genes tracking B-cell stages and accross different cytokine conditions

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(pheatmap)
library(svglite)
library(grid)

plot_marker_phase_heatmap <- function(
    seurat_obj,
    sample_name = "Sample",
    cluster_col = "custom_labels",
    phase_col = "Phase",
    assay = "RNA",
    layer = "data",
    output_file = NULL,
    width = 3.3,
    height = 7
) {
  
  ## ------------------------------------------------------------
  ## Define cluster and feature order
  ## ------------------------------------------------------------
  
  cluster_order <- c(
    "Early pro-B",
    "Late pro-B",
    "Large pre-B",
    "Small pre-B",
    "Immature B",
    "Maturing B"
  )
  
  rna_genes <- c(
    "EBF1",
    "PAX5",
    "CD19",
    "MME",
    "IL7R",
    "MKI67",
    "CD34",
    "ERG",
    "DNTT",
    "RAG1",
    "RAG2",
    "IGLL1",
    "VPREB1",
    "IKZF3",
    "SPIB",
    "FCRLA",
    "TCL1A",
    "MS4A1",
    "CD83",
    "CD40",
    "FCER2",
    "CXCR5"
  )
  
  phase_order <- c(
    "S",
    "G2M",
    "G1"
  )
  
  feature_order <- c(
    "EBF1",
    "PAX5",
    "CD19",
    "MME",
    "IL7R",
    "MKI67",
    "S phase",
    "G2M phase",
    "G1 phase",
    "CD34",
    "ERG",
    "DNTT",
    "RAG1",
    "RAG2",
    "IGLL1",
    "VPREB1",
    "IKZF3",
    "SPIB",
    "FCRLA",
    "TCL1A",
    "MS4A1",
    "CD83",
    "CD40",
    "FCER2",
    "CXCR5"
  )
  
  ## ------------------------------------------------------------
  ## Check required metadata
  ## ------------------------------------------------------------
  
  required_metadata <- c(cluster_col, phase_col)
  
  missing_metadata <- setdiff(
    required_metadata,
    colnames(seurat_obj@meta.data)
  )
  
  if (length(missing_metadata) > 0) {
    stop(
      "Missing metadata column(s): ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  
  ## ------------------------------------------------------------
  ## Check RNA genes
  ## ------------------------------------------------------------
  
  assay_genes <- rownames(seurat_obj[[assay]])
  
  rna_genes_present <- intersect(
    rna_genes,
    assay_genes
  )
  
  rna_genes_missing <- setdiff(
    rna_genes,
    rna_genes_present
  )
  
  if (length(rna_genes_present) == 0) {
    stop("None of the requested RNA genes are present in the assay.")
  }
  
  if (length(rna_genes_missing) > 0) {
    message(
      "Genes missing from the ",
      assay,
      " assay: ",
      paste(rna_genes_missing, collapse = ", ")
    )
  }
  
  ## ------------------------------------------------------------
  ## Extract normalized RNA expression
  ## ------------------------------------------------------------
  
  expr_mat <- Seurat::GetAssayData(
    seurat_obj,
    assay = assay,
    layer = layer
  )[rna_genes_present, , drop = FALSE]
  
  ## Convert expression to detected/not detected
  expr_bin <- expr_mat > 0
  
  ## ------------------------------------------------------------
  ## Calculate RNA-marker percentage by cluster
  ## ------------------------------------------------------------
  
  cluster_vector <- as.character(
    seurat_obj@meta.data[[cluster_col]]
  )
  
  names(cluster_vector) <- colnames(seurat_obj)
  
  clusters_present <- cluster_order[
    cluster_order %in% unique(cluster_vector)
  ]
  
  clusters_missing <- setdiff(
    cluster_order,
    clusters_present
  )
  
  if (length(clusters_missing) > 0) {
    message(
      "Clusters not present in this object: ",
      paste(clusters_missing, collapse = ", ")
    )
  }
  
  if (length(clusters_present) == 0) {
    stop("None of the requested B-cell clusters are present.")
  }
  
  marker_pct_rna <- sapply(
    clusters_present,
    function(cl) {
      
      cells <- names(cluster_vector)[
        cluster_vector == cl &
          !is.na(cluster_vector)
      ]
      
      if (length(cells) == 0) {
        return(rep(NA_real_, length(rna_genes_present)))
      }
      
      Matrix::rowMeans(
        expr_bin[, cells, drop = FALSE]
      ) * 100
    }
  )
  
  ## Ensure matrix structure if only one cluster is present
  if (is.null(dim(marker_pct_rna))) {
    marker_pct_rna <- matrix(
      marker_pct_rna,
      ncol = 1,
      dimnames = list(
        rna_genes_present,
        clusters_present
      )
    )
  }
  
  marker_pct_rna_long <- as.data.frame(
    marker_pct_rna,
    check.names = FALSE
  ) %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(
      cols = -feature,
      names_to = "cluster",
      values_to = "percent"
    ) %>%
    dplyr::mutate(
      feature_type = "RNA marker"
    )
  
  ## ------------------------------------------------------------
  ## Calculate cell-cycle phase percentage by cluster
  ## ------------------------------------------------------------
  
  phase_pct_long <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::transmute(
      cell = cell,
      cluster = as.character(.data[[cluster_col]]),
      Phase = as.character(.data[[phase_col]])
    ) %>%
    dplyr::filter(
      cluster %in% clusters_present,
      Phase %in% phase_order
    ) %>%
    dplyr::count(
      cluster,
      Phase,
      name = "n"
    ) %>%
    tidyr::complete(
      cluster = clusters_present,
      Phase = phase_order,
      fill = list(n = 0)
    ) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      percent = dplyr::if_else(
        sum(n) > 0,
        100 * n / sum(n),
        0
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      feature = paste0(Phase, " phase"),
      feature_type = "Cell-cycle phase"
    ) %>%
    dplyr::select(
      feature,
      cluster,
      percent,
      feature_type
    )
  
  ## ------------------------------------------------------------
  ## Combine RNA-marker and cell-cycle data
  ## ------------------------------------------------------------
  
  marker_plot_df <- dplyr::bind_rows(
    marker_pct_rna_long,
    phase_pct_long
  ) %>%
    dplyr::mutate(
      feature = factor(
        feature,
        levels = feature_order
      ),
      cluster = factor(
        cluster,
        levels = cluster_order
      )
    ) %>%
    dplyr::filter(
      !is.na(feature),
      !is.na(cluster)
    )
  
  ## ------------------------------------------------------------
  ## Construct heatmap matrix
  ## ------------------------------------------------------------
  
  marker_heatmap_mat <- marker_plot_df %>%
    dplyr::select(
      feature,
      cluster,
      percent
    ) %>%
    tidyr::pivot_wider(
      names_from = cluster,
      values_from = percent,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("feature") %>%
    as.matrix()
  
  ## Add missing requested features as zero rows
  missing_feature_rows <- setdiff(
    feature_order,
    rownames(marker_heatmap_mat)
  )
  
  if (length(missing_feature_rows) > 0) {
    
    zero_rows <- matrix(
      0,
      nrow = length(missing_feature_rows),
      ncol = ncol(marker_heatmap_mat),
      dimnames = list(
        missing_feature_rows,
        colnames(marker_heatmap_mat)
      )
    )
    
    marker_heatmap_mat <- rbind(
      marker_heatmap_mat,
      zero_rows
    )
  }
  
  ## Reorder rows and columns
  marker_heatmap_mat <- marker_heatmap_mat[
    feature_order,
    clusters_present,
    drop = FALSE
  ]
  
  ## ------------------------------------------------------------
  ## Generate heatmap
  ## ------------------------------------------------------------
  
  heatmap_plot <- pheatmap::pheatmap(
    marker_heatmap_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(
      c("white", "#000040")
    )(100),
    breaks = seq(
      0,
      100,
      length.out = 101
    ),
    border_color = "grey80",
    fontsize_row = 9,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0(
      sample_name,
      " marker and cell-cycle phase percentages"
    ),
    silent = TRUE
  )
  
  ## ------------------------------------------------------------
  ## Save as transparent SVG
  ## ------------------------------------------------------------
  
  if (!is.null(output_file)) {
    
    output_dir <- dirname(
      path.expand(output_file)
    )
    
    if (!dir.exists(output_dir)) {
      dir.create(
        output_dir,
        recursive = TRUE
      )
    }
    
    svglite::svglite(
      file = path.expand(output_file),
      width = width,
      height = height,
      bg = "transparent"
    )
    
    grid::grid.newpage()
    grid::grid.draw(heatmap_plot$gtable)
    
    grDevices::dev.off()
  }
  
  ## Display plot in RStudio
  grid::grid.newpage()
  grid::grid.draw(heatmap_plot$gtable)
  
  ## Return results
  invisible(
    list(
      plot = heatmap_plot,
      heatmap_matrix = marker_heatmap_mat,
      plot_data = marker_plot_df,
      rna_genes_present = rna_genes_present,
      rna_genes_missing = rna_genes_missing,
      clusters_present = clusters_present,
      clusters_missing = clusters_missing
    )
  )
}


il7_low_heatmap <- plot_marker_phase_heatmap(
  seurat_obj = il7_low,
  sample_name = "IL7-low hPSC",
  cluster_col = "custom_labels",
  phase_col = "Phase",
  output_file = "~/Desktop/il7_low_total_b_development.svg",
  width = 3.3,
  height = 7
)

il7_high_heatmap <- plot_marker_phase_heatmap(
  seurat_obj = il7_high,
  sample_name = "IL7-high hPSC",
  cluster_col = "custom_labels",
  phase_col = "Phase",
  output_file = "~/Desktop/il7_high_total_b_development.svg",
  width = 3.3,
  height = 7
)

bm_total_heatmap <- plot_marker_phase_heatmap(
  seurat_obj = bm_total,
  sample_name = "BM total",
  cluster_col = "custom_labels",
  phase_col = "Phase",
  output_file = "~/Desktop/BM_total_b_development.svg",
  width = 3.3,
  height = 7
)
