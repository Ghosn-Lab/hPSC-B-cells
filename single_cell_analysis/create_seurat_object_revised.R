## create_seurat_object.R
## Integrate hPSC IL7low/high and two published adult bone marrow samples
## Explore B-cell genes to determine trajectory
## Make cell-cycle phase plots, perform GO analysis for Mature B-cells
## Additional filtering based on V(D)J shown in annotate_vdj.R

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(stringr)
library(rlang)
library(pheatmap)

## -------------------------------------------------------------------------
## Load files
## -------------------------------------------------------------------------

mbg7_dir    <- "data/il7_high/cellranger_filtered_output/"
mbg51_dir   <- "data/il7_low/cellranger_filtered_output/"
bm1_dir     <- "data/public_bm1/cellranger_filtered_output/"
bm2_dir     <- "data/public_bm2/cellranger_filtered_output/"
bm_ref_file <- "data/public_bm_seurat/B_v2lowmt_noreg.rds"

## -------------------------------------------------------------------------
## Function for cleaning + non-B-cell removal
## -------------------------------------------------------------------------

clean_hpsc_sample <- function(data_dir, batch_label, barcode_suffix) {
  obj <- Read10X(data.dir = data_dir) |>
    CreateSeuratObject()
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- NormalizeData(obj)
  obj <- subset(obj, subset = percent.mt < 20)
  
  non_b_cells <- WhichCells(
    obj,
    expression = CD4 > 0 | LYZ > 0 | MPO > 0 | CDH5 > 0 |
      CD3E > 0 | ITGAM > 0 | ITGAX > 0 | CDH11 > 0
  )
  
  obj <- subset(obj, cells = setdiff(Cells(obj), non_b_cells))
  obj$batch <- batch_label
  
  obj <- RenameCells(
    obj,
    new.names = gsub("-1$", barcode_suffix, colnames(obj))
  )
  
  return(obj)
}

## -------------------------------------------------------------------------
## Load and clean mbg datasets
## -------------------------------------------------------------------------

mbg7 <- clean_hpsc_sample(
  data_dir = mbg7_dir,
  batch_label = " IL7_high_hPSC",
  barcode_suffix = "_MBG7"
)

mbg51 <- clean_hpsc_sample(
  data_dir = mbg51_dir,
  batch_label = "IL7_low_hPSC",
  barcode_suffix = "_MBG51"
)

## -------------------------------------------------------------------------
## Load and clean published bone marrow samples (GSE181543)
## -------------------------------------------------------------------------

# Published filtered barcode reference
bm_barcodes <- readRDS(bm_ref_file)
published_barcodes <- colnames(bm_barcodes)

bm1_bcs <- sub("_2$", "_BM1", published_barcodes)
bm2_bcs <- sub("_3$", "_BM2", published_barcodes)

# Load BM datasets
bm1 <- Read10X(data.dir = bm1_dir) |>
  CreateSeuratObject()
bm1 <- RenameCells(bm1, new.names = gsub("-1$", "_BM1", colnames(bm1)))
bm1$batch <- "BM1"

bm2 <- Read10X(data.dir = bm2_dir) |>
  CreateSeuratObject()
bm2 <- RenameCells(bm2, new.names = gsub("-1$", "_BM2", colnames(bm2)))
bm2$batch <- "BM2"

# Restrict to published clean barcodes
keep_bm1 <- intersect(Cells(bm1), bm1_bcs)
keep_bm2 <- intersect(Cells(bm2), bm2_bcs)

bm1 <- subset(bm1, cells = keep_bm1)
bm2 <- subset(bm2, cells = keep_bm2)

bm_total <- merge(bm1, bm2)
bm_total <- subset(bm_total, cells = published_barcodes)



## -------------------------------------------------------------------------
## Merge datasets and run Harmony integration
## -------------------------------------------------------------------------

combined <- merge(bm1, y = list(bm2, mbg7, mbg51))
combined <- JoinLayers(combined)
combined <- NormalizeData(combined)

# Remove non-B cells, then retain CD79A+ cells
non_b_cells <- WhichCells(
  combined,
  expression = CD4 > 0 | LYZ > 0 | MPO > 0 | CDH5 > 0 | CD3E > 0 |
    ITGAM > 0 | ITGAX > 0 | CDH11 > 0 | TPSAB1 > 0 | CSF1R > 0 | MMP9 > 0
)
combined <- subset(combined, cells = setdiff(Cells(combined), non_b_cells))

b_cells <- WhichCells(combined, expression = CD79A > 0)
combined <- subset(combined, cells = b_cells)

# Identify variable features and remove IG and cell cycle genes
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)

var_genes <- VariableFeatures(combined)
ig_pattern <- "^IG[HKL]"
cell_cycle_genes <- unique(c(cc.genes$s.genes, cc.genes$g2m.genes))

filtered_var_genes <- var_genes[!grepl(ig_pattern, var_genes)]
filtered_var_genes <- setdiff(filtered_var_genes, cell_cycle_genes)
VariableFeatures(combined) <- filtered_var_genes

combined <- ScaleData(combined, features = VariableFeatures(combined))
combined <- RunPCA(combined, features = VariableFeatures(combined))
combined <- RunHarmony(combined, group.by.vars = "batch")

# Cluster and visualize using Harmony embeddings
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

# DimPlot(combined, reduction = "umap", group.by = "batch", label = TRUE)
DimPlot(combined, reduction = "umap", split.by = "batch")

## -------------------------------------------------------------------------
## Create flipped UMAP, for left to right B-cell development visualization
## -------------------------------------------------------------------------

umap_coords <- Embeddings(combined, reduction = "umap")
umap_flipped <- umap_coords
umap_flipped[, 1] <- -1 * umap_flipped[, 1]
umap_flipped[, 2] <- -1 * umap_flipped[, 2]

combined[["umap_flipped"]] <- CreateDimReducObject(
  embeddings = umap_flipped,
  key = "UMAPF_",
  assay = DefaultAssay(combined)
)

## -------------------------------------------------------------------------
## Combine both BM samples into one
## -------------------------------------------------------------------------

combined$batch <- recode(
  combined$batch,
  "BM1" = "BM_total",
  "BM2" = "BM_total"
)


## -------------------------------------------------------------------------
## Explore B-cell development genes 
## -------------------------------------------------------------------------
FeaturePlot(combined, c("IGHD", "IGHM", "IGHG2", "IGHA1", "IGHG1", "IGHE",
                        "IGHG3", "IGHA2", "NR4A1", "SELL", "CR2",
                        "SLAMF7", "TNFSF12"), order = T)
DotPlot(combined, c("IGHD", "IGHM", "IGHG2", "IGHA1", "NR4A1", "SELL",
                    "CR2", "SLAMF7", "TNFSF12"))
FeaturePlot(combined, c("EBF1", "CD7", "MME", "DNTT", "ERG", "MKI67",
                        "RAG1", "IGF2", "MS4A1"),
            order = T, reduction = "umap_flipped")
FeaturePlot(combined, c("CD9", "TCL1A", "IKZF3", "FCRLA"),
            order = T, reduction = "umap_flipped", split.by = "batch")
DimPlot(combined, reduction = "umap_flipped", split.by = "batch")
VlnPlot(combined, c("MKI67", "IL7R", "RAG1", "RAG2"),
        split.by = "batch", pt.size = 0)



## -------------------------------------------------------------------------
## Save object
## -------------------------------------------------------------------------
write_rds(combined, "data/processed/hPSC_IL7low_IL7high_BM.rds")


## -------------------------------------------------------------------------
## Assign cluster labels
## -------------------------------------------------------------------------

Idents(combined) <- combined$seurat_clusters
combined$custom_labels <- as.character(Idents(combined))

combined$custom_labels[combined$custom_labels %in% c("6", "2", "3")] <- "Large pre-B"
combined$custom_labels[combined$custom_labels == "5"] <- "Late pro-B"
combined$custom_labels[combined$custom_labels %in% c("0", "4")] <- "Small pre-B"
combined$custom_labels[combined$custom_labels == "7"] <- "Immature B"
combined$custom_labels[combined$custom_labels == "1"] <- "Maturing B"
combined$custom_labels[combined$custom_labels == "9"] <- "Proliferating Mature B"
combined$custom_labels[combined$custom_labels == "8"] <- "Metabolic stress"

# Reassign ERG+ cells within clusters 2/3/6 (Proliferating cells) as Early pro-B 
cells_prolif <- WhichCells(
  combined,
  idents = c("6", "2", "3"),
  expression = ERG > 0
)
combined$custom_labels[cells_prolif] <- "Early pro-B"

new_levels <- c(
  "Early pro-B",
  "Late pro-B",
  "Large pre-B",
  "Small pre-B",
  "Immature B",
  "Maturing B",
  "Proliferating Mature B",
  "Metabolic stress"
)

combined$custom_labels <- factor(combined$custom_labels, levels = new_levels)
Idents(combined) <- combined$custom_labels

## -------------------------------------------------------------------------
## Create subset objects for plotting
## -------------------------------------------------------------------------

bm_total  <- subset(combined, subset = batch == "BM_total")
il7_high  <- subset(combined, subset = batch == "IL7_high_hPSC")
il7_low   <- subset(combined, subset = batch == "IL7_low_hPSC")

## -------------------------------------------------------------------------
## Plot UMAPs
## -------------------------------------------------------------------------

custom_colors <- c(
  "Early pro-B" = "#DDCC77",
  "Late pro-B" = "#332288",
  "Large pre-B" = "#999933",
  "Small pre-B" = "#44AA99",
  "Immature B" = "#117733",
  "Maturing B" = "#CC6677",
  "Proliferating Mature B" = "#882255",
  "Metabolic stress" = "#AA4499"
)

plot_bm_total <- DimPlot(
  bm_total,
  reduction = "umap_flipped",
  cols = custom_colors,
  pt.size = 1
) + NoLegend()

plot_il7_high <- DimPlot(
  il7_high,
  reduction = "umap_flipped",
  cols = custom_colors,
  pt.size = 1
) + NoLegend()

plot_il7_low <- DimPlot(
  il7_low,
  reduction = "umap_flipped",
  cols = custom_colors,
  pt.size = 1
) + NoLegend()

ggsave(
  filename = "~/Desktop/bm_total_umap.svg",
  plot = plot_bm_total,
  device = "svg",
  bg = "transparent",
  width = 4.5,
  height = 5,
  units = "in"
)





## -------------------------------------------------------------------------
## Assign sample origin metadata and make new dimplots for supplement
## -------------------------------------------------------------------------
## Extract sample suffix from cell barcode
combined$origin_sample <- sub(".*_", "", colnames(combined))
## Check result
table(combined$origin_sample)
table(
  batch = combined$batch,
  origin_sample = combined$origin_sample
)

bm1   <- subset(combined, subset = origin_sample == "BM1")
bm2   <- subset(combined, subset = origin_sample == "BM2")
mbg7  <- subset(combined, subset = origin_sample == "MBG7")
mbg51 <- subset(combined, subset = origin_sample == "MBG51")

custom_colors <- c(
  "Early pro-B" = "#DDCC77",
  "Late pro-B" = "#332288",
  "Large pre-B" = "#999933",
  "Small pre-B" = "#44AA99",
  "Immature B" = "#117733",
  "Maturing B" = "#CC6677",
  "Proliferating Mature B" = "#882255",
  "Metabolic stress" = "grey"
)

plot_bm1 <- DimPlot(
  bm1,
  reduction = "umap_flipped",
  cols = custom_colors,
  pt.size = 1
) + NoLegend()

plot_bm2 <- DimPlot(
  bm2,
  reduction = "umap_flipped",
  cols = custom_colors,
  pt.size = 1
) + NoLegend()

plot_mbg7 <- DimPlot(
  mbg7,
  reduction = "umap_flipped",
  cols = custom_colors,
  pt.size = 1
) + NoLegend()

plot_mbg51 <- DimPlot(
  mbg51,
  reduction = "umap_flipped",
  cols = custom_colors,
  pt.size = 1
) + NoLegend()



ggsave(
  filename = "~/Desktop/BM1_UMAP.svg",
  plot = plot_bm1,
  device = "svg",
  bg = "transparent",
  width = 4.5,
  height = 5,
  units = "in"
)

ggsave(
  filename = "~/Desktop/BM2_UMAP.svg",
  plot = plot_bm2,
  device = "svg",
  bg = "transparent",
  width = 4.5,
  height = 5,
  units = "in"
)

ggsave(
  filename = "~/Desktop/IL7_high_UMAP.svg",
  plot = plot_mbg7,
  device = "svg",
  bg = "transparent",
  width = 4.5,
  height = 5,
  units = "in"
)

ggsave(
  filename = "~/Desktop/IL7_low_UMAP.svg",
  plot = plot_mbg51,
  device = "svg",
  bg = "transparent",
  width = 4.5,
  height = 5,
  units = "in"
)



## -------------------------------------------------------------------------
## Assign cell cycle scores to combined
## -------------------------------------------------------------------------

# Seurat has built-in human cell cycle gene sets
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes

# Keep only genes present in your object
s.genes.use <- intersect(s.genes, rownames(combined))
g2m.genes.use <- intersect(g2m.genes, rownames(combined))

length(s.genes.use)
length(g2m.genes.use)

combined <- CellCycleScoring(
  combined,
  s.features = s.genes.use,
  g2m.features = g2m.genes.use,
  set.ident = FALSE
)

# Check assigned phases
table(combined$Phase)

## -------------------------------------------------------------------------
## Plot cell cycle phase on UMAP
## -------------------------------------------------------------------------

plot_cell_cycle <- DimPlot(
  combined,
  reduction = "umap_flipped",
  group.by = "Phase",
  pt.size = 0.7
)

plot_cell_cycle

ggsave(
  filename = "~/Desktop/cell_cycle_UMAP.svg",
  plot = plot_cell_cycle,
  device = "svg",
  bg = "transparent",
  width = 5.5,
  height = 6,
  units = "in"
)



## -------------------------------------------------------------------------
## FeaturePlots: ERG and SPIB
## -------------------------------------------------------------------------

plot_erg <- FeaturePlot(
  combined,
  features = "ERG",
  reduction = "umap_flipped",
  order = T,
  pt.size = 0.7
)

plot_spib <- FeaturePlot(
  combined,
  features = "SPIB",
  reduction = "umap_flipped",
  order = T,
  pt.size = 0.7
)

plot_erg
plot_spib

ggsave(
  filename = "~/Desktop/combined_ERG_UMAP.svg",
  plot = plot_erg,
  device = "svg",
  bg = "transparent",
  width = 6,
  height = 5.5,
  units = "in"
)

ggsave(
  filename = "~/Desktop/combined_SPIB_UMAP.svg",
  plot = plot_spib,
  device = "svg",
  bg = "transparent",
  width = 6,
  height = 5.5,
  units = "in"
)




## -------------------------------------------------------------------------
## Compare mature ABM to mature hPSC (IL7 low)
## Top 10 gene DEG
## Pseudobulk expression comparison for Mature B cells
## -------------------------------------------------------------------------

bm_total  <- subset(combined, subset = batch == "BM_total")
mbg51 <- subset(combined, subset = origin_sample == "MBG51")

DimPlot(bm_total)

deg_mature_abm <- FindMarkers(bm_total, ident.1 = "Maturing B")

mature_b <- subset(
  combined,
  subset = custom_labels == "Maturing B"
)

mature_b <- subset(
  mature_b,
  subset = batch %in% c("BM_total", "IL7_low_hPSC")
)


top50_genes <- deg_mature_abm %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(
    !stringr::str_detect(gene, "^RP[SL]"),
    !stringr::str_detect(gene, "^MT-")
  ) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::pull(gene)

top50_genes

plot_mature_dotplot <- DotPlot(
  mature_b,
  features = rev(top50_genes),
  group.by = "batch",
  assay = "RNA",
  scale = FALSE
) +
  scale_color_gradient2(
    low = "navy",
    mid = "lightblue",
    high = "firebrick",
    midpoint = 3,
    name = "Average\nexpression"
  ) +
  coord_flip() +
  RotatedAxis() +
  theme_classic() +
  labs(
    title = "Top 50 Mature B DEG genes by batch",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(hjust = 0.5)
  )

plot_mature_dotplot



## -------------------------------------------------------------------------
## Gene ontology analysis
## -------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

table(Idents(combined))
table(combined$batch)

samples_use <- c(
  "BM_total",
  "IL7_high_hPSC",
  "IL7_low_hPSC"
)

stages_use <- c(
  "Early pro-B",
  "Late pro-B",
  "Large pre-B",
  "Small pre-B",
  "Immature B",
  "Maturing B"
)

bm_sample <- "BM_total"

## Make sure custom_labels has your original B-cell stage labels
table(combined$custom_labels)

combined$GO_stage_group <- dplyr::case_when(
  combined$custom_labels %in% c("Early pro-B", "Large pre-B") ~ "Proliferation",
  combined$custom_labels %in% c("Late pro-B", "Small pre-B") ~ "Rearrangement",
  combined$custom_labels %in% c("Immature B", "Maturing B") ~ "Maturation",
  TRUE ~ NA_character_
)

combined$GO_stage_group <- factor(
  combined$GO_stage_group,
  levels = c("Proliferation", "Rearrangement", "Maturation")
)

table(combined$GO_stage_group, useNA = "ifany")
table(combined$batch, combined$GO_stage_group)

Idents(combined) <- "GO_stage_group"

levels(Idents(combined))
table(Idents(combined))

bm_sample <- "BM_total"

bm_grouped_stages <- c(
  "Proliferation",
  "Rearrangement",
  "Maturation"
)

bm_go_results_grouped <- list()

for (stage_name in bm_grouped_stages) {
  
  result_name <- paste(bm_sample, stage_name, sep = "_")
  
  bm_go_results_grouped[[result_name]] <- run_go_for_stage(
    seurat_obj = combined,
    sample_name = bm_sample,
    stage_name = stage_name,
    sample_col = "batch"
  )
}

bm_go_summary_grouped <- tibble::tibble(
  result_name = names(bm_go_results_grouped),
  n_markers = sapply(bm_go_results_grouped, function(x) {
    if (is.null(x)) return(NA_integer_)
    nrow(x$markers)
  }),
  n_go_terms = sapply(bm_go_results_grouped, function(x) {
    if (is.null(x) || is.null(x$go)) return(NA_integer_)
    nrow(as.data.frame(x$go))
  })
)

bm_go_summary_grouped

bm_go_df_grouped <- purrr::map_dfr(names(bm_go_results_grouped), function(result_name) {
  
  result <- bm_go_results_grouped[[result_name]]
  
  if (is.null(result)) return(NULL)
  if (is.null(result$go)) return(NULL)
  
  go_df <- as.data.frame(result$go)
  
  if (nrow(go_df) == 0) return(NULL)
  
  go_df %>%
    dplyr::mutate(
      sample = result$sample,
      GO_stage_group = result$stage,
      result_name = result_name,
      neg_log10_padj = -log10(p.adjust)
    )
})

dim(bm_go_df_grouped)
head(bm_go_df_grouped)

bm_go_top10_grouped <- bm_go_df_grouped %>%
  dplyr::group_by(GO_stage_group) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    GO_stage_group = factor(
      GO_stage_group,
      levels = c("Proliferation", "Rearrangement", "Maturation")
    ),
    Description_clean = stringr::str_wrap(Description, width = 45)
  )

go_order <- bm_go_top10_grouped %>%
  dplyr::group_by(Description_clean) %>%
  dplyr::summarise(
    max_sig = max(neg_log10_padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(max_sig) %>%
  dplyr::pull(Description_clean)

bm_go_top10_grouped <- bm_go_top10_grouped %>%
  dplyr::mutate(
    Description_clean = factor(Description_clean, levels = go_order)
  )

plot_bm_go_single_group <- function(
    bm_go_df,
    group_name,
    top_n = 15,
    save_name = NULL
) {
  
  plot_df <- bm_go_df %>%
    dplyr::filter(GO_stage_group == group_name) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      neg_log10_padj = -log10(p.adjust),
      Description_clean = stringr::str_wrap(Description, width = 45),
      Description_clean = forcats::fct_reorder(
        Description_clean,
        neg_log10_padj
      )
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = neg_log10_padj,
      y = Description_clean,
      size = Count,
      color = neg_log10_padj
    )
  ) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(
      low = "blue",
      high = "red",
      name = "-log10(adj. p)"
    ) +
    scale_size(
      range = c(2, 7),
      name = "Gene count"
    ) +
    theme_classic() +
    labs(
      title = paste0("BM_total GO BP: ", group_name),
      x = "-log10 adjusted p-value",
      y = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
  
  if (!is.null(save_name)) {
    ggsave(
      filename = paste0("~/Desktop/", save_name, ".svg"),
      plot = p,
      device = "svg",
      bg = "transparent",
      width = 6.2,
      height = 6.5,
      units = "in"
    )
  }
  
  return(p)
}


plot_bm_go_proliferation <- plot_bm_go_single_group(
  bm_go_df = bm_go_df_grouped,
  group_name = "Proliferation",
  top_n = 10,
  save_name = "BM_total_GO_Proliferation"
)

plot_bm_go_proliferation

plot_bm_go_rearrangement <- plot_bm_go_single_group(
  bm_go_df = bm_go_df_grouped,
  group_name = "Rearrangement",
  top_n = 10,
  save_name = "BM_total_GO_Rearrangement"
)

plot_bm_go_rearrangement

plot_bm_go_maturation <- plot_bm_go_single_group(
  bm_go_df = bm_go_df_grouped,
  group_name = "Maturation",
  top_n = 15,
  save_name = "BM_total_GO_Maturation"
)

plot_bm_go_maturation


bm_reference_go_terms_grouped <- bm_go_df_grouped %>%
  dplyr::group_by(GO_stage_group) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    GO_stage_group,
    ID,
    Description,
    p.adjust,
    Count,
    GeneRatio,
    BgRatio
  ) %>%
  dplyr::distinct()

bm_reference_go_terms_grouped

write.csv(
  bm_reference_go_terms_grouped,
  "~/Desktop/BM_total_reference_GO_terms_grouped_top15.csv",
  row.names = FALSE
)



il7low_sample <- "IL7_low_hPSC"

il7low_grouped_stages <- c(
  "Proliferation",
  "Rearrangement",
  "Maturation"
)

il7low_go_results_grouped <- list()

for (stage_name in il7low_grouped_stages) {
  
  result_name <- paste(il7low_sample, stage_name, sep = "_")
  
  il7low_go_results_grouped[[result_name]] <- run_go_for_stage(
    seurat_obj = combined,
    sample_name = il7low_sample,
    stage_name = stage_name,
    sample_col = "batch"
  )
}

il7low_go_summary_grouped <- tibble::tibble(
  result_name = names(il7low_go_results_grouped),
  n_markers = sapply(il7low_go_results_grouped, function(x) {
    if (is.null(x)) return(NA_integer_)
    nrow(x$markers)
  }),
  n_go_terms = sapply(il7low_go_results_grouped, function(x) {
    if (is.null(x) || is.null(x$go)) return(NA_integer_)
    nrow(as.data.frame(x$go))
  })
)

il7low_go_summary_grouped

il7low_go_df_grouped <- purrr::map_dfr(names(il7low_go_results_grouped), function(result_name) {
  
  result <- il7low_go_results_grouped[[result_name]]
  
  if (is.null(result)) return(NULL)
  if (is.null(result$go)) return(NULL)
  
  go_df <- as.data.frame(result$go)
  
  if (nrow(go_df) == 0) return(NULL)
  
  go_df %>%
    dplyr::mutate(
      sample = result$sample,
      GO_stage_group = result$stage,
      result_name = result_name,
      neg_log10_padj = -log10(p.adjust)
    )
})

dim(il7low_go_df_grouped)
head(il7low_go_df_grouped)

il7low_projected_go <- bm_reference_go_terms_grouped %>%
  dplyr::select(
    GO_stage_group,
    ID,
    BM_Description = Description
  ) %>%
  dplyr::left_join(
    il7low_go_df_grouped %>%
      dplyr::select(
        GO_stage_group,
        ID,
        Description,
        p.adjust,
        Count,
        GeneRatio,
        BgRatio,
        neg_log10_padj
      ),
    by = c("GO_stage_group", "ID")
  ) %>%
  dplyr::mutate(
    Description = dplyr::coalesce(Description, BM_Description),
    neg_log10_padj = dplyr::if_else(
      is.na(neg_log10_padj),
      0,
      neg_log10_padj
    ),
    Count = dplyr::if_else(
      is.na(Count),
      0L,
      as.integer(Count)
    )
  )

head(il7low_projected_go)

plot_projected_go_single_group <- function(
    projected_go_df,
    group_name,
    sample_label = "IL7_low_hPSC",
    save_name = NULL
) {
  
  plot_df <- projected_go_df %>%
    dplyr::filter(GO_stage_group == group_name) %>%
    dplyr::mutate(
      Description_clean = stringr::str_wrap(Description, width = 45),
      Description_clean = forcats::fct_reorder(
        Description_clean,
        neg_log10_padj
      )
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = neg_log10_padj,
      y = Description_clean,
      size = Count,
      color = neg_log10_padj
    )
  ) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(
      low = "blue",
      high = "red",
      name = "-log10(adj. p)"
    ) +
    scale_size(
      range = c(1, 7),
      name = "Gene count"
    ) +
    theme_classic() +
    labs(
      title = paste0(sample_label, " projected GO BP: ", group_name),
      x = "-log10 adjusted p-value",
      y = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
  
  if (!is.null(save_name)) {
    ggsave(
      filename = paste0("~/Desktop/", save_name, ".svg"),
      plot = p,
      device = "svg",
      bg = "transparent",
      width = 6.2,
      height = 6.5,
      units = "in"
    )
  }
  
  return(p)
}

plot_il7low_go_proliferation <- plot_projected_go_single_group(
  projected_go_df = il7low_projected_go,
  group_name = "Proliferation",
  sample_label = "IL7_low_hPSC",
  save_name = "IL7_low_hPSC_projected_GO_Proliferation"
)

plot_il7low_go_proliferation

plot_il7low_go_rearrangement <- plot_projected_go_single_group(
  projected_go_df = il7low_projected_go,
  group_name = "Rearrangement",
  sample_label = "IL7_low_hPSC",
  save_name = "IL7_low_hPSC_projected_GO_Rearrangement"
)

plot_il7low_go_rearrangement

plot_il7low_go_maturation <- plot_projected_go_single_group(
  projected_go_df = il7low_projected_go,
  group_name = "Maturation",
  sample_label = "IL7_low_hPSC",
  save_name = "IL7_low_hPSC_projected_GO_Maturation"
)

plot_il7low_go_maturation






## -------------------------------------------------------------------------
## Reconstruct shared BM_total + IL7_low_hPSC GO dataframe
## using all shared GO terms
## -------------------------------------------------------------------------

plot_selected_go_group <- function(
    shared_df,
    group_name,
    go_ids,
    save_name = NULL
) {
  
  id_order_df <- tibble::tibble(
    ID = go_ids,
    go_order = seq_along(go_ids)
  )
  
  plot_df <- shared_df %>%
    dplyr::filter(
      GO_stage_group == group_name,
      ID %in% go_ids
    ) %>%
    dplyr::left_join(id_order_df, by = "ID") %>%
    dplyr::mutate(
      Description_clean = paste0(
        stringr::str_wrap(Description, width = 45),
        " [", ID, "]"
      )
    )
  
  ## Define y-axis order using your supplied GO ID order
  go_order <- plot_df %>%
    dplyr::filter(sample == "BM_total") %>%
    dplyr::arrange(go_order) %>%
    dplyr::distinct(ID, Description_clean, go_order) %>%
    dplyr::pull(Description_clean)
  
  ## Reverse so first GO ID appears at top
  go_order <- rev(go_order)
  
  # plot_df <- plot_df %>%
  #   dplyr::mutate(
  #     Description_clean = factor(
  #       Description_clean,
  #       levels = go_order
  #     )
  #   )
  
  ## Least significant adjusted p-value among plotted terms
  # This is to show the minimum p-value for our plot
  max_padj_bm <- plot_df %>%
    dplyr::filter(sample == "BM_total") %>%
    dplyr::summarise(max_padj = max(p.adjust, na.rm = TRUE)) %>%
    dplyr::pull(max_padj)
  
  max_padj_il7low <- plot_df %>%
    dplyr::filter(sample == "IL7_low_hPSC") %>%
    dplyr::summarise(max_padj = max(p.adjust, na.rm = TRUE)) %>%
    dplyr::pull(max_padj)
  
  caption_text <- paste0(
    "Least significant adjusted p-value among plotted terms: ",
    "BM_total = ", signif(max_padj_bm, 3),
    "; IL7_low_hPSC = ", signif(max_padj_il7low, 3)
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = Description_clean,
      y = sample,
      size = Count
    )
  ) +
    geom_point(color = "black", alpha = 0.85) +
    scale_size(
      range = c(2, 7),
      name = "Gene count"
    ) +
    theme_classic() +
    labs(
      title = paste0("Shared GO BP enrichment: ", group_name),
      x = NULL,
      y = NULL,
      caption = caption_text
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  
  if (!is.null(save_name)) {
    ggsave(
      filename = paste0("~/Desktop/", save_name, ".svg"),
      plot = p,
      device = "svg",
      bg = "transparent",
      width = 8,
      height = 4,
      units = "in"
    )
  }
  
  return(p)
}

shared_go_plot_df_full <- dplyr::bind_rows(
  shared_go_terms_full %>%
    dplyr::transmute(
      sample = "BM_total",
      GO_stage_group,
      ID,
      Description = BM_Description,
      p.adjust = BM_padj,
      Count = BM_Count,
      neg_log10_padj = BM_neg_log10_padj
    ),
  shared_go_terms_full %>%
    dplyr::transmute(
      sample = "IL7_low_hPSC",
      GO_stage_group,
      ID,
      Description = BM_Description,  # use BM labels for matching y-axis
      p.adjust = hPSC_padj,
      Count = hPSC_Count,
      neg_log10_padj = hPSC_neg_log10_padj
    )
) %>%
  dplyr::mutate(
    sample = factor(sample, levels = c("BM_total", "IL7_low_hPSC")),
    GO_stage_group = factor(
      GO_stage_group,
      levels = c("Proliferation", "Rearrangement", "Maturation")
    )
  )

shared_go_plot_df_full %>%
  dplyr::count(sample, GO_stage_group)

maturation_go_ids <- c(
  "GO:0019886",
  "GO:0042113",
  "GO:0002449",
  "GO:0002460",
  "GO:0016064",
  "GO:0019724",
  "GO:0002377",
  "GO:0050864",
  "GO:0050853",
  "GO:0030183"
)

plot_shared_maturation_selected <- plot_selected_go_group(
  shared_df = shared_go_plot_df_full,
  group_name = "Maturation",
  go_ids = maturation_go_ids,
  save_name = "Shared_GO_BMtotal_IL7low_Maturation_selected_black"
)

plot_shared_maturation_selected


