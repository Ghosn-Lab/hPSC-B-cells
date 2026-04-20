## create_seurat_object.R
## Integrate hPSC IL7low/high and two published adult bone marrow samples

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(stringr)
library(rlang)
library(pheatmap)

## -------------------------------------------------------------------------
## File paths (update as needed)
## -------------------------------------------------------------------------

mbg7_dir    <- "data/il7_high/cellranger_filtered_output/"
mbg51_dir   <- "data/il7_low/cellranger_filtered_output/"
bm1_dir     <- "data/public_bm1/cellranger_filtered_output/"
bm2_dir     <- "data/public_bm2/cellranger_filtered_output/"
bm_ref_file <- "data/public_bm_seurat/B_v2lowmt_noreg.rds"

## -------------------------------------------------------------------------
## Helper: preprocessing + non-B-cell removal
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
  batch_label = "il7_high",
  barcode_suffix = "_MBG7"
)

mbg51 <- clean_hpsc_sample(
  data_dir = mbg51_dir,
  batch_label = "il7_low",
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
## Explore canonical B-cell development genes 
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

# Reassign DNTT+ cells within clusters 2/3/6 as Early pro-B
cells_prolif <- WhichCells(
  combined,
  idents = c("6", "2", "3"),
  expression = DNTT > 0
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
il7_high  <- subset(combined, subset = batch == "il7_high")
il7_low   <- subset(combined, subset = batch == "il7_low")

## -------------------------------------------------------------------------
## Plot UMAPs
## -------------------------------------------------------------------------

custom_colors <- c(
  "Early pro-B" = "#332288",
  "Late pro-B" = "#DDCC77",
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
  filename = "figures/for_legend_UMAP.svg",
  plot = plot_il7_low,
  device = "svg",
  bg = "transparent",
  width = 4.5,
  height = 5,
  units = "in"
)


