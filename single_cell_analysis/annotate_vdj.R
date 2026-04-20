## annotate_vdj.R
## Annotate and plod VDJ groups

library(Seurat)
library(harmony)
library(dplyr)
library(readr)
library(stringr)
library(rlang)
library(pheatmap)



## -------------------------------------------------------------------------
## Add VDJ metadata and generate HC/LC plots for each group
## -------------------------------------------------------------------------

annotate_vdj <- function(seurat_obj, vdj_df,
                         barcode_col = "cellranger_indexed_barcode",
                         chain_col   = "chain_type",
                         func_col    = "v_domain_functionality") {
  
  req <- c(barcode_col, chain_col, func_col)
  miss <- setdiff(req, colnames(vdj_df))
  if (length(miss) > 0) {
    stop("annotate_vdj: missing columns: ", paste(miss, collapse = ", "))
  }
  
  b  <- rlang::sym(barcode_col)
  ch <- rlang::sym(chain_col)
  fn <- rlang::sym(func_col)
  
  vdj <- vdj_df %>%
    dplyr::mutate(
      chain_low = tolower(!!ch),
      func_low  = tolower(!!fn),
      is_prod   = stringr::str_detect(func_low, "^prod"),
      is_unprod = stringr::str_detect(func_low, "unprod")
    )
  
  ## Heavy chain
  igh_summ <- vdj %>%
    dplyr::filter(chain_low == "igh") %>%
    dplyr::group_by(!!b) %>%
    dplyr::summarise(
      IGH_prod   = any(is_prod, na.rm = TRUE),
      IGH_unprod = any(is_unprod, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      IGH_category = dplyr::case_when(
        IGH_prod   ~ "IGH_productive",
        IGH_unprod ~ "IGH_unproductive",
        TRUE       ~ "No_IGH"
      )
    ) %>%
    dplyr::select(!!b, IGH_category)
  
  ## Light chain
  igl_summ <- vdj %>%
    dplyr::filter(chain_low == "igl") %>%
    dplyr::group_by(!!b) %>%
    dplyr::summarise(
      IGL_prod   = any(is_prod, na.rm = TRUE),
      IGL_unprod = any(is_unprod, na.rm = TRUE),
      .groups = "drop"
    )
  
  igk_summ <- vdj %>%
    dplyr::filter(chain_low == "igk") %>%
    dplyr::group_by(!!b) %>%
    dplyr::summarise(
      IGK_prod   = any(is_prod, na.rm = TRUE),
      IGK_unprod = any(is_unprod, na.rm = TRUE),
      .groups = "drop"
    )
  
  light_summ <- dplyr::full_join(igl_summ, igk_summ, by = barcode_col) %>%
    dplyr::mutate(
      has_igl = dplyr::coalesce(IGL_prod, FALSE) | dplyr::coalesce(IGL_unprod, FALSE),
      has_igk = dplyr::coalesce(IGK_prod, FALSE) | dplyr::coalesce(IGK_unprod, FALSE),
      Light_category = dplyr::case_when(
        has_igl & dplyr::coalesce(IGL_prod, FALSE)   ~ "IGL_productive",
        has_igl & dplyr::coalesce(IGL_unprod, FALSE) ~ "IGL_unproductive",
        !has_igl & has_igk & dplyr::coalesce(IGK_prod, FALSE)   ~ "IGK_productive",
        !has_igl & has_igk & dplyr::coalesce(IGK_unprod, FALSE) ~ "IGK_unproductive",
        TRUE ~ "No_Light"
      )
    ) %>%
    dplyr::select(!!b, Light_category)
  
  ## Combine labels
  vdj_labels <- dplyr::full_join(igh_summ, light_summ, by = barcode_col) %>%
    dplyr::mutate(
      IGH_category   = dplyr::coalesce(IGH_category, "No_IGH"),
      Light_category = dplyr::coalesce(Light_category, "No_Light"),
      VDJ_category   = paste0(IGH_category, " + ", Light_category)
    )
  
  ## Add default labels for cells without VDJ data
  all_cells <- data.frame(Cells = Seurat::Cells(seurat_obj), stringsAsFactors = FALSE)
  colnames(all_cells) <- barcode_col
  
  vdj_labels_all <- dplyr::left_join(all_cells, vdj_labels, by = barcode_col) %>%
    dplyr::mutate(
      IGH_category   = dplyr::coalesce(IGH_category, "No VDJ"),
      Light_category = dplyr::coalesce(Light_category, "No VDJ"),
      VDJ_category   = dplyr::coalesce(VDJ_category, "No VDJ")
    )
  
  meta_to_add <- vdj_labels_all %>%
    dplyr::select(-!!b) %>%
    as.data.frame()
  
  rownames(meta_to_add) <- vdj_labels_all %>%
    dplyr::pull(!!b)
  
  Seurat::AddMetaData(seurat_obj, metadata = meta_to_add)
}

## -------------------------------------------------------------------------
## Helper functions
## -------------------------------------------------------------------------

read_and_format_vdj <- function(file, barcode_pattern, barcode_replacement, cells_keep) {
  vdj <- readr::read_csv(file, show_col_types = FALSE)
  vdj$cellranger_indexed_barcode <- gsub(
    barcode_pattern,
    barcode_replacement,
    vdj$cellranger_indexed_barcode
  )
  vdj <- vdj[vdj$cellranger_indexed_barcode %in% cells_keep, ]
  vdj
}

plot_hc <- function(obj, file_out) {
  plot <- DimPlot(
    obj,
    group.by = "IGH_category",
    reduction = "umap_flipped",
    pt.size = 0.75,
    cols = c(
      "No_IGH"           = "grey80",
      "No VDJ"           = "grey80",
      "IGH_productive"   = "springgreen4",
      "IGH_unproductive" = "grey80"
    ),
    order = c("IGH_productive", "IGH_unproductive", "No_IGH", "No VDJ"),
    shuffle = FALSE
  ) + NoLegend()
  
  ggsave(
    filename = file_out,
    plot = plot,
    device = "svg",
    bg = "transparent",
    width = 4.5,
    height = 5.5,
    units = "in"
  )
}

plot_lc <- function(obj, file_out) {
  plot <- DimPlot(
    obj,
    group.by = "Light_category",
    reduction = "umap_flipped",
    pt.size = 0.75,
    cols = c(
      "No_Light"         = "grey80",
      "No VDJ"           = "grey80",
      "IGL_unproductive" = "grey80",
      "IGK_unproductive" = "grey80",
      "IGK_productive"   = "purple3",
      "IGL_productive"   = "darkorange3"
    ),
    order = rev(c(
      "No_Light",
      "No VDJ",
      "IGK_unproductive",
      "IGL_unproductive",
      "IGL_productive",
      "IGK_productive"
    )),
    shuffle = FALSE
  ) + NoLegend()
  
  ggsave(
    filename = file_out,
    plot = plot,
    device = "svg",
    bg = "transparent",
    width = 4.5,
    height = 5.5,
    units = "in"
  )
}

## -------------------------------------------------------------------------
## Create group-specific subset objects from combined
## -------------------------------------------------------------------------

bm_total <- subset(combined, subset = batch == "BM_total")
il7_high <- subset(combined, subset = batch == "il7_high")
il7_low  <- subset(combined, subset = batch == "il7_low")

## -------------------------------------------------------------------------
## Read VDJ data and annotate each group
## Replace file names with actual project-relative paths as needed
## -------------------------------------------------------------------------

vdj_il7_low <- read_and_format_vdj(
  file = "data/vdj/il7_low/imgt_processed_output.csv",
  barcode_pattern = "-170$",
  barcode_replacement = "_MBG51",
  cells_keep = colnames(il7_low)
)
il7_low <- annotate_vdj(il7_low, vdj_il7_low)

vdj_il7_high <- read_and_format_vdj(
  file = "data/vdj/il7_high/imgt_processed_output.csv",
  barcode_pattern = "-158$",
  barcode_replacement = "_MBG7",
  cells_keep = colnames(il7_high)
)
il7_high <- annotate_vdj(il7_high, vdj_il7_high)

vdj_bm1 <- read_and_format_vdj(
  file = "data/vdj/bm1/imgt_processed_output.csv",
  barcode_pattern = "-57$",
  barcode_replacement = "_BM1",
  cells_keep = grep("_BM1$", colnames(bm_total), value = TRUE)
)

vdj_bm2 <- read_and_format_vdj(
  file = "data/vdj/bm2/imgt_processed_output.csv",
  barcode_pattern = "-58$",
  barcode_replacement = "_BM2",
  cells_keep = grep("_BM2$", colnames(bm_total), value = TRUE)
)

vdj_bm_total <- dplyr::bind_rows(vdj_bm1, vdj_bm2)
bm_total <- annotate_vdj(bm_total, vdj_bm_total)

## -------------------------------------------------------------------------
## Add VDJ metadata back to combined
## -------------------------------------------------------------------------

meta_cols <- c("IGH_category", "Light_category", "VDJ_category")

get_vdj_meta <- function(obj, cols) {
  md <- obj@meta.data[, cols, drop = FALSE]
  md <- md[!duplicated(rownames(md)), , drop = FALSE]
  md
}

vdj_meta_all <- dplyr::bind_rows(
  get_vdj_meta(il7_low, meta_cols),
  get_vdj_meta(il7_high, meta_cols),
  get_vdj_meta(bm_total, meta_cols)
)

combined <- AddMetaData(combined, metadata = vdj_meta_all)

## -------------------------------------------------------------------------
## Save HC and LC plots for each group
## -------------------------------------------------------------------------

plot_hc(bm_total,  "figures/BM_total_HC.svg")
plot_lc(bm_total,  "figures/BM_total_LC.svg")

plot_hc(il7_high,  "figures/il7_high_HC.svg")
plot_lc(il7_high,  "figures/il7_high_LC.svg")

plot_hc(il7_low,   "figures/il7_low_HC.svg")
plot_lc(il7_low,   "figures/il7_low_LC.svg")