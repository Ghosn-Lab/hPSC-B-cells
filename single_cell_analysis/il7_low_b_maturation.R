## IL7_low_B_maturation.R
## Dotplot and heatmap for IL-7 low condition hPSC --> B-cell matuaration

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(rlang)
library(purrr)
library(tibble)
library(grid)

## -------------------------------------------------------------------------
## B-cell maturation in IL-7 low hPSC
## -------------------------------------------------------------------------

mbg51_sub <- subset(combined, subset = batch == "il7_low")

plot <- DotPlot(
  mbg51_sub,
  features = c("CCR6", "CR2", "CCR10", "SLAMF7", "IGHG3"),
  idents = c(
    "Early pro-B",
    "Late pro-B",
    "Large pre-B",
    "Small pre-B",
    "Immature B",
    "Maturing B",
    "Proliferating Mature B"
  )
) +
  scale_size(range = c(-3, 8)) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.15, "cm")
  )

ggsave(
  filename = "figures/MBG51_dotplot_adjusted.svg",
  plot = plot,
  device = "svg",
  bg = "transparent",
  width = 3,
  height = 2.5,
  units = "in"
)

## -------------------------------------------------------------------------
## Heatmap for ADT stored in meta.data
## -------------------------------------------------------------------------
make_adt_heatmap_meta <- function(
    obj,
    features,
    group.by = "seurat_clusters", 
    caps = NULL,                    
    value_scale = c("feature_z","global_minmax","none"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    palette = c("black","white"),
    fontsize = 10
) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(pheatmap); library(rlang); library(purrr)
  })
  value_scale <- match.arg(value_scale)
  
  # choose grouping column
  if (group.by == "ident") {
    grp <- Idents(obj)
    grp_name <- "active_ident"
    df_meta <- obj@meta.data
    df_meta[[grp_name]] <- grp
  } else {
    grp_name <- group.by
    if (!grp_name %in% colnames(obj@meta.data))
      stop("group.by '", grp_name, "' not found in obj@meta.data")
    df_meta <- obj@meta.data
  }
  
  # check features exist in metadata
  missing <- setdiff(features, colnames(df_meta))
  if (length(missing)) stop("Missing ADT meta columns: ", paste(missing, collapse = ", "))
  
  # long table (group + feature + value)
  df <- df_meta %>%
    dplyr::select(all_of(c(grp_name, features))) %>%
    tidyr::pivot_longer(cols = all_of(features),
                        names_to = "feature", values_to = "value")
  
  # ---- robust caps handling (no rowwise / no if inside mutate) ----
  if (!is.null(caps)) {
    # build a tiny lookup table of caps per feature
    caps_df <- tibble::tibble(
      feature = names(caps),
      cap_min = map_dbl(caps, ~ .x[1]),
      cap_max = map_dbl(caps, ~ .x[2])
    )
    # left-join caps and clamp when caps exist
    df <- df %>%
      dplyr::left_join(caps_df, by = "feature") %>%
      dplyr::mutate(
        value = dplyr::if_else(
          !is.na(cap_min) & !is.na(cap_max),
          pmin(pmax(value, cap_min), cap_max),
          value
        )
      ) %>%
      dplyr::select(-cap_min, -cap_max)
  }
  
  # average by group x feature
  avg_mat <- df %>%
    dplyr::group_by(.data[[grp_name]], feature) %>%
    dplyr::summarise(avg = mean(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = .data[[grp_name]], values_from = avg) %>%
    as.data.frame()
  
  rownames(avg_mat) <- avg_mat$feature
  avg_mat$feature <- NULL
  mat <- as.matrix(avg_mat[features, , drop = FALSE])
  
  # ----- scale for color -----
  if (value_scale == "feature_z") {
    row_mu <- rowMeans(mat, na.rm = TRUE)
    row_sd <- apply(mat, 1, sd, na.rm = TRUE)
    z <- sweep(mat, 1, row_mu, "-")
    z <- sweep(z, 1, ifelse(is.na(row_sd) | row_sd == 0, 1, row_sd), "/")
    z[is.na(z)] <- 0
    z <- pmin(pmax(z, -3), 3)
    mat_col <- (z + 3) / 6   # map [-3,3] -> [0,1]
    pal <- colorRampPalette(palette)(200)
    breaks <- seq(0, 1, length.out = 201)
  } else {
    rng <- range(mat, na.rm = TRUE)
    if (!is.finite(diff(rng)) || diff(rng) == 0) {
      mat_col <- matrix(0.5, nrow = nrow(mat), ncol = ncol(mat),
                        dimnames = dimnames(mat))
    } else {
      mat_col <- (mat - rng[1]) / diff(rng)
    }
    pal <- colorRampPalette(palette)(200)
    breaks <- seq(0, 1, length.out = 201)
  }
  
  # preserve group order if grouping column is a factor
  if (is.factor(df_meta[[grp_name]])) {
    keep_cols <- intersect(levels(df_meta[[grp_name]]), colnames(mat_col))
    mat_col <- mat_col[, keep_cols, drop = FALSE]
  }
  
  pheatmap(
    mat_col,
    color = pal,
    breaks = breaks,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    border_color = NA,
    fontsize = fontsize,
    angle_col = "45",
    labels_row = rownames(mat),       
    labels_col = colnames(mat_col),
    legend = TRUE,
    legend_breaks = c(0, 0.5, 1),
    legend_labels = if (value_scale == "feature_z")
      c("≤ -3 z", "0 z", "≥ +3 z") else c("Low", "Mid", "High")
  )
}

adt_feats <- c(
  "ADT_CD10","ADT_IgM", "ADT_IgD", "ADT_Ig_light_chain_lambda",
  "ADT_Ig_light_chain_kappa", "ADT_CD44")

caps_list <- list(
  ADT_CD10 = c(0.5, 10),
  ADT_CD44 = c(3, 10),
  ADT_IgM  = c(0, 9),
  ADT_IgD  = c(0, 5),
  ADT_Ig_light_chain_lambda = c(0, 9),
  ADT_Ig_light_chain_kappa  = c(0, 4)
)

make_adt_heatmap_meta(
  mbg51_sub, adt_feats,
  group.by = "ident",           # use active identities
  caps = caps_list,
  value_scale = "feature_z",
  palette = c("#2166AC","white","#B2182B"),
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

svg("~/Desktop/MBG51_ADT_heatmap.svg", width = 4.3, height = 2.5, bg = "transparent")
ComplexHeatmap::draw(plot)
dev.off()

ggsave(
  filename = "~/Desktop/MBG51_ADT_heatmap.svg",
  plot     = plot,
  device   = "svg",
  bg       = "transparent",
  width    = 4.3,
  height   = 2.5,
  units    = "in"
)
