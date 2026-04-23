## vdj_usage_plots.R
# Make plots to show variety of VDJ genes being used among hPSC derived B-cells

library(readr)
library(dplyr)
library(tidyr)
library(rlang)
library(ggplot2)
library(scales)

## -------------------------------------------------------------------------
## Load heavy-chain annotations
## -------------------------------------------------------------------------

bm_ann <- read_csv("data/vdj/bm_total/imgt_processed_output.csv", show_col_types = FALSE)
bm_ann <- bm_ann %>%
  filter(chain_type == "igh")

il7_low_ann <- read_csv("data/vdj/il7_low/imgt_processed_output.csv", show_col_types = FALSE)
il7_low_ann <- il7_low_ann %>%
  filter(chain_type == "igh")

## -------------------------------------------------------------------------
## Compare V-gene usage between two groups
## -------------------------------------------------------------------------

plot_v_gene_usage_two <- function(
    df1, df2,
    name1 = "BM_total",
    name2 = "IL7_low_hPSC",
    v_col = "v_gene",
    order = c("none", "alpha", "total", "genomic"),
    flip = TRUE,
    palette = NULL,
    break_at = 0.10,
    compress = 0.25
) {
  order <- match.arg(order)
  
  # 5' to 3' order from IMGT
  # https://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=human&group=IGHVr
  imgt_pos <- c(
    "ighv7-81" = 2,
    "ighv4-39" = 70, "ighv4-38-2" = 79, "ighv4-34" = 86, "ighv4-31" = 93,
    "ighv4-30-4" = 99, "ighv4-30-2" = 105, "ighv4-28" = 113,
    "ighv4-61" = 37, "ighv4-59" = 40, "ighv4-4" = 156,
    "ighv3-74" = 12, "ighv3-73" = 13, "ighv3-72" = 14, "ighv3-64d" = 31,
    "ighv3-66" = 29, "ighv3-64" = 32, "ighv3-53" = 47, "ighv3-49" = 54, "ighv3-48" = 55,
    "ighv3-43" = 65, "ighv3-23" = 64, "ighv3-43d" = 76,
    "ighv3-33" = 89, "ighv3-30-5" = 96, "ighv3-30" = 109,
    "ighv3-21" = 128, "ighv3-20" = 130, "ighv3-15" = 138, "ighv3-13" = 141,
    "ighv3-11" = 144, "ighv3-9" = 146, "ighv3-7" = 150,
    "ighv2-70" = 16, "ighv2-26" = 117, "ighv2-5" = 154,
    "ighv1-69" = 21, "ighv1-69-2" = 18, "ighv1-58" = 41, "ighv1-46" = 59, "ighv1-45" = 60,
    "ighv1-24" = 120, "ighv1-18" = 132, "ighv1-8" = 147, "ighv1-3" = 157, "ighv1-2" = 159,
    "ighv5-51" = 51, "ighv5-10-1" = 50,
    "ighv7-4-1" = 155,
    "ighv6-1" = 163
  )
  
  d1 <- df1 %>%
    filter(!is.na(.data[[v_col]])) %>%
    count(!!sym(v_col), name = "n") %>%
    mutate(source = name1)
  
  d2 <- df2 %>%
    filter(!is.na(.data[[v_col]])) %>%
    count(!!sym(v_col), name = "n") %>%
    mutate(source = name2)
  
  all_v <- union(d1[[v_col]], d2[[v_col]])
  
  d1 <- d1 %>%
    complete(!!sym(v_col) := all_v, fill = list(n = 0)) %>%
    mutate(source = name1)
  
  d2 <- d2 %>%
    complete(!!sym(v_col) := all_v, fill = list(n = 0)) %>%
    mutate(source = name2)
  
  comb <- bind_rows(d1, d2) %>%
    group_by(source) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  ord_levels <- switch(
    order,
    "none" = all_v,
    "alpha" = sort(unique(comb[[v_col]])),
    "total" = comb %>%
      group_by(!!sym(v_col)) %>%
      summarise(total_prop = sum(prop), .groups = "drop") %>%
      arrange(desc(total_prop)) %>%
      pull(!!sym(v_col)),
    "genomic" = {
      v_lower <- tolower(as.character(unique(comb[[v_col]])))
      pos <- imgt_pos[v_lower]
      ordered <- data.frame(
        v = unique(comb[[v_col]]),
        v_lower = v_lower,
        pos = pos,
        stringsAsFactors = FALSE
      ) %>%
        arrange(if_else(is.na(pos), Inf, pos))
      ordered$v
    }
  )
  
  comb[[v_col]] <- factor(comb[[v_col]], levels = ord_levels)
  levels(comb[[v_col]]) <- gsub("^IGHV", "", levels(comb[[v_col]]), ignore.case = TRUE)
  
  palette <- if (is.null(palette)) {
    setNames(c("#D537B7", "turquoise4"), c(name1, name2))
  } else {
    setNames(palette, c(name1, name2))
  }
  
  broken_percent_trans <- scales::trans_new(
    name = "broken_percent",
    transform = function(y) {
      y <- as.numeric(y)
      out <- y
      idx <- !is.na(y) & (y > break_at)
      out[idx] <- break_at + (y[idx] - break_at) * compress
      out
    },
    inverse = function(t) {
      t <- as.numeric(t)
      out <- t
      idx <- !is.na(t) & (t > break_at)
      out[idx] <- break_at + (t[idx] - break_at) / compress
      out
    }
  )
  
  low_breaks <- seq(0, break_at, by = 0.02)
  high_breaks <- c(0.15, 0.20, 0.25, 0.30)
  y_breaks <- unique(c(low_breaks, high_breaks))
  
  p <- ggplot(comb, aes(x = .data[[v_col]], y = prop, fill = source)) +
    geom_col(
      position = position_dodge(width = 0.85),
      width = 0.8,
      colour = NA
    ) +
    scale_y_continuous(
      trans = broken_percent_trans,
      breaks = y_breaks,
      labels = percent_format(accuracy = 1)
    ) +
    scale_fill_manual(values = palette) +
    labs(
      x = "V gene",
      y = "Proportion of sequences",
      fill = "Source"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        size = 9,
        angle = if (flip) 0 else 90,
        vjust = 0.5,
        hjust = if (flip) 0.5 else 1
      )
    )
  
  if (flip) {
    p <- p + coord_flip()
  }
  
  p
}

plot <- plot_v_gene_usage_two(
  bm_ann,
  il7_low_ann,
  name1 = "BM_total",
  name2 = "IL7_low_hPSC",
  v_col = "v_gene",
  order = "genomic",
  flip = FALSE
)

ggsave(
  filename = "figures/v_gene_usage_bm_vs_il7_low.svg",
  plot = plot,
  device = "svg",
  bg = "transparent",
  width = 7.5,
  height = 3,
  units = "in"
)

## -------------------------------------------------------------------------
## Compare D- or J-gene usage between two groups
## -------------------------------------------------------------------------

plot_d_gene_usage_two <- function(
    df1, df2,
    name1 = "BM_total",
    name2 = "IL7_low_hPSC",
    d_col = "d_gene",
    flip = TRUE,
    palette = NULL,
    break_at = 0.10,
    compress = 0.25,
    high_breaks = c(0.15, 0.20, 0.25, 0.30, 0.35)
) {
  d1 <- df1 %>%
    filter(!is.na(.data[[d_col]])) %>%
    count(!!sym(d_col), name = "n") %>%
    mutate(source = name1)
  
  d2 <- df2 %>%
    filter(!is.na(.data[[d_col]])) %>%
    count(!!sym(d_col), name = "n") %>%
    mutate(source = name2)
  
  all_d <- union(d1[[d_col]], d2[[d_col]])
  
  d1 <- d1 %>%
    complete(!!sym(d_col) := all_d, fill = list(n = 0)) %>%
    mutate(source = name1)
  
  d2 <- d2 %>%
    complete(!!sym(d_col) := all_d, fill = list(n = 0)) %>%
    mutate(source = name2)
  
  comb <- bind_rows(d1, d2) %>%
    group_by(source) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  comb[[d_col]] <- factor(comb[[d_col]], levels = all_d)
  levels(comb[[d_col]]) <- gsub("^IGHD", "", levels(comb[[d_col]]), ignore.case = TRUE)
  
  palette <- if (is.null(palette)) {
    setNames(c("#D537B7", "turquoise4"), c(name1, name2))
  } else {
    setNames(palette, c(name1, name2))
  }
  
  broken_percent_trans <- scales::trans_new(
    name = "broken_percent",
    transform = function(y) {
      y <- as.numeric(y)
      out <- y
      idx <- !is.na(y) & (y > break_at)
      out[idx] <- break_at + (y[idx] - break_at) * compress
      out
    },
    inverse = function(t) {
      t <- as.numeric(t)
      out <- t
      idx <- !is.na(t) & (t > break_at)
      out[idx] <- break_at + (t[idx] - break_at) / compress
      out
    }
  )
  
  low_breaks <- seq(0, break_at, by = 0.02)
  y_breaks <- unique(c(low_breaks, high_breaks))
  
  p <- ggplot(comb, aes(x = .data[[d_col]], y = prop, fill = source)) +
    geom_col(
      position = position_dodge(width = 0.85),
      width = 0.8,
      colour = NA
    ) +
    scale_y_continuous(
      trans = broken_percent_trans,
      breaks = y_breaks,
      labels = percent_format(accuracy = 1),
      limits = c(0, max(high_breaks, na.rm = TRUE))
    ) +
    scale_fill_manual(values = palette) +
    labs(
      x = "Gene",
      y = "Proportion of sequences",
      fill = "Source"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        size = 9,
        angle = if (flip) 0 else 90,
        vjust = 0.5,
        hjust = if (flip) 0.5 else 1
      )
    )
  
  if (flip) {
    p <- p + coord_flip()
  }
  
  p
}

plot <- plot_d_gene_usage_two(
  bm_ann,
  il7_low_ann,
  name1 = "BM_total",
  name2 = "IL7_low_hPSC",
  d_col = "primary_j_allele",
  high_breaks = c(0.10, 0.20, 0.30, 0.40, 0.50)
)

ggsave(
  filename = "figures/j_gene_usage_bm_vs_il7_low.svg",
  plot = plot,
  device = "svg",
  bg = "transparent",
  width = 5,
  height = 2,
  units = "in"
)

## -------------------------------------------------------------------------
## Pie chart for top gene usage
## -------------------------------------------------------------------------

plot_pie_proportion <- function(
    df,
    column_name,
    top_n = 12,
    collapse_other = TRUE,
    show_labels_top = 8,
    other_label = "Other",
    other_color = "grey90",
    palette = c("turquoise4", "turquoise"),
    title = NULL
) {
  if (!column_name %in% names(df)) {
    stop("Column ", column_name, " not found in dataframe")
  }
  
  col_sym <- sym(column_name)
  
  prop_df <- df %>%
    count(!!col_sym, name = "n") %>%
    mutate(prop = n / sum(n)) %>%
    arrange(desc(prop))
  
  if (collapse_other && nrow(prop_df) > top_n) {
    head_df <- prop_df %>% slice_head(n = top_n)
    tail_df <- prop_df %>%
      slice((top_n + 1):n()) %>%
      summarise(
        !!column_name := other_label,
        n = sum(n),
        prop = sum(prop)
      )
    plot_df <- bind_rows(head_df, tail_df)
  } else {
    plot_df <- prop_df
  }
  
  gene_levels <- as.character(plot_df[[column_name]])
  if (other_label %in% gene_levels) {
    gene_levels <- c(setdiff(gene_levels, other_label), other_label)
  }
  
  plot_df[[column_name]] <- factor(plot_df[[column_name]], levels = gene_levels)
  
  n_explicit <- nrow(plot_df) - as.integer(any(plot_df[[column_name]] == other_label))
  grad_cols <- colorRampPalette(palette)(max(n_explicit, 1))
  
  fill_vals <- grad_cols
  if (any(plot_df[[column_name]] == other_label)) {
    fill_vals <- c(grad_cols, other_color)
  }
  names(fill_vals) <- gene_levels
  
  ggplot(plot_df, aes(x = "", y = prop, fill = !!col_sym)) +
    geom_bar(stat = "identity", width = 1, color = "grey80") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    labs(
      title = title %||% paste("Proportion of", column_name),
      fill = column_name,
      x = NULL,
      y = NULL
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

plot <- plot_pie_proportion(
  il7_low_ann,
  "j_gene",
  top_n = 5,
  title = "Top 5 J genes in IL7_low_hPSC"
)

ggsave(
  filename = "figures/top5_j_gene_piechart_il7_low.svg",
  plot = plot,
  device = "svg",
  bg = "transparent",
  width = 7,
  height = 4,
  units = "in"
)