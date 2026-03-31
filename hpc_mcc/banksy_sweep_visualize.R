#!/usr/bin/env Rscript
#
# banksy_sweep_visualize.R
#
# Reads banksy_sweep.qs2 produced by banksy_parameter_sweep.R and writes one
# PDF per lambda × resolution combination (16 total).  Each PDF is a single
# 4-panel page:
#   top-left   – UMAP coloured by domain
#   top-right  – spatial plot of all cells coloured by domain (Polychrome palette)
#   bottom-left  – cell-type composition per domain (stacked bar)
#   bottom-right – treatment composition per domain (stacked bar)
#
# Usage:
#   Rscript hpc_mcc/banksy_sweep_visualize.R \
#     --out_dir spatial/banksy_sweep


user_lib <- path.expand("~/R/libs")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(Polychrome)
})

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  args[[idx + 1]]
}

out_dir  <- parse_arg(args, "--out_dir", default = "spatial/banksy_sweep")
plot_dir <- file.path(out_dir, "plots")

message("out_dir : ", out_dir)
message("plot_dir: ", plot_dir)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load object
# ---------------------------------------------------------------------------
message("[", Sys.time(), "] Loading banksy_sweep.qs2...")
obj  <- qs_read(file.path(out_dir, "banksy_sweep.qs2"))
meta <- as_tibble(obj@meta.data, rownames = "cell")

message(sprintf("  %d cells loaded", ncol(obj)))

# ---------------------------------------------------------------------------
# Parameter grid  (must match banksy_parameter_sweep.R exactly)
# ---------------------------------------------------------------------------
lambdas     <- c(0.2, 0.4, 0.6, 0.8)
resolutions <- c(0.1, 0.3, 0.5, 1.0)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
lam_tag <- function(lam) sprintf("%02d", as.integer(round(lam * 10)))
res_tag <- function(res) sprintf("%02d", as.integer(round(res * 10)))

# ---------------------------------------------------------------------------
# Plotting loop
# ---------------------------------------------------------------------------
message("[", Sys.time(), "] Generating PDFs...")

for (lam in lambdas) {
  umap_key <- paste0("umap_banksy_lam", lam_tag(lam))

  for (res in resolutions) {
    meta_col <- sprintf("lambda%.1f_resolution%.1f", lam, res)
    message(sprintf("  %s", meta_col))

    # Attach domain labels to object and metadata tibble
    obj$banksy_domain  <- factor(obj[[meta_col, drop = TRUE]])
    Idents(obj)        <- "banksy_domain"
    meta$banksy_domain <- factor(meta[[meta_col]])

    # Colour palette: one colour per domain (Polychrome kelly palette, max 22;
    # fall back to createPalette for higher-resolution runs with >22 clusters)
    domain_levels <- sort(levels(obj$banksy_domain))
    n_domains     <- length(domain_levels)
    raw_cols <- if (n_domains <= 22) {
      kelly.colors(n_domains)
    } else {
      createPalette(n_domains, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
    }
    domain_cols <- setNames(raw_cols, domain_levels)

    subtitle_text <- sprintf("lambda = %.1f  /  resolution = %.1f", lam, res)

    # -- Plot 1: UMAP --------------------------------------------------------
    p_umap <- DimPlot(
      obj,
      reduction  = umap_key,
      group.by   = "banksy_domain",
      label      = TRUE,
      label.size = 4,
      pt.size    = 0.3,
      cols       = domain_cols
    ) +
      labs(
        title    = "BANKSY tissue domains",
        subtitle = subtitle_text,
        colour   = "Domain"
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")

    # -- Plot 2: Spatial per sample (as in 20_banksy_tissue_domains.qmd) -------
    sample_ids <- sort(unique(meta$sample_ID))
    spatial_plots <- map2(sample_ids, seq_along(sample_ids), \(sid, i) {
      df <- meta |> filter(sample_ID == sid)
      ggplot(df, aes(x = x, y = y, colour = banksy_domain)) +
        geom_point(size = 0.2, alpha = 0.7, show.legend = (i == 1)) +
        scale_colour_manual(
          values = domain_cols,
          limits = domain_levels,
          drop   = FALSE
        ) +
        coord_equal() +
        scale_y_reverse() +
        labs(title = sid, colour = "Domain") +
        theme_void(base_size = 10) +
        theme(legend.position = if (i == 1) "left" else "none")
    })
    p_spatial <- wrap_plots(spatial_plots, nrow = 2) +
      plot_annotation(title = "BANKSY tissue domains per brain") &
      guides(colour = guide_legend(override.aes = list(size = 4)))

    # -- Plot 3: Cell-type composition per domain ----------------------------
    ct_frac <- meta |>
      count(banksy_domain, annotatedclusters) |>
      group_by(banksy_domain) |>
      mutate(frac = n / sum(n)) |>
      ungroup()

    p_domain_comp <- ggplot(ct_frac,
                            aes(x = banksy_domain, y = frac,
                                fill = annotatedclusters)) +
      geom_col() +
      labs(
        title = "Cell-type composition per domain",
        x     = "Domain",
        y     = "Fraction of cells",
        fill  = "Cell type"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "right",
        axis.text.x     = element_text(angle = 45, hjust = 1)
      )

    # -- Plot 4: Treatment composition per domain ----------------------------
    trt_frac <- meta |>
      count(banksy_domain, Treatment.Group) |>
      group_by(banksy_domain) |>
      mutate(frac = n / sum(n)) |>
      ungroup()

    p_trt_comp <- ggplot(trt_frac,
                         aes(x = banksy_domain, y = frac,
                             fill = Treatment.Group)) +
      geom_col() +
      scale_fill_manual(values = c(Adu = "#E07B39", IgG = "#7B68AE")) +
      labs(
        title = "Treatment composition per domain",
        x     = "Domain",
        y     = "Fraction of cells",
        fill  = "Treatment"
      ) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # -- Assemble panel and save PDF -----------------------------------------
    # Top row: UMAP (small) + two composition bars; bottom row: spatial (tall)
    top_row <- (p_umap | p_domain_comp | p_trt_comp)
    panel   <- top_row / p_spatial +
      plot_layout(heights = c(1, 2)) +
      plot_annotation(
        title = sprintf("BANKSY sweep - lambda = %.1f / resolution = %.1f", lam, res)
      )

    pdf_name <- sprintf("banksy_lam%s_res%s.pdf", lam_tag(lam), res_tag(res))
    pdf_path <- file.path(plot_dir, pdf_name)

    ggsave(pdf_path, panel, width = 18, height = 16, device = "pdf")
    message(sprintf("    -> %s", pdf_path))
  }
}

message("[", Sys.time(), "] Done. 16 PDFs written to ", plot_dir)
