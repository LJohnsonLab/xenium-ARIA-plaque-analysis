#!/usr/bin/env Rscript
#
# banksy_parameter_sweep.R
#
# BANKSY lambda × resolution parameter sweep for HPC.
#
# For each lambda in {0.2, 0.4, 0.6, 0.8}:
#   - Runs RunBanksy → RunPCA → RunUMAP → FindNeighbors
#   For each resolution in {0.1, 0.3, 0.5, 1.0}:
#     - Runs FindClusters → FindAllMarkers
#
# Outputs (written to --out_dir):
#   markers_all_combinations.csv  — all 16 runs stacked; key cols: lambda, resolution, cluster, gene
#   sweep_cluster_assignments.csv — one row per cell; cols: cell, lambda0.2_resolution0.1, ...
#   banksy_sweep.qs2              — aria_sct_split with all 16 combo columns added to metadata
#
# Usage:
#   Rscript R/banksy_parameter_sweep.R \
#     --input  /path/to/object.qs2 \
#     --out_dir spatial/banksy_sweep


user_lib <- path.expand("~/R/libs")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(Banksy)
  library(qs2)
  library(dplyr)
  library(readr)
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

input_file <- parse_arg(args, "--input",
  default = "20260320_fullobj_correctedChloe_sct_SPLIT.qs2")
out_dir    <- parse_arg(args, "--out_dir", default = "spatial/banksy_sweep")

message("input  : ", input_file)
message("out_dir: ", out_dir)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parameter grid
# ---------------------------------------------------------------------------
lambdas     <- c(0.2,0.4,0.6,0.8)
resolutions <- c(0.1,0.3,0.5,1)

# ---------------------------------------------------------------------------
# Load data once
# ---------------------------------------------------------------------------
message("[", Sys.time(), "] Loading Seurat object...")
t0_load <- proc.time()
aria_sct_split <- qs_read(input_file)
DefaultAssay(aria_sct_split) <- "SCT_SPLIT"

message(sprintf("  Loaded %d cells in %.1f s", ncol(aria_sct_split),
                (proc.time() - t0_load)[["elapsed"]]))

# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
all_markers <- vector("list", length(lambdas) * length(resolutions))
combo_idx <- 1L

for (lam in lambdas) {

  message(sprintf("\n[%s] === lambda = %s ===", Sys.time(), lam))
  t0_banksy <- proc.time()

  obj <- RunBanksy(
    aria_sct_split,
    lambda      = lam,
    assay       = "SCT_SPLIT",
    slot        = "data",
    features    = "all",
    k_geom      = 15,
    dimx        = "x",
    dimy        = "y",
    use_agf     = TRUE,
    group       = "Slide",
    split.scale = TRUE
  )

  message(sprintf("  RunBanksy done in %.1f s", (proc.time() - t0_banksy)[["elapsed"]]))

  # Active assay is now BANKSY; do NOT call ScaleData
  obj <- RunPCA(obj,
                assay    = "BANKSY",
                features = rownames(obj),  # all BANKSY features
                npcs     = 50)

  # Lambda-specific names prevent key conflicts when reductions are copied to aria_sct_split
  lam_tag   <- gsub("\\.", "", as.character(lam))  # "02", "04", "06", "08"
  umap_name <- paste0("umap_banksy_lam", lam_tag)  # e.g. "umap_banksy_lam02"
  umap_key  <- paste0("BKSY", lam_tag, "UMAP_")    # e.g. "BKSY02UMAP_"

  obj <- RunUMAP(obj, dims = 1:30, seed.use = 42,
                 reduction.name = umap_name,
                 reduction.key  = umap_key)

  obj <- FindNeighbors(obj, dims = 1:30)

  # Copy UMAP reduction to base object so all lambdas are accessible for visualization
  aria_sct_split[[umap_name]] <- obj[[umap_name]]

  message(sprintf("  RunPCA + RunUMAP + FindNeighbors done in %.1f s",
                  (proc.time() - t0_banksy)[["elapsed"]]))

  for (res in resolutions) {

    message(sprintf("[%s]   resolution = %s", Sys.time(), res))
    t0_res <- proc.time()

    obj <- FindClusters(obj, graph.name = "BANKSY_snn", resolution = res, random.seed = 42)

    # Column name written by FindClusters
    clust_col <- paste0("BANKSY_snn_res.", res)
    Idents(obj) <- clust_col
    n_clust <- length(levels(obj))

    DefaultAssay(obj) <- "SCT_SPLIT"
    markers <- FindAllMarkers(
      obj,
      assay           = "SCT_SPLIT",
      only.pos        = TRUE,
      min.pct         = 0.1,
      logfc.threshold = 0.25,
      test.use        = "wilcox"
    )

    elapsed <- (proc.time() - t0_res)[["elapsed"]]
    message(sprintf("    %d clusters, %d markers in %.1f s",
                    n_clust, nrow(markers), elapsed))

    # Annotate markers with sweep parameters
    markers <- markers |>
      mutate(lambda = lam, resolution = res, .before = 1)
    all_markers[[combo_idx]] <- markers

    # Write cluster assignments back to the base object as a metadata column
    meta_col <- sprintf("lambda%.1f_resolution%.1f", lam, res)
    aria_sct_split[[meta_col]] <- obj[[clust_col]][Cells(aria_sct_split), , drop = FALSE]

    combo_idx <- combo_idx + 1L
  }

  rm(obj)
  gc()
}

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
message(sprintf("\n[%s] Writing outputs...", Sys.time()))

markers_combined <- bind_rows(all_markers)
write_csv(markers_combined, file.path(out_dir, "markers_all_combinations.csv"))

# Export per-cell cluster assignments as CSV
sweep_cols <- grep("^lambda", colnames(aria_sct_split@meta.data), value = TRUE)
cluster_assignments <- aria_sct_split@meta.data[, sweep_cols, drop = FALSE] |>
  tibble::rownames_to_column("cell")
write_csv(cluster_assignments, file.path(out_dir, "sweep_cluster_assignments.csv"))

# Save Seurat object with all 16 metadata columns attached
qs_save(aria_sct_split, file.path(out_dir, "banksy_sweep.qs2"))

message("Done.")
message("  ", file.path(out_dir, "markers_all_combinations.csv"),
        " (", nrow(markers_combined), " rows)")
message("  ", file.path(out_dir, "sweep_cluster_assignments.csv"),
        " (", ncol(cluster_assignments) - 1L, " combo columns, ", nrow(cluster_assignments), " cells)")
message("  ", file.path(out_dir, "banksy_sweep.qs2"))
