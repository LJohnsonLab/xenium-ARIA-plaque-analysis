# setup.R — shared setup for all analysis notebooks
# Source with: source("R/setup.R")

# ── Core libraries ────────────────────────────────────────────────────────────
library(readxl)
library(tidyverse)
library(kableExtra)
library(Seurat)
library(RANN)
options(scipen = 999)

# ── Paths ─────────────────────────────────────────────────────────────────────
seurat_path <- "/Users/arbones/Dropbox/SyncBriefcase/LAB/UK/Xenium/ARIA_Akhil/aria_analysis/seurat_objects/"

# ── Sample dictionary ─────────────────────────────────────────────────────────
sample_dict <- tribble(
  ~brain,   ~sample,   ~treatment, ~slide,
  "Brain1", "KK4_465", "IgG",      1L,
  "Brain2", "KK4_492", "Adu",      1L,
  "Brain3", "KK4_504", "IgG",      1L,
  "brain4", "KK4_496", "IgG",      2L,
  "brain5", "KK4_502", "Adu",      2L,
  "brain6", "KK4_464", "Adu",      2L
)

# ── CAA labels (all 6 brains in one file) ─────────────────────────────────────
caa_path <- "HALO_thioS_Akhil/20260217_CAAquant_AP.xlsx"

caa_akhil <- excel_sheets(caa_path) |>
  set_names() |>
  map(\(sheet) read_excel(caa_path, sheet = sheet)) |>
  bind_rows(.id = "sheet") |>
  mutate(
    brain = str_remove(sheet, "_AP"),
    index = paste0(tolower(brain), "_", Id)
  ) |>
  pull(index)

# ── HALO plaque CSVs ──────────────────────────────────────────────────────────
thios_all <- list.files(
  "HALO_thioS_Akhil",
  pattern    = "^20260220.*AP\\.csv$",
  full.names = TRUE
) |>
  map(\(x) read_csv(x, locale = locale(encoding = "latin1"), show_col_types = FALSE)) |>
  bind_rows() |>
  rename_with(\(x) str_remove(x, "\\s*\\(.*\\)$")) |>
  dplyr::rename(brain = "Analysis Region", objectID = "Object Id") |>
  mutate(
    index       = paste0(tolower(brain), "_", objectID),
    plaque_type = if_else(index %in% caa_akhil, "CAA", "Parenchymal")
  ) |>
  left_join(sample_dict, by = "brain")

thios_1 <- thios_all |> filter(slide == 1)
thios_2 <- thios_all |> filter(slide == 2)

# ── Plaque centroids (pixel to µm, factor = 0.2125) ──────────────────────────
make_centroids <- function(thios) {
  thios |>
    mutate(
      x_plaque = ((XMin + XMax) / 2) * 0.2125,
      y_plaque = ((YMin + YMax) / 2) * 0.2125,
      plaque_id = index
    ) |>
    dplyr::select(plaque_id, brain, sample, treatment, plaque_type, x_plaque, y_plaque)
}

plaque_centroids_1 <- make_centroids(thios_1)
plaque_centroids_2 <- make_centroids(thios_2)

caa_coords_1   <- plaque_centroids_1 |> filter(plaque_type == "CAA")
caa_coords_2   <- plaque_centroids_2 |> filter(plaque_type == "CAA")
paren_coords_1 <- plaque_centroids_1 |> filter(plaque_type == "Parenchymal")
paren_coords_2 <- plaque_centroids_2 |> filter(plaque_type == "Parenchymal")

# ── Cell-to-plaque distances (k-d tree) ───────────────────────────────────────
compute_cell_distances <- function(seurat_obj, caa_coords, paren_coords) {
  coords_df <- GetTissueCoordinates(seurat_obj) |>
    as.data.frame() |>
    dplyr::select(cell, x, y)

  cell_mat  <- coords_df |> dplyr::select(x, y) |> as.matrix()
  caa_mat   <- caa_coords   |> dplyr::select(x_plaque, y_plaque) |> as.matrix()
  paren_mat <- paren_coords |> dplyr::select(x_plaque, y_plaque) |> as.matrix()

  caa_nn   <- nn2(data = caa_mat,   query = cell_mat, k = 1)
  paren_nn <- nn2(data = paren_mat, query = cell_mat, k = 1)

  coords_df |>
    mutate(
      dist_caa         = caa_nn$nn.dist[, 1],
      nearest_caa_id   = caa_coords$plaque_id[caa_nn$nn.idx[, 1]],
      dist_paren       = paren_nn$nn.dist[, 1],
      nearest_paren_id = paren_coords$plaque_id[paren_nn$nn.idx[, 1]]
    )
}

# ── Load Seurat + compute distances ───────────────────────────────────────────

###### new full object from Chloe after catchin a bug in slide assignation ######
# aria <- qs2::qs_read(paste0(seurat_path, "20260320_fullobj_correctedChloe.qs2"))
# slide1 <- subset(aria, subset = sample_ID%in%unique(thios_1$sample))
# slide2 <- subset(aria, subset = sample_ID%in%unique(thios_2$sample))
# qs2::qs_save(slide1,file = paste0(seurat_path, "20260320_fullobj_correctedChloe_slide1.qs2"))
# qs2::qs_save(slide2,file = paste0(seurat_path, "20260320_fullobj_correctedChloe_slide2.qs2"))

###### new microglia object from Chloe after catchin a bug in slide assignation ######
# mg <- qs2::qs_read(paste0(seurat_path, "20260323_updated_mg_kai_correctedChloe.qs2"))
# mg <- UpdateSeuratObject(mg)
# slide1_mg <- subset(mg, subset = sample_ID%in%unique(thios_1$sample))
# slide2_mg <- subset(mg, subset = sample_ID%in%unique(thios_2$sample))
# qs2::qs_save(slide1_mg,file = paste0(seurat_path, "20260323_updated_mg_kai_correctedChloe.qs2_slide1.qs2"))
# qs2::qs_save(slide2_mg,file = paste0(seurat_path, "20260323_updated_mg_kai_correctedChloe.qs2_slide2.qs2"))

load_slides_and_distances <- function(type = c("full", "microglia")) {
  type <- match.arg(type)

  if (type == "full") {
    s1 <- qs2::qs_read(paste0(seurat_path, "20260320_fullobj_correctedChloe_slide1.qs2"))
    s2 <- qs2::qs_read(paste0(seurat_path, "20260320_fullobj_correctedChloe_slide2.qs2"))
  } else {
    s1 <- qs2::qs_read(paste0(seurat_path, "20260323_updated_mg_kai_correctedChloe.qs2_slide1.qs2"))
    s2 <- qs2::qs_read(paste0(seurat_path, "20260323_updated_mg_kai_correctedChloe.qs2_slide2.qs2"))
  }

  d1 <- compute_cell_distances(s1, caa_coords_1, paren_coords_1)
  d2 <- compute_cell_distances(s2, caa_coords_2, paren_coords_2)

  d_all <- bind_rows(
    d1 |> mutate(slide = 1L),
    d2 |> mutate(slide = 2L)
  ) |>
    mutate(dist_any = pmin(dist_caa, dist_paren))

  treatment_lookup <- bind_rows(
    s1@meta.data |> as_tibble(rownames = "cell") |> dplyr::select(cell, sample_ID, Treatment.Group),
    s2@meta.data |> as_tibble(rownames = "cell") |> dplyr::select(cell, sample_ID, Treatment.Group)
  )

  d_all <- d_all |> left_join(treatment_lookup, by = "cell")

  list(
    slide1 = s1, slide2 = s2,
    distances_1 = d1, distances_2 = d2,
    distances_all = d_all
  )
}

# ── Gene-distance Spearman correlation ────────────────────────────────────────
# Vectorised: ranks all genes in one BLAS cor() call, then applies the
# t-approximation used internally by cor.test(exact = FALSE).
# Default min_expressing_cells = 100 for all-cell-type notebooks (6, 10);
# pass min_expressing_cells = 5 for module-level notebooks (10.2).
compute_gene_distance_correlation <- function(
    seurat_list,
    cell_distances,
    distance_col         = "dist_any",
    min_expressing_cells = 100
) {
  count_list   <- map(seurat_list, \(obj) LayerData(obj, assay = DefaultAssay(obj), layer = "counts"))
  shared_genes <- Reduce(intersect, map(count_list, rownames))
  counts       <- do.call(cbind, map(count_list, \(m) m[shared_genes, ]))

  shared_cells <- intersect(colnames(counts), cell_distances$cell)
  counts   <- counts[, shared_cells]
  dist_vec <- cell_distances |>
    filter(cell %in% shared_cells) |>
    arrange(match(cell, shared_cells)) |>
    pull(!!sym(distance_col))

  keep   <- Matrix::rowSums(counts > 0) >= min_expressing_cells
  counts <- counts[keep, ]

  n          <- ncol(counts)
  expr_ranks <- apply(as.matrix(counts), 1, rank)
  dist_rank  <- rank(dist_vec)
  rho_vec    <- cor(expr_ranks, dist_rank)[, 1]
  t_stat     <- rho_vec * sqrt((n - 2) / (1 - rho_vec^2))
  pval_vec   <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)

  tibble(gene = names(rho_vec), rho = rho_vec, pval = pval_vec) |>
    mutate(
      padj     = p.adjust(pval, method = "BH"),
      rho_flip = -rho,
      rank     = rank(-rho_flip, ties.method = "first"),
      sig      = case_when(
        padj < 0.05 & rho_flip > 0 ~ "enriched_near",
        padj < 0.05 & rho_flip < 0 ~ "enriched_far",
        TRUE                        ~ "ns"
      )
    ) |>
    arrange(rank)
}
