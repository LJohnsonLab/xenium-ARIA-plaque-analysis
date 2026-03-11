########################################################################
# Name: CCL_ST_v10_ARIA_02.R
# Project: low dose ARIA (Akhil's project)
# Purpose: Analysis of low dose ARIA Xenium Spatial Transcriptomics data
#         - Part 1: Merging files, attaching metadata and pre-processing 
#         - Part 2: RCTD and SPLIT
#         - Part 3: manually annotating all clusters 
#         - Part 4: subclustering, cleaning and manually annotating microglia 
#         - Part 5: code for publishable figures 
# Input Files: rds files from part 1 and single cell reference dataset (Allen Brain Cell Atlas 10xv3)
# Output Files: RCTD and SPLIT object
# UPDATE FOR SOURCE CODE:
#             - change all "SaveSeuratRds" to "qsave" (saves files as qs rather than rds)
# Date created: 2/20/26
# Last updated: 2/20/27
# Author: Chloe Lucido
########################################################################

# Load Libraries ----
BiocManager::install("spacexr")
library(spacexr)
library(Seurat)
library(SeuratDisk)
library(future)
library(ggplot2)
library(arrow)
library(hdf5r)
library(presto)
library(glmGamPoi)
library(readr)
library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
library(patchwork)
library(sceasy)
library(reticulate)
library(SPLIT)
library(qs2) 
library(RColorBrewer)
library(Polychrome)
library(purrr)
options(future.globals.maxSize = 100 *1024^3)

### 06. RCTD ----
# vignette: https://github.com/bdsc-tds/SPLIT/blob/main/vignettes/Run_RCTD_and_SPLIT_on_Xenium.Rmd

# load ABC 10xv3 reference dataset
ABC_ref <- readRDS("/Users/cclu223/Desktop/ABC_reference/WMB_10xv3_data/downsampled_objs/20260123_WMB_10xv3_FINAL.rds")

# load xenium.aria obj
xenium.aria <- readRDS("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260126_mergedobj_QC.rds")

## create "spatial" reduction for xenium.aria ----
img_names <- names(xenium.aria@images)

# extracting coordinates from images 
coord_list <- lapply(img_names, function(img) {
  coords <- GetTissueCoordinates(xenium.aria, image = img, which = "centroids")
  coords$image <- img
  coords
})

# binding all the coordinates together 
sp_coords <- do.call(rbind, coord_list)

rownames(sp_coords) <- sp_coords$cell

sp_coords <- sp_coords[, c("x", "y")]

# common cells between coords and Cells from xenium.aria
common_cells <- intersect(Cells(xenium.aria), rownames(sp_coords))

# SANITY CHECK
length(common_cells)

# subsetting out only the common cells 
xenium.aria <- subset(xenium.aria, cells = common_cells)
sp_coords   <- sp_coords[common_cells, ]

# SANITY CHECK: ensuring the rownames of the coords and the cell in xenium.aria match 
stopifnot(all(rownames(sp_coords) == Cells(xenium.aria)))

# attaching "ST_" prefix 
colnames(sp_coords) <- paste0("ST_", seq_len(ncol(sp_coords)))


# creating spatial reduction in seurat obj
xenium.aria[["spatial"]] <- CreateDimReducObject(
  embeddings = as.matrix(sp_coords),
  key = "SPATIAL_",
  assay = DefaultAssay(xenium.aria)
)

# attached the coordinates to xenium.aria metadata
xenium.aria$x <- sp_coords[,1]
xenium.aria$y <- sp_coords[,2]


## update ABC_ref ----
ABC_ref <- UpdateSeuratObject(ABC_ref)

Idents(ABC_ref) <- "Subcluster"

# merging L5_ET into L5_IT cluster because L5_ET has <25 cells
ABC_ref$Subcluster[ABC_ref$Subcluster == "L5_ET"] <- "L5_IT"

# broad cluster dataframe 
subcluster_to_broadcluster <- c(
  # "subcluster" = "broad_cluster"
  "NK Cells" = "NK Cells",
  "T Cells" = "T Cells",
  "Microglia" = "Microglia",
  "BAM" = "BAM",
  "Monocytes" = "Monocytes",
  "DC" = "DC",
  "ABC" = "ABC",
  "VLMC" = "VLMC",
  "Pericytes" = "Pericytes",
  "SMC" = "SMC",
  "Endothelial Cells" = "Endothelial Cells",
  "OPC" = "OPC",
  "Oligodendrocytes" = "Oligodendrocytes",
  "Astrocytes" = "Astrocytes",
  "Astroependymal" = "Astroependymal",
  "Tanycytes" = "Tanycytes",
  "Ependymal" = "Ependymal",
  "Hypendymal" = "Hypendymal",
  "CP" = "CP",
  "Dopaminergic_Neurons" = "Dopaminergic_Neurons",
  "Serotonergic_Neurons" = "Serotonergic_Neurons",
  "LAMP5" = "GABAergic Neuron",
  "PVALB" = "GABAergic Neuron",
  "SNCG" = "GABAergic Neuron",
  "SST" = "GABAergic Neuron",
  "VIP" = "GABAergic Neuron",
  "L2_3_IT" = "Glutamatergic Neuron",
  "L5_IT" = "Glutamatergic Neuron",
  "L6_IT" = "Glutamatergic Neuron",
  #"L5_ET" = "Glutamatergic Neuron",
  "L6_CT" = "Glutamatergic Neuron"
)

## defining colors for subclusters DO I NEED THIS? ----
cell_types <- unique(ABC_ref$Subcluster) %>% sort()

colors <- brewer.pal(n = max(3, min(length(cell_types), 12)), name = "Set3")

# Recycle colors if not enough
colors <- rep(colors, length.out = length(cell_types))
pal <- setNames(colors, cell_types)

## defining reference labels as subclusters ----
ref_labels <- as.factor(ABC_ref$Subcluster)

# drop unused levels
ref_labels <- droplevels(ref_labels)

# FILTER BEFORE CONVERTING TO DATAFRAME!
subcluster_to_broadcluster <- subcluster_to_broadcluster[names(subcluster_to_broadcluster) %in% levels(ref_labels)]

# SANITY CHECK: make sure there are no differences between ref_labels levels and subcluster_to_broadcluster names 
missing_types <- setdiff(levels(ref_labels), names(subcluster_to_broadcluster)) # should get character(0) meaning there is no difference 

# coerce Broad_Cluster_df to dataframe
class_df <- data.frame(class = subcluster_to_broadcluster, row.names = names(subcluster_to_broadcluster))

# common genes between ref dataset and xenium ARIA dataset
common_genes <- intersect(rownames(xenium.aria), rownames(ABC_ref))  

## initializing Ref obj and test obj for RCTD ----
# Reference object
ref_counts <- GetAssayData(ABC_ref, assay = "RNA", layer = "counts")[common_genes, , drop = FALSE]
ref_counts <- round(ref_counts) %>% as(., "dgCMatrix")
ref.obj <- Reference(ref_counts, cell_types = ref_labels, min_UMI = 10, require_int = TRUE)

# Test object
test_counts <- GetAssayData(xenium.aria, assay = "Xenium", layer = "counts")[common_genes, , drop = FALSE]
test_counts <- round(test_counts) %>% as(., "dgCMatrix")

# Align cells
common_cells <- intersect(rownames(sp_coords), colnames(test_counts))
coords <- sp_coords[common_cells, ]
test_counts <- test_counts[, common_cells]

# Create SpatialRNA 
test.obj <- SpatialRNA(coords = coords, counts = test_counts, require_int = FALSE)

# SANITY CHECK:

dim(xenium.aria)

nrow(test.obj@coords) # SHOULD HAVE SAME NUMBER OF COLUMNS OUTPUT ABOVE

# Create RCTD
RCTD <- create.RCTD(test.obj, ref.obj, max_cores = 8, class_df = class_df,
                    UMI_min = 10, counts_MIN = 10, CELL_MIN_INSTANCE = 3)

# run RCTD
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

## save RCTD obj ----
SaveSeuratRds(RCTD, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260204_RCTDobj_bothslides.rds")



### 07. SPLIT----
# github: https://github.com/bdsc-tds/SPLIT.git 

#RCTD <- LoadSeuratRds("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Files/20260127_RCTDobj_correctprotocol.rds")

# post process RCTD output
RCTD <- SPLIT::run_post_process_RCTD(RCTD)

# add RCTD results to xenium.aria metadata
xenium.aria <- AddMetaData(xenium.aria, RCTD@results$results_df)

# visualize RCTD annnotation
xenium.aria <- xenium.aria %>% SCTransform(assay = "Xenium", verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30, verbose = FALSE)

## 07a. Run SPLIT purification ON DOUBLETS ONLY ----
res_split <- SPLIT::purify(
  counts = GetAssayData(xenium.aria, assay = 'Xenium', layer = 'counts'), # or any gene x cells counts matrix
  rctd = RCTD,
  DO_purify_singlets = FALSE # Optional. If TRUE, singlets with an available secondary type are purified the same way as doublets_certain; otherwise, left unchanged.
)

# Create a purified Seurat object
xenium_pur <- CreateSeuratObject(
  counts = res_split$purified_counts,
  meta.data = res_split$cell_meta,
  assay = "Xenium"
)



## 07b. Run spatially aware SPLIT ----

# computing spatial neighborhood for each cell
sp_nw <- SPLIT::build_spatial_network(
  xenium.aria, 
  reduction = "spatial",
  dims = 1:2, 
  DO_prune = TRUE, 
  rad_pruning = 15, # remove connections further than 15um
  k_knn = 20
)

sp_nw <- SPLIT::add_spatial_metric(spatial_neighborhood = sp_nw, rctd = RCTD)

sp_neigh_df <- SPLIT::neighborhood_analysis_to_metadata(sp_nw)

# adding to xenium.aria metadata
xenium.aria <- AddMetaData(xenium.aria, sp_neigh_df)

## 07c. purifying cells with secondary signal in their spatial neighborhood (e.g., `neighborhood_weights_second_type`) and keep other cells unchanged ----
xenium_purified_balanced_score <- SPLIT::balance_raw_and_purified_data_by_score(
  xe_raw = xenium.aria,
  xe_purified = xenium_pur,
  default_assay = "Xenium", # should be param, but can wait 
  spot_class_key = "spot_class",
  threshold = 0.05, # lower -> more cells will be purified
  score_name = "neighborhood_weights_second_type"
)



## 07d. SPLIT-SHIFT ----

# compute transcriptomics neighborhood
tr_nw <- build_transcriptomics_network(
  xenium.aria,
  DO_prune = FALSE,
  k_knn = 100
)
tr_nw <- add_transcriptomics_metric(transcriptomics_neighborhood = tr_nw, rctd = RCTD) 
tr_neigh_df <- neighborhood_analysis_to_metadata(tr_nw)
xenium.aria <- AddMetaData(xenium.aria, tr_neigh_df)

rm(tr_nw, tr_neigh_df)

# re-run SPLIT::balance_raw_and_purified_data_by_score but with DO_swap_labels = TRUE

xe_split_shift <- SPLIT::balance_raw_and_purified_data_by_score(
  xe_raw = xenium.aria,
  xe_purified = xenium_pur,
  default_assay = "Xenium",
  spot_class_key = "spot_class",
  threshold = 0.05, # to be consistent with spatially-aware SPLIT results
  score_name = "neighborhood_weights_second_type",
  DO_swap_lables = TRUE
)

# optional post-processing: filter, normalize and visualize
xe_split_shift <- subset(xe_split_shift, subset = nCount_Xenium > 5)
xe_split_shift <- xe_split_shift %>%
  SCTransform(assay = "Xenium", verbose = FALSE) %>%
  RunPCA(npcs = 30, features = rownames(xenium.aria), verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE) # DOUBLE CHECK RESOLUTION, usually do 30-50


### 08. attaching res_split metadata to xenium.aria obj ----
## !!!FOR NEXT TIME, ATTACH METADATA FROM ORIGINAL OBJ TO SPLIT OBJ (easier and less risky)!!!

# only missing purified counts from xenium.aria metadata

# creating column for purification status metadata to go into 
col <- "purification_status"

# intializing the metadata column
xenium.aria$purification_status <- NA

# align barcodes for smooth transfer
cells_split <- colnames(xenium_pur)

xenium.aria@meta.data[
  cells_split,
  "purification_status"
] <- xenium_pur@meta.data[ 
  cells_split,
  "purification_status"
]


# transferring metadata from split shift to og obj 
xenium.aria <- AddMetaData(
  xenium.aria,
  xe_split_shift@meta.data
) 

# transferring split shift assays to og obj

xenium.aria[["Xenium_SPLIT"]] <- xe_split_shift[["Xenium"]]

xenium.aria[["SCT_SPLIT"]] <- xe_split_shift[["SCT"]]

# transferring split shift reductions to og obj

xenium.aria[["pca_SPLIT"]]  <- xe_split_shift[["pca"]]

xenium.aria[["umap_SPLIT"]] <- xe_split_shift[["umap"]]

### 09. save xenium.aria obj with all RCTD and SPLIT shift data attached ----

SaveSeuratRds(xenium.aria, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260209_aria_RCTDSPLIT.rds") 

