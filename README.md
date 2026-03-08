# Xenium ARIA Plaque Analysis

Co-registration of HALO ThioS plaque annotations with 10x Xenium spatial transcriptomics in 6 mouse brains (3 Aducanumab, 3 IgG control). Plaques are classified as CAA or parenchymal based on single-annotator ThioS segmentation (Akhil Palleria).

## Results

**Website:** [https://ljohnsonlab.github.io/xenium-ARIA-plaque-analysis/](https://ljohnsonlab.github.io/xenium-ARIA-plaque-analysis/)

## Notebooks

| Notebook | Description |
|----------|-------------|
| `2_coregistration_Akhil.qmd` | HALO co-registration and plaque classification |
| `6_gene_distance_correlation.qmd` | Gene–distance Spearman correlations (all cell types) |
| `6.1_cell_type_distance_ridges.qmd` | Cell-to-plaque distance ridge plots |
| `7_gene_distance_correlation_microglia.qmd` | Microglia-specific correlations |
| `8_microglia_proximal_treatment_DE.qmd` | Pseudobulk DESeq2: Adu vs IgG in proximal microglia |
| `9_multi_plaque_neighbourhood.qmd` | Multi-plaque neighbourhood and density effects |
| `10_global_spatial_autocorrelation.qmd` | Global Moran's I and spatial gene structure |
| `11_neighborhood_enrichment.qmd` | Cell-type proximity enrichment and niche composition |

## Data

- **Xenium**: 6 mouse brains across 2 slides (`aria_analysis/seurat_objects/`)
- **HALO annotations**: `HALO_thioS_Akhil/` (CAA Excel + parenchymal CSVs)

## Rendering

```bash
quarto render          # renders all notebooks to docs/
```
