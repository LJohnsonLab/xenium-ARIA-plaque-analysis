# Co-registration of HALO ThioS Plaque Annotations with Xenium Spatial Transcriptomics

## Summary

This project co-registers HALO-segmented ThioS-positive amyloid plaques with 10x Xenium single-cell spatial transcriptomics data from 6 mouse brains treated with either Aducanumab (Adu, n = 3) or IgG control (n = 3). Plaques are classified as **cerebral amyloid angiopathy (CAA)** or **parenchymal** based on manual annotation by a single observer (Akhil Palleria). The analysis examines spatial relationships between cell types and plaques, gene expression gradients as a function of plaque proximity, and treatment-dependent transcriptional changes in plaque-associated microglia.

## Key Findings

- **Plaque burden** is comparable between Adu and IgG groups (no significant difference in total plaque area, CAA area, or parenchymal area).
- **Gene-distance correlations** (all cell types) show neuronal markers (Lamp5, Stard8, Cux2) enriched near CAA plaques and immune/microglial markers (Cst7, Trem2, Cd68, Lyz2) enriched near parenchymal plaques.
- **Microglia-specific correlations** confirm immune activation genes (Cst7, Cd68, Bcl2a1b) are upregulated in microglia proximal to both plaque types, consistent with disease-associated microglia (DAM) signatures.
- **Treatment DE** in CAA-proximal microglia identifies **Il1b** as upregulated under Adu treatment (log2FC = 1.16, p = 0.005), suggesting a treatment-dependent inflammatory response at vascular amyloid deposits.
- **Ridge plot analysis** reveals a treatment-dependent spatial shift: Adu-treated cells tend to sit closer to CAA plaques than IgG controls, while parenchymal plaque distances are largely unaffected.

## Repository Structure

```
├── HALO_thioS_Akhil/          # HALO plaque segmentation CSVs and CAA Excel
├── 2_coregistration_Akhil.qmd # Co-registration pipeline and plaque classification
├── 6_gene_distance_correlation.qmd        # Gene-distance Spearman correlations (all cells)
├── 6.1_cell_type_distance_ridges.qmd      # Distance ridge plots per cell type
├── 7_gene_distance_correlation_microglia.qmd  # Microglia gene-distance + treatment DE
├── 8_microglia_proximal_treatment_DE.qmd  # Pseudobulk DESeq2: Adu vs IgG proximal microglia
├── CLAUDE.md                  # Claude Code project instructions
└── README.md                  # This file
```

## Future Lines of Work

### 1. Expand treatment DE beyond microglia

The current pseudobulk DE is restricted to microglia. Other cell types that cluster near plaques (astrocytes, oligodendrocytes, OPCs) may also show treatment-dependent transcriptional changes. Running Adu vs IgG DE for each cell type proximal to CAA and parenchymal plaques would reveal whether the Il1b signal is microglia-specific or part of a broader inflammatory programme.

### 2. Plaque-size-stratified analysis

The current pipeline treats all plaques equally regardless of size. Stratifying by plaque area (e.g., small vs large tertiles) could reveal whether gene expression gradients and treatment effects scale with plaque burden. Larger plaques may recruit more activated microglia or elicit stronger DAM signatures.

### 3. Multi-plaque neighbourhood effects

Each cell is currently assigned to its single nearest plaque. Cells surrounded by multiple plaques (high local plaque density) may behave differently from cells near an isolated plaque. Computing a local plaque density metric (e.g., number of plaques within 100 um) and incorporating it as a covariate could disentangle proximity from density effects.

### 4. Spatial autocorrelation and niche analysis

Moran's I or similar spatial statistics could test whether gene expression is spatially autocorrelated beyond what plaque proximity alone explains. Combining with cell neighbourhood graphs (e.g., which cell types co-localise within 30 um of a plaque) would move from univariate distance metrics toward a niche-level description.

### 5. Pathway and gene-set enrichment

The Xenium panel (~370 genes) limits conventional pathway analysis, but curated gene sets (DAM signature, complement cascade, cytokine signalling) could be tested for coordinated enrichment near plaques. A gene-set-level correlation with distance would be more robust than single-gene tests given the panel size.

### 6. Interaction modelling (distance x treatment)

The current analyses test distance and treatment separately. A formal interaction model (gene ~ distance * treatment) would identify genes whose distance gradient differs between Adu and IgG, capturing treatment effects that are specifically spatial rather than global.

### 7. Validation with larger cohorts

With n = 3 per group, statistical power is limited. Replicating the Il1b finding and other top hits in an independent cohort or with additional biological replicates would strengthen the conclusions. The pseudobulk framework is ready to accommodate additional samples.

### 8. CAA subtype characterisation

CAA plaques were annotated as a single category. If morphological or molecular subtypes exist (e.g., leptomeningeal vs cortical arteriolar), further stratification could reveal subtype-specific microglial responses and treatment sensitivities.
