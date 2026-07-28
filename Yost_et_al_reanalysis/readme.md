### Overview

This repository contains an R workflow for comparing gene expression between complete responders (CR) and non-responders (NR) in pre-treatment CD8 T cells from the Yost melanoma single-cell RNA-sequencing dataset.

The analysis:

- Loads a preprocessed Seurat object containing pre-treatment CD8 T cells.
- Retains cells from patients classified as CR or NR.
- Evaluates expression of selected T-cell regulatory and effector genes.
- Performs two-sided Wilcoxon rank-sum tests.
- Calculates group-specific arithmetic mean expression.
- Generates stacked and gene-specific violin plots.
- Exports cell-level expression values for reproducibility and downstream analysis.

The genes included in the current workflow are:

- `BCL11B`
- `BACH2`
- `SATB1`
- `PRF1`
- `GZMB`

## Repository Structure

```text
.
├── README.md
├── scripts/
│   └── Yost_CD8_Pre_CR_vs_NR_violin.R
├── data/
│   └── CD8_Pre_Yost_Tcell_combined_seurat_20260720.rds
├── figures/
└── results/
```

The input Seurat object may be stored elsewhere, but the file path in the R script must be updated accordingly.

## Input Data

The workflow requires the following RDS file:

```text
CD8_Pre_Yost_Tcell_combined_seurat_20260720.rds
```

The file must contain a Seurat object with:

1. Normalized gene-expression data.
2. Gene symbols stored as feature names.
3. A metadata column named `outcome`.
4. Outcome labels containing `CR` and `NR`.

The expected outcome definitions are:

- `CR`: complete responder
- `NR`: non-responder

Cells with outcome labels other than `CR` or `NR` are excluded.

## Software Requirements

### R

The workflow is intended for R version 4.3 or later.

### Required R Packages

```r
Seurat
harmony
scCustomize
ggplot2
dplyr
patchwork
ggbeeswarm
svglite
```

Install CRAN packages with:

```r
install.packages(c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "patchwork",
  "ggbeeswarm",
  "svglite",
  "remotes"
))
```

Install Harmony with:

```r
install.packages("harmony")
```

Install `scCustomize` with:

```r
remotes::install_github("samuel-marsh/scCustomize")
```

Although `harmony` and `scCustomize` are loaded in the current script, the CR-versus-NR violin-plot section does not directly call functions from these packages. They are retained to match the analysis environment used to generate the input Seurat object.

## Analysis Workflow

### 1. Confirm Gene Availability

The script checks whether each requested gene is present in the Seurat object:

```r
genes_present <- genes[
  genes %in% rownames(CD8_Pre_Yost_Tcell_combined_seurat)
]

genes_missing <- setdiff(genes, genes_present)
print(genes_missing)
```

Missing genes are reported and excluded from subsequent analyses.

### 2. Restrict the Analysis to CR and NR Cells

```r
CD8_Pre_CRNR <- subset(
  CD8_Pre_Yost_Tcell_combined_seurat,
  subset = outcome %in% c("CR", "NR")
)
```

The plotting order is set to NR followed by CR:

```r
CD8_Pre_CRNR$outcome <- factor(
  CD8_Pre_CRNR$outcome,
  levels = c("NR", "CR")
)
```

### 3. Statistical Analysis

For each gene, expression values are extracted with `FetchData()` and compared between CR and NR cells using:

```r
wilcox.test(expression ~ outcome, data = df)
```

The script reports:

- Arithmetic mean expression in CR cells
- Arithmetic mean expression in NR cells
- Unadjusted Wilcoxon rank-sum test p-value

Although the script defines a `geo_mean()` function, the current pre-treatment analysis uses the arithmetic mean:

```r
cr_mean <- mean(df$expression[df$outcome == "CR"])
nr_mean <- mean(df$expression[df$outcome == "NR"])
```

The analysis is conducted at the cell level. Therefore, cells from the same patient are treated as separate observations. For patient-level inference, pseudobulk or patient-level summary analyses should be considered.

### 4. Visualization

For each gene, the script generates a Seurat violin plot with cell-level points added using `geom_quasirandom()`.

Each plot displays:

- Gene name
- Wilcoxon p-value
- Mean expression in NR cells
- Mean expression in CR cells

The individual plots are vertically combined using:

```r
wrap_plots(plot_list, ncol = 1)
```

## Outputs

### Stacked Violin Plot

The workflow produces one stacked SVG figure containing all analyzed genes:

```text
figures/Yost_et_al_CD8_Pre_CR_vs_NR_stacked_violin.svg
```

### Gene-specific Violin Plots

One SVG file is generated for each gene:

```text
figures/Yost_et_al_CD8_Pre_BCL11B_CR_vs_NR_violin.svg
figures/Yost_et_al_CD8_Pre_BACH2_CR_vs_NR_violin.svg
figures/Yost_et_al_CD8_Pre_SATB1_CR_vs_NR_violin.svg
figures/Yost_et_al_CD8_Pre_PRF1_CR_vs_NR_violin.svg
figures/Yost_et_al_CD8_Pre_GZMB_CR_vs_NR_violin.svg
```

The updated plotting section may append `_v2` to the file name:

```text
figures/Yost_et_al_CD8_Pre_<GENE>_CR_vs_NR_violin_v2.svg
```

### Cell-level Expression Tables

One CSV file is generated for each gene:

```text
results/Yost_et_al_CD8_Pre_<GENE>_CR_vs_NR_raw_values.csv
```

Each file contains:

| Column | Description |
|---|---|
| `Cell` | Seurat cell barcode or cell identifier |
| `outcome` | Clinical outcome category, CR or NR |
| `<GENE>` | Normalized expression value for the specified gene |

### Statistical Summary

The in-memory `pvalue_table` object contains:

| Column | Description |
|---|---|
| `Gene` | Gene symbol |
| `CR_GMean` | Arithmetic mean expression in CR cells |
| `NR_GMean` | Arithmetic mean expression in NR cells |
| `Pvalue` | Unadjusted Wilcoxon p-value |
| `Plot_Title` | Formatted plot title |

For clarity, the names `CR_GMean` and `NR_GMean` are retained from the original script, but the values are arithmetic means rather than geometric means.

The table can be exported with:

```r
write.csv(
  pvalue_table,
  file = "results/Yost_et_al_CD8_Pre_CR_vs_NR_statistics.csv",
  row.names = FALSE
)
```
# Reproducibility

To ensure reproducibility, record your software environment:

```r
R version 4.5.1 (2025-06-13)
Platform: x86_64-apple-darwin20
Running under: macOS Tahoe 26.5.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] shiny_1.11.1                DT_0.33                     effectsize_1.0.3            emmeans_2.0.2               car_3.1-3                  
 [6] carData_3.0-5               uwot_0.2.3                  pheatmap_1.0.13             org.Hs.eg.db_3.21.0         AnnotationDbi_1.70.0       
[11] GEOquery_2.77.1             biomaRt_2.64.0              future_1.67.0               harmony_2.0.2               Rcpp_1.1.0                 
[16] Nebulosa_1.18.0             scCustomize_3.2.4           lubridate_1.9.4             forcats_1.0.1               purrr_1.1.0                
[21] tidyr_1.3.1                 tidyverse_2.0.0             reshape_0.8.10              reshape2_1.4.4              reticulate_1.43.0          
[26] SeuratWrappers_0.4.0        stringr_1.5.1               ggpubr_0.6.1                UCell_2.12.0                ggrepel_0.9.6              
[31] cowplot_1.2.0               ggplot2_3.5.2               SingleCellExperiment_1.30.1 SingleR_2.10.0              SeuratData_0.2.2.9002      
[36] tibble_3.3.0                SeuratDisk_0.0.0.9021       fields_17.1                 RColorBrewer_1.1-3          viridisLite_0.4.2          
[41] spam_2.11-1                 Matrix_1.7-3                celldex_1.18.0              SummarizedExperiment_1.38.1 Biobase_2.68.0             
[46] GenomicRanges_1.60.0        GenomeInfoDb_1.44.1         IRanges_2.42.0              S4Vectors_0.46.0            BiocGenerics_0.54.0        
[51] generics_0.1.4              MatrixGenerics_1.20.0       matrixStats_1.5.0           DoubletFinder_2.0.6         readr_2.1.5                
[56] patchwork_1.3.1             Seurat_5.5.1                SeuratObject_5.2.0          sp_2.2-0                    dplyr_1.2.1                

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2         vroom_1.6.5               progress_1.2.3            goftest_1.2-3             Biostrings_2.76.0        
  [6] HDF5Array_1.36.0          TH.data_1.1-5             vctrs_0.7.3               spatstat.random_3.4-5     digest_0.6.37            
 [11] png_0.1-8                 shape_1.4.6.1             gypsum_1.4.0              bayestestR_0.18.1         deldir_2.0-4             
 [16] parallelly_1.45.1         MASS_7.3-65               httpuv_1.6.16             withr_3.0.2               ggrastr_1.0.2            
 [21] survival_3.8-3            memoise_2.0.1             ggbeeswarm_0.7.3          janitor_2.2.1             parameters_0.29.2        
 [26] zoo_1.8-14                GlobalOptions_0.1.2       pbapply_1.7-4             R.oo_1.27.1               Formula_1.2-5            
 [31] prettyunits_1.2.0         datawizard_1.3.1          rematch2_2.1.2            KEGGREST_1.48.1           promises_1.3.3           
 [36] httr_1.4.7                rstatix_0.7.2             globals_0.18.0            fitdistrplus_1.2-4        rhdf5filters_1.20.0      
 [41] rhdf5_2.52.1              rstudioapi_0.17.1         UCSC.utils_1.4.0          miniUI_0.1.2              curl_6.4.0               
 [46] h5mread_1.0.1             polyclip_1.10-7           GenomeInfoDbData_1.2.14   ExperimentHub_2.16.1      SparseArray_1.8.1        
 [51] xtable_1.8-4              pracma_2.4.4              S4Arrays_1.8.1            BiocFileCache_2.16.1      hms_1.1.3                
 [56] irlba_2.3.5.1             colorspace_2.1-1          filelock_1.0.3            hdf5r_1.3.12              ROCR_1.0-11              
 [61] spatstat.data_3.1-9       magrittr_2.0.3            lmtest_0.9-40             snakecase_0.11.1          later_1.4.2              
 [66] lattice_0.22-7            glmGamPoi_1.20.0          spatstat.geom_3.7-3       future.apply_1.20.0       scattermore_1.2          
 [71] XML_3.99-0.18             RcppAnnoy_0.0.22          pillar_1.11.0             nlme_3.1-168              compiler_4.5.1           
 [76] beachmat_2.24.0           RSpectra_0.16-2           stringi_1.8.7             tensor_1.5.1              plyr_1.8.9               
 [81] crayon_1.5.3              abind_1.4-8               locfit_1.5-9.12           bit_4.6.0                 sandwich_3.1-1           
 [86] multcomp_1.4-29           codetools_0.2-20          crosstalk_1.2.1           bslib_0.9.0               alabaster.ranges_1.8.0   
 [91] paletteer_1.6.0           plotly_4.11.0             mime_0.13                 splines_4.5.1             circlize_0.4.16          
 [96] fastDummies_1.7.5         dbplyr_2.5.0              sparseMatrixStats_1.20.0  blob_1.2.4                utf8_1.2.6               
[101] BiocVersion_3.21.1        listenv_0.9.1             DelayedMatrixStats_1.30.0 estimability_1.5.1        ggsignif_0.6.4           
[106] statmod_1.5.0             tzdb_0.5.0                pkgconfig_2.0.3           tools_4.5.1               cachem_1.1.0             
[111] RhpcBLASctl_0.23-42       RSQLite_2.4.2             DBI_1.2.3                 fastmap_1.2.0             scales_1.4.0             
[116] grid_4.5.1                ica_1.0-3                 sass_0.4.10               broom_1.0.8               AnnotationHub_3.16.1     
[121] coda_0.19-4.1             FNN_1.1.4.1               insight_1.5.2             ggprism_1.0.7             BiocManager_1.30.26      
[126] dotCall64_1.2             RANN_2.6.2                alabaster.schemas_1.8.0   farver_2.1.2              yaml_2.3.10              
[131] cli_3.6.5                 lifecycle_1.0.5           rsconnect_1.5.0           mvtnorm_1.3-3             backports_1.5.0          
[136] BiocParallel_1.42.1       timechange_0.3.0          gtable_0.3.6              ggridges_0.5.7            progressr_0.17.0         
[141] parallel_4.5.1            limma_3.64.3              jsonlite_2.0.0            edgeR_4.6.3               RcppHNSW_0.6.0           
[146] bit64_4.6.0-1             Rtsne_0.17                alabaster.matrix_1.8.0    spatstat.utils_3.2-2      BiocNeighbors_2.2.0      
[151] jquerylib_0.1.4           alabaster.se_1.8.0        spatstat.univar_3.1-7     R.utils_2.13.0            lazyeval_0.2.2           
[156] alabaster.base_1.8.1      htmltools_0.5.8.1         sctransform_0.4.2         rappdirs_0.3.3            glue_1.8.0               
[161] httr2_1.2.1               XVector_0.48.0            mclust_6.1.2              ks_1.15.1                 gridExtra_2.3            
[166] igraph_2.1.4              R6_2.6.1                  labeling_0.4.3            cluster_2.1.8.1           pkgload_1.4.0            
[171] Rhdf5lib_1.30.0           mcprogress_0.1.1          DelayedArray_0.34.1       tidyselect_1.2.1          vipor_0.4.7              
[176] maps_3.4.3                xml2_1.3.8                rsvd_1.0.5                KernSmooth_2.23-26        data.table_1.17.8        
[181] htmlwidgets_1.6.4         rlang_1.2.0               spatstat.sparse_3.1-0     spatstat.explore_3.8-0    remotes_2.5.0            
[186] rentrez_1.2.4             Cairo_1.7-0               beeswarm_0.4.0           
```

Package versions can be frozen using:

```r
renv::snapshot()
```

Using the same input data and package versions should reproduce all analyses and figures.

---

# Data Availability

The pipeline was developed using publicly available melanoma single-cell RNA-seq data from:

Yost KE *et al.* Nature Medicine (2019).


---

# Code Availability

All scripts required to reproduce the analyses described in the associated manuscript are provided in this repository. The software depends exclusively on publicly available R packages and is distributed under an open-source license.


---

# License

This project is distributed under the MIT License.
