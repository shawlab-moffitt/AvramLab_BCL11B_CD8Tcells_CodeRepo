library(Seurat)
library(harmony)
library(scCustomize)
library(ggplot2)
Yost_Tcell_combined_seurat <- readRDS("/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/Yost_Tcell_combined_seurat_20260707.rds")

#Yost_combined_seurat <- SplitObject(Yost_combined_seurat, split.by = "orig.ident")

# Human
Yost_Tcell_combined_seurat[["percent.mt"]] <- PercentageFeatureSet(Yost_Tcell_combined_seurat, pattern = "^MT-")
Yost_Tcell_combined_seurat[["percent.ribo"]] <- PercentageFeatureSet(Yost_Tcell_combined_seurat, pattern = "^RPL|^RPS")
Yost_Tcell_combined_seurat[["percent.hb"]] <- PercentageFeatureSet(Yost_Tcell_combined_seurat, pattern = "^HB[AB]")


Yost_Tcell_combined_seurat$All <- "All"

VlnPlot(
  Yost_Tcell_combined_seurat,
  group.by = "All",
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

FeatureScatter(Yost_Tcell_combined_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Yost_Tcell_combined_seurat <- subset(
  Yost_Tcell_combined_seurat,
  subset =
    nFeature_RNA > 200 
    #nFeature_RNA < 8000 &
)

DefaultAssay(Yost_Tcell_combined_seurat) <- "RNA"

Yost_Tcell_combined_seurat <- NormalizeData(Yost_Tcell_combined_seurat)

Yost_Tcell_combined_seurat <- FindVariableFeatures(
  Yost_Tcell_combined_seurat,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

hvgs <- VariableFeatures(Yost_Tcell_combined_seurat)

# Important: scale only HVGs, not all genes
Yost_Tcell_combined_seurat <- ScaleData(
  Yost_Tcell_combined_seurat,
  assay = "RNA",
  #features = hvgs,
  verbose = FALSE
)


Yost_Tcell_combined_seurat <- RunPCA(
  Yost_Tcell_combined_seurat,
  assay = "RNA",
  reduction.name = "rna_pca_v1",
  reduction.key = "RNAPCV1_",
  verbose = FALSE
)

Yost_Tcell_combined_seurat@meta.data
#unique(integrated_darwin_t_cell_seurat$genotype)
Yost_Tcell_combined_seurat <- RunHarmony(
  object = Yost_Tcell_combined_seurat,
  group.by.vars = "donor",
  reduction.use = "rna_pca_v1",
  dims.use = 1:30,
  reduction.save = "harmony"
)

Yost_Tcell_combined_seurat@reductions
# UMAP using SCT PCA
Yost_Tcell_combined_seurat <- RunUMAP(
  Yost_Tcell_combined_seurat,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 50,
  min.dist = 0.05,
  reduction.name = "rna_umap_v1",
  reduction.key = "RNAUMAPV1_"
)
Yost_Tcell_combined_seurat@meta.data$patient = Yost_Tcell_combined_seurat@meta.data$donor
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", group.by = "donor")
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", group.by = "treatment")
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", group.by = "outcome")
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", group.by = "cell type")
FeaturePlot(Yost_Tcell_combined_seurat, features = c("CD3E", "CD8A", "CD4", "MKI67"))
FeaturePlot(Yost_Tcell_combined_seurat, features = c("TCF7", "CCR7", "PRF1", "PDCD1"))
Yost_Tcell_combined_seurat@meta.data
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", group.by = "donor")

# Neighbors using SCT PCA
Yost_Tcell_combined_seurat <- FindNeighbors(
  Yost_Tcell_combined_seurat,
  reduction = "harmony",
  dims = 1:30,
  graph.name = "rna_nn_v1",
  verbose = FALSE
)

# Clustering using SCT neighbor graph
Yost_Tcell_combined_seurat <- FindClusters(
  Yost_Tcell_combined_seurat,
  graph.name = "rna_nn_v1",
  cluster.name = "rna_clusters_orig_v1",
  verbose = FALSE
)

# Lower-resolution clustering
Yost_Tcell_combined_seurat <- FindClusters(
  Yost_Tcell_combined_seurat,
  graph.name = "rna_nn_v1",
  cluster.name = "rna_clusters_v1_0.2",
  resolution = 0.2,
  verbose = FALSE
)
# Lower-resolution clustering
Yost_Tcell_combined_seurat <- FindClusters(
  Yost_Tcell_combined_seurat,
  graph.name = "rna_nn_v1",
  cluster.name = "rna_clusters_v1_0.1",
  resolution = 0.1,
  verbose = FALSE
)
# Lower-resolution clustering
Yost_Tcell_combined_seurat <- FindClusters(
  Yost_Tcell_combined_seurat,
  graph.name = "rna_nn_v1",
  cluster.name = "rna_clusters_v1_0.3",
  resolution = 0.3,
  verbose = FALSE
)
# Lower-resolution clustering
Yost_Tcell_combined_seurat <- FindClusters(
  Yost_Tcell_combined_seurat,
  graph.name = "rna_nn_v1",
  cluster.name = "rna_clusters_v1_0.5",
  resolution = 0.5,
  verbose = FALSE
)
# Lower-resolution clustering
Yost_Tcell_combined_seurat <- FindClusters(
  Yost_Tcell_combined_seurat,
  graph.name = "rna_nn_v1",
  cluster.name = "rna_clusters_v1_1.0",
  resolution = 1.0,
  verbose = FALSE
)




gene_groups1 <- list(
  "Immune" = c("PTPRC"),
  "T" = c("CD3D", "CD3E", "BCL11B", "TCF12", "TCF3", "TRAC"),
  "Myeloid" = c("LYZ", "ITGAX", "ITGAM"),
  "CD8" = c("CD8A", "CD8B", "RUNX3"),
  "CD4" = c("CD4", "ZBTB7B"),
  "gdT" = c("TRDV2", "TRGV9", "TRDV1", "TRDV3", "CX3CR1", "TRDC", "TRGC1"),
  "Mito" = c("MT-ND4L", "MT-ND5", "MT-CYB"),
  "Th1Th2h17Tfh" = c("TBX21", "GATA3", "RORC", "BCL6"),
  "nT" = c("PECAM1"),
  "cmT" = c("FAS"),
  "Treg" = c("IL2RA", "FOXP3", "CTLA4", "TNFRSF18"),
  "NK" = c("NCAM1", "NCR1", "KLRK1", "KLRD1", "KLRC1"),
  "Cycling" = c("MKI67", "CDKN2A"),
  "Effector" = c("PRF1", "GZMA", "GZMB", "GZMK", "XCL1", "XCL2", "CRTAM"),
  "Exhaustion" = c("TOX", "PDCD1", "LAG3", "TIGIT", "HAVCR2"),
  "Residency" = c("ITGAE", "ITGA1", "CD69", "ZNF683", "PRDM1"),
  "Trafficking" = c("CXCR6", "CXCR4", "S1PR1", "KLF2", "LTB", "S1PR4"),
  "Stemness" = c("TCF7", "SELL", "IL7R", "SLAMF6", "CCR7", "CCR5", "CXCR3", "BCL2", "EOMES", "BACH2", "LEF1")
)

dotplot <- DotPlot(Yost_Tcell_combined_seurat, features = gene_groups1, group.by = "rna_clusters_v1_0.5") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "DotPlot with Grouped Genes")

dotplot

DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", label=T, group.by = "donor")
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", label=T, group.by = "rna_clusters_v1_0.5")
table(Yost_Tcell_combined_seurat$donor, Yost_Tcell_combined_seurat$rna_clusters_v1_0.3)
table(Yost_Tcell_combined_seurat$rna_clusters_v1_1.0, Yost_Tcell_combined_seurat$patient)



Yost_Tcell_combined_seurat$CD45RO_Proportion = Yost_Tcell_combined_seurat$CD45RO_Junction1_umis / (Yost_Tcell_combined_seurat$CD45RO_Junction1_umis + Yost_Tcell_combined_seurat$CD45RA_Ex3to4_umis + Yost_Tcell_combined_seurat$CD45RB_Ex3to5_umis + Yost_Tcell_combined_seurat$CD45RC_Ex3to6_umis)
Yost_Tcell_combined_seurat$CD45RA_Proportion = Yost_Tcell_combined_seurat$CD45RA_Ex3to4_umis / (Yost_Tcell_combined_seurat$CD45RO_Junction1_umis + Yost_Tcell_combined_seurat$CD45RA_Ex3to4_umis + Yost_Tcell_combined_seurat$CD45RB_Ex3to5_umis + Yost_Tcell_combined_seurat$CD45RC_Ex3to6_umis)
Yost_Tcell_combined_seurat$CD45RB_Proportion = Yost_Tcell_combined_seurat$CD45RB_Ex3to5_umis / (Yost_Tcell_combined_seurat$CD45RO_Junction1_umis + Yost_Tcell_combined_seurat$CD45RA_Ex3to4_umis + Yost_Tcell_combined_seurat$CD45RB_Ex3to5_umis + Yost_Tcell_combined_seurat$CD45RC_Ex3to6_umis)
Yost_Tcell_combined_seurat$CD45RC_Proportion = Yost_Tcell_combined_seurat$CD45RC_Ex3to6_umis / (Yost_Tcell_combined_seurat$CD45RO_Junction1_umis + Yost_Tcell_combined_seurat$CD45RA_Ex3to4_umis + Yost_Tcell_combined_seurat$CD45RB_Ex3to5_umis + Yost_Tcell_combined_seurat$CD45RC_Ex3to6_umis)

Yost_Tcell_combined_seurat$CD74withoutExon7_Proportion = Yost_Tcell_combined_seurat$CD74withoutExon7_umis / (Yost_Tcell_combined_seurat$CD74withoutExon7_umis + Yost_Tcell_combined_seurat$CD74withExon7Up_umis)
Yost_Tcell_combined_seurat$CD74withExon7_Proportion = Yost_Tcell_combined_seurat$CD74withExon7Up_umis / (Yost_Tcell_combined_seurat$CD74withoutExon7_umis + Yost_Tcell_combined_seurat$CD74withExon7Up_umis)
Yost_Tcell_combined_seurat$CCL5withoutExon3_Proportion = Yost_Tcell_combined_seurat$CCL5withoutExon3_umis / (Yost_Tcell_combined_seurat$CCL5withoutExon3_umis + Yost_Tcell_combined_seurat$CCL5withExon3up_umis)
Yost_Tcell_combined_seurat$CCL5withExon3_Proportion = Yost_Tcell_combined_seurat$CCL5withExon3up_umis / (Yost_Tcell_combined_seurat$CCL5withoutExon3_umis + Yost_Tcell_combined_seurat$CCL5withExon3up_umis)
Yost_Tcell_combined_seurat$XBP1s_Proportion = Yost_Tcell_combined_seurat$XBP1s_umis / (Yost_Tcell_combined_seurat$XBP1s_umis + Yost_Tcell_combined_seurat$XBP1ex1to2_umis)
Yost_Tcell_combined_seurat$XBP1noframeshift_Proportion = Yost_Tcell_combined_seurat$XBP1ex1to2_umis / (Yost_Tcell_combined_seurat$XBP1s_umis + Yost_Tcell_combined_seurat$XBP1ex1to2_umis)


Idents(Yost_Tcell_combined_seurat) <- "rna_clusters_v1_0.5"
CD8_Yost_Tcell_combined_seurat <- subset(Yost_Tcell_combined_seurat, idents = c("0", "4", "6", "7"))
CD8_Yost_Tcell_combined_seurat@meta.data$outcome



CD8_Pre_Yost_Tcell_combined_seurat <- subset(
  Yost_Tcell_combined_seurat,
  idents = "0",
  subset = treatment == "pre"
)

# Post-treatment CD8 cells
CD8_Post_Yost_Tcell_combined_seurat <- subset(
  Yost_Tcell_combined_seurat,
  idents = "0",
  subset = treatment == "post"
)

# ---- write out seurat object to file ----
saveRDS(
  CD8_Yost_Tcell_combined_seurat,
  file = "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/CD8_Yost_Tcell_combined_seurat_20260720.rds"
)

saveRDS(
  CD8_Pre_Yost_Tcell_combined_seurat,
  file = "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/CD8_Pre_Yost_Tcell_combined_seurat_20260720.rds"
)

saveRDS(
  CD8_Post_Yost_Tcell_combined_seurat,
  file = "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/CD8_Post_Yost_Tcell_combined_seurat_20260720.rds"
)

saveRDS(
  Yost_Tcell_combined_seurat,
  file = "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/Yost_Tcell_combined_seurat_20260720.rds"
)

# ---- generate violin plot comparing outcome in PRE treatment ----

genes <- c(
 "BCL11B","TCF7", "BACH1", "SATB1", "CXCR3", "TOX", "LAG3",
  "NR4A1", "NR4A2", "PDCD1", "CTLA4", "HAVCR2", "TIGIT",
  "PRF1", "GZMB", "TNF"
)

# Check genes present
genes_present <- genes[genes %in% rownames(CD8_Pre_Yost_Tcell_combined_seurat)]
genes_missing <- setdiff(genes, genes_present)
print(genes_missing)

# Check outcome labels
table(CD8_Pre_Yost_Tcell_combined_seurat$outcome)

# Restrict to CR and NR only
CD8_Pre_CRNR <- subset(
  CD8_Pre_Yost_Tcell_combined_seurat,
  subset = outcome %in% c("CR", "NR")
)

# Ensure desired order
CD8_Pre_CRNR$outcome <- factor(
  CD8_Pre_CRNR$outcome,
  levels = c("CR", "NR")
)


# ---- PRE Treatment ----

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggbeeswarm)

CD8_Pre_CRNR$outcome <- factor(CD8_Pre_CRNR$outcome, levels = c("NR", "CR"))

genes <- c(
  "BCL11B","BACH2", "SATB1", "PRF1", "GZMB"
  #, "CXCR3", "TOX", "LAG3",
  #"NR4A1", "NR4A2", "PDCD1", "CTLA4", "HAVCR2", "TIGIT",
  #"PRF1", "GZMB"
)
genes_present <- genes[genes %in% rownames(CD8_Pre_CRNR)]

# Function to calculate geometric mean
geo_mean <- function(x) {
  exp(mean(log(x + 1))) - 1
}

# Calculate geometric means + Wilcoxon p-values
pvalue_table <- lapply(genes_present, function(gene) {
  df <- FetchData(CD8_Pre_CRNR, vars = c(gene, "outcome"))
  colnames(df) <- c("expression", "outcome")
  df <- na.omit(df)
  
  wt <- wilcox.test(expression ~ outcome, data = df)
  
  cr_gmean <- mean(df$expression[df$outcome == "CR"])
  nr_gmean <- mean(df$expression[df$outcome == "NR"])
  
  p_text <- ifelse(
    wt$p.value < 0.001,
    "<0.001",
    sprintf("%.2f", wt$p.value)
  )
  
  data.frame(
    Gene = gene,
    CR_GMean = cr_gmean,
    NR_GMean = nr_gmean,
    Pvalue = wt$p.value,
    Plot_Title = paste0(
      gene,

      "(p=", p_text, ")\n",
      "(NR Mean=", round(nr_gmean, 2), "; ",
      "  CR Mean=", round(cr_gmean, 2), ")"
      
    )
  )
}) %>%
  bind_rows()

# Generate individual violin plots
plot_list <- lapply(seq_along(genes_present), function(i) {
  gene <- genes_present[i]
  title_i <- pvalue_table$Plot_Title[pvalue_table$Gene == gene]
  
  VlnPlot(
    object = CD8_Pre_CRNR,
    features = gene,
    group.by = "outcome",
    pt.size = 0,
    cols = c("#ff8686", "white")
  ) +
  #  geom_boxplot(
  #    width = 0.05,
  #    fill = "white",
  #    color = "black",
  #    outlier.shape = NA,
  #    alpha = 0.4
  #  ) +
    geom_quasirandom(
      width = 0.2,
      varwidth = TRUE,
      size = 0.15,
      alpha = 0.3,
      color = "black"
    ) +
    ggtitle(title_i) +
    xlab(NULL) +
    ylab(NULL) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, hjust = 0.5),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank()
    )
})

# Stack vertically, similar to Stacked_VlnPlot
p <- wrap_plots(plot_list, ncol = 1)

p
ggsave(
  filename = "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/Yost_et_al_CD8_Pre_CR_vs_NR_stacked_violin.svg",
  plot = p,
  device = svglite,
  width = 4,
  height = length(genes_present) * 1.2,
  units = "in"
)

# Generate and save one SVG per gene
plot_list <- lapply(seq_along(genes_present), function(i) {
  gene <- genes_present[i]
  title_i <- pvalue_table$Plot_Title[pvalue_table$Gene == gene]
  
  p <- VlnPlot(
    object = CD8_Pre_CRNR,
    features = gene,
    group.by = "outcome",
    pt.size = 0,
    cols = c("#ff8686", "white")
  ) +
    geom_quasirandom(
      width = 0.2,
      varwidth = TRUE,
      size = 0.15,
      alpha = 0.3,
      color = "black"
    ) +
    ggtitle(title_i) +
    xlab(NULL) +
    ylab(gene) +      # Show gene name on the y-axis
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, hjust = 0.5),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold")
    )
  
  # Save SVG per gene
  ggsave(
    filename = file.path(paste0("/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/Yost_et_al_CD8_Pre_", gene, "_CR_vs_NR_violin.svg")),
    plot = p,
    device = "svg",
    width = 4,
    height = 4
  )
  
  return(p)
})

# ---- updated one  ----
# Generate and save one SVG and one CSV per gene
plot_list <- lapply(seq_along(genes_present), function(i) {
  
  gene <- genes_present[i]
  title_i <- pvalue_table$Plot_Title[pvalue_table$Gene == gene]
  
  # Extract raw values for this gene
  expr_df <- FetchData(
    object = CD8_Pre_CRNR,
    vars = c(gene, "outcome")
  )
  
  expr_df$Cell <- rownames(expr_df)
  
  # Reorder columns
  expr_df <- expr_df[, c("Cell", "outcome", gene)]
  
  # Save raw values
  write.csv(
    expr_df,
    file = file.path(
      "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal",
      paste0("Yost_et_al_CD8_Pre_", gene, "_CR_vs_NR_raw_values.csv")
    ),
    row.names = FALSE
  )
  
  # Generate violin plot
  p <- VlnPlot(
    object = CD8_Pre_CRNR,
    features = gene,
    group.by = "outcome",
    pt.size = 0,
    cols = c("#ff8686", "white")
  ) +
    geom_quasirandom(
      width = 0.2,
      varwidth = TRUE,
      size = 0.15,
      alpha = 0.3,
      color = "black"
    ) +
    ggtitle(title_i) +
    xlab(NULL) +
    ylab(gene) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, hjust = 0.5),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold")
    )
  
  # Save SVG
  ggsave(
    filename = file.path(
      "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal",
      paste0("Yost_et_al_CD8_Pre_", gene, "_CR_vs_NR_violin_v2.svg")
    ),
    plot = p,
    device = "svg",
    width = 4,
    height = 4
  )
  
  return(p)
})

# ---- POST Treatment ----

# Check outcome labels
table(CD8_Post_Yost_Tcell_combined_seurat$outcome)

# Restrict to CR and NR only
CD8_Post_CRNR <- subset(
  CD8_Post_Yost_Tcell_combined_seurat,
  subset = outcome %in% c("CR", "NR")
)

# Ensure desired order
CD8_Post_CRNR$outcome <- factor(
  CD8_Post_CRNR$outcome,
  levels = c("CR", "NR")
)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggbeeswarm)
CD8_Post_CRNR$outcome <- factor(CD8_Post_CRNR$outcome, levels = c("CR", "NR"))

genes <- c(
  "BCL11B","TCF7", "BACH2", "SATB1", "CXCR3", "TOX", "LAG3",
  "NR4A1", "NR4A2", "PDCD1", "CTLA4", "HAVCR2", "TIGIT",
  "PRF1", "GZMB"
)
genes_present <- genes[genes %in% rownames(CD8_Post_CRNR)]

# Function to calculate geometric mean
geo_mean <- function(x) {
  exp(mean(log(x + 1))) - 1
}

# Calculate geometric means + Wilcoxon p-values
pvalue_table <- lapply(genes_present, function(gene) {
  df <- FetchData(CD8_Post_CRNR, vars = c(gene, "outcome"))
  colnames(df) <- c("expression", "outcome")
  df <- na.omit(df)
  
  wt <- wilcox.test(expression ~ outcome, data = df)
  
  cr_gmean <- geo_mean(df$expression[df$outcome == "CR"])
  nr_gmean <- geo_mean(df$expression[df$outcome == "NR"])
  
  p_text <- ifelse(
    wt$p.value < 0.001,
    "<0.001",
    sprintf("%.2f", wt$p.value)
  )
  
  data.frame(
    Gene = gene,
    CR_GMean = cr_gmean,
    NR_GMean = nr_gmean,
    Pvalue = wt$p.value,
    Plot_Title = paste0(
      gene,
      "  CR GM=", round(cr_gmean, 2),
      " | NR GM=", round(nr_gmean, 2),
      " | p=", p_text
    )
  )
}) %>%
  bind_rows()

# Generate individual violin plots
plot_list <- lapply(seq_along(genes_present), function(i) {
  gene <- genes_present[i]
  title_i <- pvalue_table$Plot_Title[pvalue_table$Gene == gene]
  
  VlnPlot(
    object = CD8_Post_CRNR,
    features = gene,
    group.by = "outcome",
    pt.size = 0,
    cols = c("white", "#ff8686")
  ) +
    #  geom_boxplot(
    #    width = 0.05,
    #    fill = "white",
    #    color = "black",
    #    outlier.shape = NA,
    #    alpha = 0.4
    #  ) +
    geom_quasirandom(
      width = 0.2,
      varwidth = TRUE,
      size = 0.15,
      alpha = 0.3,
      color = "black"
    ) +
    ggtitle(title_i) +
    xlab(NULL) +
    ylab(NULL) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 7, hjust = 0.5),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 5),
      axis.title = element_blank()
    )
})

# Stack vertically, similar to Stacked_VlnPlot
p <- wrap_plots(plot_list, ncol = 1)

p

library(svglite)

ggsave(
  filename = "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/Yost_et_al_CD8_Post_CR_vs_NR_stacked_violin.svg",
  plot = p,
  device = svglite,
  width = 4,
  height = length(genes_present) * 1.2,
  units = "in"
)
# ---- END POST Treatment plot ----


colnames(Yost_Tcell_combined_seurat@meta.data)
juncscope_meta = Yost_Tcell_combined_seurat@meta.data
# Create matrix of CD45 isoform proportions
cd45_props <- Yost_Tcell_combined_seurat@meta.data[, c(
  "CD45RO_Proportion",
  "CD45RA_Proportion",
  "CD45RB_Proportion",
  "CD45RC_Proportion"
)]

# Assign highest isoform for each cell
Yost_Tcell_combined_seurat@meta.data$CD45_Dominant_Isoform <- colnames(cd45_props)[
  max.col(cd45_props, ties.method = "first")
]

# Optional: remove "_Proportion" suffix
Yost_Tcell_combined_seurat@meta.data$CD45_Dominant_Isoform <- gsub(
  "_Proportion",
  "",
  Yost_Tcell_combined_seurat@meta.data$CD45_Dominant_Isoform
)

Yost_Tcell_combined_seurat$response = "No"
Yost_Tcell_combined_seurat@meta.data[which(Yost_Tcell_combined_seurat$`pathological response` == "CR"), "response"] = "Yes"
Yost_Tcell_combined_seurat@meta.data[which(Yost_Tcell_combined_seurat$`pathological response` == "PR"), "response"] = "Yes"
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", label=T, group.by = "rna_clusters_v1_1.0")
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", label=T, group.by = "patient")
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", label=T, group.by = "CD45_Dominant_Isoform")

DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", label=T, group.by = "treatment")
DimPlot(Yost_Tcell_combined_seurat, reduction = "rna_umap_v1", label=T, group.by = "pathological response")

FeaturePlot(Yost_Tcell_combined_seurat, feature = c("CD45RA_Proportion", "CD45RO_Proportion", "CD74withoutExon7_Proportion", "CD74withExon7_Proportion", 
                                                    "CCL5withoutExon3_Proportion", "CCL5withExon3_Proportion", "XBP1s_Proportion", "XBP1noframeshift_Proportion"))
table(Yost_Tcell_combined_seurat$rna_clusters_v1_1.0, Yost_Tcell_combined_seurat$patient)
table(Yost_Tcell_combined_seurat$response, Yost_Tcell_combined_seurat$CD45_Dominant_Isoform)


CD8_Yost_Tcell_combined_seurat <- subset(Yost_Tcell_combined_seurat, idents = c("0", "2", "4", "5", "6", "7", "8", "9", "11", "12", "15", "17", "18", "19"))
table(CD8_Yost_Tcell_combined_seurat$response, CD8_Yost_Tcell_combined_seurat$CD45_Dominant_Isoform)


tab <- table(
  Yost_Tcell_combined_seurat$rna_clusters_v1_1.0,
  Yost_Tcell_combined_seurat$patient
)

prop.tab <- prop.table(tab, margin = 2)

round(prop.tab, 3)
prop.df <- as.data.frame(prop.tab)

colnames(prop.df) <- c("Cluster", "Patient", "Proportion")
Yost_Tcell_combined_seurat$Histology
prop.df
patient.response <- unique(
  Yost_Tcell_combined_seurat@meta.data[, c("patient", "pathological response", "treatment", "")]
)

prop.df <- left_join(
  prop.df,
  patient.response,
  by = c("Patient" = "patient")
)

ggplot(
  subset(prop.df, Cluster == "2"),
  aes(x = `pathological response`, y = Proportion)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3) +
  stat_compare_means(method = "kruskal.test") +
  theme_classic()

# ---- compare Response vs No Response ----

patient.df <- CD8_Yost_Tcell_combined_seurat@meta.data %>%
  group_by(patient, treatment) %>%
  summarize(
    mean_CD45RA = mean(CD45RA_Proportion, na.rm = TRUE),
    median_CD45RA = median(CD45RA_Proportion, na.rm = TRUE),
    mean_CD45RB = mean(CD45RB_Proportion, na.rm = TRUE),
    median_CD45RB = median(CD45RB_Proportion, na.rm = TRUE),
    mean_CD45RO = mean(CD45RO_Proportion, na.rm = TRUE),
    median_CD45RO = median(CD45RO_Proportion, na.rm = TRUE),
    
    n_cells = n(),
    .groups = "drop"
  )

patient.df

ggplot(patient.df,
       aes(x = treatment,
           y = mean_CD45RA,
           color = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_jitter(width = 0.1, size = 3) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  theme_classic() +
  labs(
    x = "treatment",
    y = "Mean CD45RA Proportion"
  )

ggplot(patient.df,
       aes(response, mean_CD45RA, fill = response)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3) +
  stat_compare_means(method = "wilcox.test") +
  theme_classic()

