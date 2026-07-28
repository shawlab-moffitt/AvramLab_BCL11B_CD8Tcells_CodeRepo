library(Seurat)
library(harmony)
library(scCustomize)
library(ggplot2)
Yost_Tcell_combined_seurat <- readRDS("Yost_Tcell_combined_seurat_20260707.rds")


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

#saveRDS(
#  Yost_Tcell_combined_seurat,
#  file = "/Volumes/Shaw_Timothy/TcellSubtypes/external_data/melanoma/Yostetal/Yost_Tcell_combined_seurat_20260720.rds"
#)

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


