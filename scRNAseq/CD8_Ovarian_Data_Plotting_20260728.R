# Load the libraries ------------------------------------------------------


library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(presto)
library(ggrepel)

# Load the .rds files -----------------------------------------------------
seu_obj = readRDS("CD8_Ovarian_TILS_Bcl11bKOvsWT_Annalysis_Finalised_20262707.Rds")


custom_colors <- c(
  "Tpex1"   = "#00B050",  
  "Tinex1"  = "#FF6F61",  
  "Tinex3" = "#FF00FF",  
  "Tpex2"   = "#00CED1", 
  "Tinex2"  = "#BF9000",  
  "Ttex"    = "#000080"  
)

# Subset KO and WT cells
seu_ko <- subset(seu_obj, subset = condition == "KO")
seu_wt <- subset(seu_obj, subset = condition == "WT")

# Get same PCA axis limits from full object
pca_embed <- Embeddings(seu_obj, reduction = "pca_RNA")

xlims <- range(pca_embed[, 3], na.rm = TRUE)
ylims <- range(pca_embed[, 1], na.rm = TRUE)

# Plot again with same axis limits
PCA <- DimPlot(
  seu_obj,
  reduction = "pca_RNA",
  group.by = "renamed_cluster",
  dims = c(3, 1),
  pt.size = 2.8,
  cols = custom_colors,
  label = FALSE
) +
  xlim(xlims) +
  ylim(ylims) +
  ggtitle("WT")+
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold", size = 19),
    axis.text = element_text(face = "bold", size = 18),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

PCA

set.seed(123)  # For reproducibility


WT <- DimPlot(
  seu_wt,
  reduction = "pca_RNA",
  group.by = "renamed_cluster",
  dims = c(3, 1),
  pt.size = 2.7,
  cols = custom_colors,
  label = FALSE
) +
  xlim(xlims) +
  ylim(ylims) +
  ggtitle("WT")+
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold", size = 26),
    axis.text = element_text(face = "bold", size = 22),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22)
  )

# Plot KO
KO <- DimPlot(
  seu_ko,
  reduction = "pca_RNA",
  group.by = "renamed_cluster",
  dims = c(3, 1),
  cols = custom_colors,
  pt.size = 2.5,
  label = F
) +
  xlim(xlims) +
  ylim(ylims) +
  ggtitle("KO")+
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold", size = 26),
    axis.text = element_text(face = "bold", size = 22),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22)
  )

WT
KO

#ggsave("~/Bcl11b_KO_CD8_OVARIAN_Mice/Fig_PCA_UMAP/CD8_Ovarian_PCAbased_Integrated_20260706.svg", plot = PCA, width = 11, height = 8, dpi = 400)
#ggsave("~/Bcl11b_Project/CD8_Ovarian_mouse/Figures/Version1/CD8_Ovarian_PCAbased_Condition_20250915_2.svg", plot = cond, width = 10, height = 8, dpi = 400)

#ggsave("~/Bcl11b_KO_CD8_OVARIAN_Mice/Fig_PCA_UMAP/CD8_Ovarian_PCAbased_WT_20260706.svg", plot = WT, width = 15, height = 10, dpi = 400)
#ggsave("~/Bcl11b_KO_CD8_OVARIAN_Mice/Fig_PCA_UMAP/CD8_Ovarian_PCAbased_KO_20260706.svg", plot = KO, width = 10, height = 8, dpi = 400)


# Umaps -------------------------------------------------------------------


WT = DimPlot(seu_wt, group.by = "renamed_cluster",reduction = "umap" ,cols = custom_colors,label = FALSE, pt.size = 2.5) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold", size = 19),
    axis.text = element_text(face = "bold", size = 18),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    #legend.text = element_text(size = 16),
    #legend.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )+
  ggtitle("WT")

WT



KO = DimPlot(seu_ko, group.by = "renamed_cluster",reduction = "umap" ,cols = custom_colors,label = FALSE, pt.size = 2.5) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold", size = 19),
    axis.text = element_text(face = "bold", size = 18),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    #legend.text = element_text(size = 14),
    #legend.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )+
  ggtitle("KO")

KO



#  ================== Calculate Cluster Proportions  ==================

# 1. Extract metadata
metadata <- seu_obj@meta.data

# 2. Calculate cluster proportions
cluster_counts <- metadata %>%
  group_by(condition, renamed_cluster) %>%
  summarise(n = n()) %>%
  group_by(condition) %>%
  mutate(freq = n / sum(n) * 100)

# 2. Set condition order
position = position_dodge(width = 0.8)

cluster_counts$condition <- factor(
  cluster_counts$condition,
  levels = c("WT", "KO")
)

cluster_counts$renamed_cluster <- factor(
  cluster_counts$renamed_cluster,
  levels = c(
    "Tpex1",
    "Tpex2",
    "Tinex1",
    "Tinex2",
    "Tinex3",
    "Ttex"
    
  )
)

# 3. Plot
prop_plot <- ggplot(cluster_counts, aes(x = renamed_cluster, y = freq, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("WT" = "#6eb2ff", "KO" = "#FF6F61")) +  # Cyan for WT, Coral for KO
  labs(x = NULL, y = "Cluster proportion [%]") +
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

# 4. Show plot
prop_plot


# 4. Plot with black x and y axis lines
prop_plot <- ggplot(cluster_counts, aes(x = renamed_cluster, y = freq, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("WT" = "#6eb2ff", "KO" = "#FF6F61")) +  # Cyan for WT, Coral for KO
  labs(x = NULL, y = "Cluster proportion [%]") +
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold", angle = 0, hjust = 1),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 22, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 1, linewidth = 1.5)   # << add X and Y axis lines
  ) +
  guides(fill = guide_legend(override.aes = list(size = 8)))

# 5. View the final plot
prop_plot




#  ================== Dotplot  ==================

Idents(seu_obj) = "condition_Annotation"

# Define ordered gene list by functional category
genes_ordered <- c(
  # Stemness
  "Tcf7", "Lef1", "Bach2", "Satb1", "Sell", "Runx1", "Il7r", "Icos", "Bcl2",
  # Exhaustion
  "Nr4a1", "Nr4a2", "Pdcd1", "Lag3", "Havcr2", "Tigit", "Ctla4",
  # Effector
  "Gzma","Gzmb", "Gzmc", "Gzmd", "Gzme", "Gzmf", "Gzmg", "Prf1",
  # Innate
  "Ncr1", "Klrd1", "Klrb1c", "Klre1", "Klrg1", "Klrk1", "Ly6c2", "Fcer1g", "Fcgr3"
)

# Reverse gene order for coord_flip (so stemness at top)
genes_reversed <- rev(genes_ordered)

# Set desired cluster order: alternating WT and KO by number
new_order <- c("WT_Tpex1", 
               "KO_Tpex1",
               "WT_Tpex2",
               "KO_Tpex2" ,
               "WT_Tinex1", 
               "KO_Tinex1", 
               "WT_Tinex2",
               "KO_Tinex2",
               "WT_Tinex3", 
               "KO_Tinex3",
               "WT_Ttex", 
               "KO_Ttex")

# Ensure the metadata column used for grouping exists and is correctly ordered
seu_obj$condition_Annotation <- factor(seu_obj$condition_Annotation, levels = new_order)

# Generate the DotPlot
P = DotPlot(seu_obj, features = genes_reversed, group.by = "condition_Annotation") +
  scale_size(range = c(3, 10)) +
  scale_color_gradient2(
    low = "#1f78b4", mid = "white", high = "#b30000", midpoint = 0
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 14)
  ) +
  coord_flip()

P




# Ucell Signature Score Calculation and Visualization---------------------------------------


# Load required libraries
library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)

# 1. Define stemness gene set
stemness_genes <- c("Tcf7", "Lef1", "Bach2", "Satb1", "Sell", "Runx1", "Il7r", "Icos", "Bcl2")
Exhaustion_genes = c("Nr4a1", "Nr4a2", "Tox2" , "Pdcd1", "Lag3", "Havcr2", "Tigit", "Ctla4")
Innate_genes = c("Ncr1", "Klrd1", "Klrb1c", "Klre1", "Klrg1", "Klrk1", "Ly6c2", "Fcer1g", "Fcgr3")
Effector_genes = c("Gzmb", "Gzmc", "Gzmd", "Gzme", "Gzmf", "Gzmg", "Prf1")

# 2. Extract normalized expression matrix (from RNA assay)
expr_mat <- GetAssayData(seu_obj, assay = "RNA", layer = "data")  # Updated for Seurat v5


# 3. Prepare signature list
signatures_stem <- list(Stemness = stemness_genes)
signatures_exh <- list(Exhaustion = Exhaustion_genes)
signatures_inn <- list(Innate = Innate_genes)
signatures_eff <- list(Effector = Effector_genes)

# 4. Calculate UCell scores
u_scores_stem <- ScoreSignatures_UCell(expr_mat, features = signatures_stem)
u_scores_exh <- ScoreSignatures_UCell(expr_mat, features = signatures_exh)
u_scores_eff <- ScoreSignatures_UCell(expr_mat, features = signatures_eff)
u_scores_inn <- ScoreSignatures_UCell(expr_mat, features = signatures_inn)

# 5. Assign scores to Seurat metadata (ensure rownames match)
u_scores_stem <- as.data.frame(u_scores_stem)
u_scores_exh <- as.data.frame(u_scores_exh)
u_scores_eff <- as.data.frame(u_scores_eff)
u_scores_inn <- as.data.frame(u_scores_inn)

seu_obj@meta.data$Stemness_UCell <- u_scores_stem[rownames(seu_obj@meta.data), "Stemness_UCell"]
seu_obj@meta.data$Exhaustion_UCell <- u_scores_exh[rownames(seu_obj@meta.data), "Exhaustion_UCell"]
seu_obj@meta.data$Innate_UCell <- u_scores_inn[rownames(seu_obj@meta.data), "Innate_UCell"]
seu_obj@meta.data$Effector_UCell <- u_scores_eff[rownames(seu_obj@meta.data), "Effector_UCell"]


# 6. Define desired order of groups
custom_order <- c(
  "WT_Tpex1", "KO_Tpex1",
  "WT_Tpex2", "KO_Tpex2",
  "WT_Tinex1", "KO_Tinex1",
  "WT_Tinex2", "KO_Tinex2",
  "WT_Tinex3", "KO_Tinex3",
  "WT_Ttex", "KO_Ttex"
)

# 7. Prepare data frame for plotting
plot_df <- seu_obj@meta.data %>%
  mutate(group = factor(condition_Annotation, levels = custom_order)) %>%
  filter(!is.na(Stemness_UCell))

# 8. Prepare data frame for plotting
plot_df <- seu_obj@meta.data %>%
  mutate(group = factor(condition_Annotation, levels = custom_order)) %>%
  filter(!is.na(Exhaustion_UCell))

# 9. Prepare data frame for plotting
plot_df <- seu_obj@meta.data %>%
  mutate(group = factor(condition_Annotation, levels = custom_order)) %>%
  filter(!is.na(Effector_UCell))

# 10. Prepare data frame for plotting
plot_df <- seu_obj@meta.data %>%
  mutate(group = factor(condition_Annotation, levels = custom_order)) %>%
  filter(!is.na(Innate_UCell))


# 11. Plot violin plot of UCell score
p <- ggplot(plot_df, aes(x = group, y = Innate_UCell, fill = gsub("_.*", "", group))) +
  geom_violin(trim = TRUE, scale = "width", color = "black") +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = c("WT" = "white", "KO" = "salmon")) +
  labs(x = NULL, y = "Effector UCell Signature Score") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.y = element_text(face = "bold", size = 18),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

# 12. Display plot
p

# Save the plot
#ggsave("CD8_Ovarian_Stemness_UCellSignatureScore_violin_plot_Version2_20250721.svg", plot = p, width = 10, height = 6, dpi = 300)




