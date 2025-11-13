# ================================================================
# PCA, Clustering, and QC Visualization for GSE98421 (PCOS Study)
# ================================================================

# ----- Load Required Packages -----
# These packages handle data processing, plotting, and clustering
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(factoextra)

# ----- Define Sample Groups from Phenotype Data -----
# Extract sample group (PCOS vs Control) using the 'title' column in phenotype data
groups <- ifelse(grepl("PCOS", phenotype_data_a$title, ignore.case = TRUE),
                 "PCOS", "Control")

# Convert group names to factor type for grouping in plots
groups <- factor(groups, levels = c("Control", "PCOS"))

# View how many samples are in each group
table(groups)

# Match group names to sample IDs (columns in the expression matrix)
names(groups) <- colnames(data)

# ================================================================
# Principal Component Analysis (PCA)
# ================================================================

# PCA helps visualize sample clustering and detect overall variation patterns
pca_res <- prcomp(t(data), scale. = TRUE)

# Extract first two principal components for plotting
pca_df <- data.frame(pca_res$x[, 1:2], Group = groups)
rownames(pca_df) <- colnames(data)

# ----- PCA Visualization -----
# Scatter plot showing clustering of PCOS vs Control samples
ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = rownames(pca_df)), size = 3, max.overlaps = Inf) +
  theme_minimal() +
  labs(
    title = "PCA of PCOS vs Control (GSE98421)",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")
  ) +
  scale_color_manual(values = c("Control" = "#1F77B4", "PCOS" = "#D62728")) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# ================================================================
# Scree Plot: Variance Explained by Principal Components
# ================================================================
fviz_eig(pca_res, addlabels = TRUE,
         barfill = "#3182bd", barcolor = "black") +
  ggtitle("Variance Explained by Principal Components")

# ================================================================
# Outlier Detection
# ================================================================
# Samples far from the cluster center are potential outliers
distances <- sqrt(pca_res$x[, 1]^2 + pca_res$x[, 2]^2)
outliers <- names(distances[distances > mean(distances) + 2 * sd(distances)])
cat("Potential outliers:", outliers, "\n")

# ================================================================
# Heatmap of Top Variable Genes
# ================================================================
# Select top 30 genes with highest expression variance
top_genes <- names(sort(apply(data, 1, var), decreasing = TRUE))[1:30]
heatmap_data <- data[top_genes, ]

# Add group color annotation (Control vs PCOS)
annotation_col <- data.frame(Group = groups)
rownames(annotation_col) <- colnames(heatmap_data)

# Plot heatmap
pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 30 Most Variable Genes (PCOS vs Control)",
  fontsize_row = 5,
  fontsize_col = 8,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 45,
  treeheight_row = 35,
  annotation_colors = list(Group = c(Control = "#1F77B4", PCOS = "#D62728"))
)

# ================================================================
# Hierarchical Clustering of Samples
# ================================================================
# Measures sample similarity and clusters based on expression profiles
dist_matrix <- dist(t(data))
hc <- hclust(dist_matrix, method = "complete")

plot(hc, labels = groups,
     main = "Hierarchical Clustering of Samples",
     sub = "", xlab = "")
