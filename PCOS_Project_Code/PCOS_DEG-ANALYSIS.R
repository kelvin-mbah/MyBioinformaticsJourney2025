# ================================================================
# Differential Expression Analysis using Limma (GSE98421 - PCOS vs Control)
# ================================================================

# ----- Load Required Libraries -----
# These packages handle statistical modeling, data manipulation, and visualization
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)


# ----- Define Sample Groups -----
# The variable 'groups' contain our sample labels (Control and PCOS)
table(groups)       # Check how many samples are in each group
levels(groups)      # Confirm the factor levels
groups <- factor(groups, levels = c("Control", "PCOS"))   # Ensure the right order


# ----- Create Design Matrix -----
# The design matrix represents our experimental setup
# Here we model the two conditions (Control and PCOS) without an intercept (~0)
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)    # Name columns as 'Control' and 'PCOS'
design                                # View the matrix to confirm


# ----- Fit Linear Model -----
# Fit a linear model for each gene to estimate expression differences
fit <- lmFit(data, design)

# Define the contrast: comparing PCOS vs Control
contrast_matrix <- makeContrasts(PCOS_vs_Control = PCOS - Control, levels = design)

# Apply the contrast to the fitted model
fit2 <- contrasts.fit(fit, contrast_matrix)


# ----- Apply Empirical Bayes Moderation -----
# This stabilizes variance across genes and improves statistical reliability
fit2 <- eBayes(fit2)


# ----- Extract Differentially Expressed Genes (DEGs) -----
# Get all genes with statistics (logFC, p-value, adj.P.Val, etc.)
deg_results <- topTable(fit2,
                        coef = "PCOS_vs_Control",
                        number = Inf,
                        adjust.method = "BH")   # BH = Benjamini–Hochberg FDR correction


# ----- Apply Thresholds -----
# Define cutoffs for fold-change and FDR significance
logFC_cutoff <- 0.7    # Genes changing more than ±0.7 log2 fold are considered strong
fdr_cutoff <- 0.05     # FDR < 0.05 = statistically significant

# Categorize genes based on their significance and direction
deg_results <- deg_results %>%
  mutate(Regulation = case_when(
    adj.P.Val < fdr_cutoff & logFC > logFC_cutoff ~ "Upregulated",
    adj.P.Val < fdr_cutoff & logFC < -logFC_cutoff ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# ----- Summary -----
# Count how many genes fall into each regulation category
table(deg_results$Regulation)


# ----- Save DEG Tables -----
# Save both full results and significant DEGs to CSV files
write.csv(deg_results, "DEG_FullResults_GSE98421.csv", row.names = TRUE)
write.csv(filter(deg_results, Regulation != "Not Significant"),
          "DEG_Significant_GSE98421.csv", row.names = TRUE)


# ----- Volcano Plot -----
# Visualize fold-change vs. significance of each gene
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: PCOS vs Control (GSE98421)",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value",
       color = "Regulation") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linetype = "dashed", color = "black")


# ----- Heatmap of Top DEGs -----
# Select the top 30 most significant DEGs for visualization
top_genes <- rownames(deg_results[deg_results$Regulation != "Not Significant", ])[1:30]
heatmap_data <- data[top_genes, ]

# Annotate columns with sample group info
annotation_col <- data.frame(Group = groups)
rownames(annotation_col) <- colnames(data)

# Plot the heatmap
pheatmap(heatmap_data,
         scale = "row",     # Normalize expression per gene
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 30 Differentially Expressed Genes",
         fontsize_row = 6,
         fontsize_col = 9)
