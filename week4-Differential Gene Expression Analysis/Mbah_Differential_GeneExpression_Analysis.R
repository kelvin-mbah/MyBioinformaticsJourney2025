# Install and Load Required Packages----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Load Bioconductor packages
library(AnnotationDbi)   # for annotation and probe–gene mapping
library(limma)           # Performs linear modeling and differential expression
library(dplyr)           # For data manipulation 
library(tibble)          
library(ggplot2)         # Used for creating plots and visualization
library(pheatmap)        # Generates heatmap for gene expression data

# Probe-to-Gene Mapping using AnnotationDbi----
#NB: Another database that can be used in place for AnnotationDbi includes;
     # 1. Ensembl Genome Browser 
     # 2. DAVID Functional Annotation Tool

# check annotation slot of the dataset being used

annotation(raw_data)
raw_data

# Install hta20transcriptcluster.db required for pd.hta.2.0
BiocManager::install("hta20transcriptcluster.db")
library(hta20transcriptcluster.db)
library(org.Hs.eg.db)

# Display objects available in the annotation package
ls("package:hta20transcriptcluster.db")
columns(hta20transcriptcluster.db)
keytypes(hta20transcriptcluster.db)

# Extract probe IDs from processed microarray data----
probe_ids <- rownames(processed_data)

# Map probe IDs to gene symbols using the platform annotation database
gene_symbols <- mapIds(
  hta20transcriptcluster.db,     #Database used for mapping
  keys = probe_ids,              #Input probe IDs
  keytype = 'PROBEID',           #PROBE ID keytype
  column = 'SYMBOL',             #Desired annotation column should be (gene symbols)
  multiVals = 'first'            #Return First match if multiple exist
  )

# Convert mapping to a data frame and rename columns
mapped_gene <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2)


# Handle multiple probes mapping to a single gene----
# Several strategies exist:
# 1. Retain probe with highest expression or variance
# 2. Average or summarize probe signals
# 3. Remove duplicate probes to maintain one row per gene


# Summarize number of probes per gene symbol
duplicate_summary <- mapped_gene %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

# Identify genes associated with multiple probes
duplicate_genes <- duplicate_summary %>%
  filter(probes_per_gene > 1)

sum(duplicate_genes$probes_per_gene)

# Merge annotation mapping with expression data----

# Verify if probe IDs in mapping correspond to expression data
all(mapped_gene$PROBEID == row.names(processed_data))

#Merge annotation (SYMBOL) with expression matrix
processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)

# Remove probes without valid gene symbol annotation
processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

# Select only numeric expression columns
expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

# Collapse multiple probes per gene using average expression----


# limma::avereps() computes the average for probes representing the same gene
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

# Example to demonstrate how avereps works
x <- matrix(rnorm(8*3), 8, 3)
colnames(x) <- c("S1", "S2", "S3")
rownames(x) <- c("b", "a", "a", "c", "c", "b", "b", "b")
head(x)
avereps(x)  # Collapses duplicated row names by averaging

dim(averaged_data)

# Convert averaged expression data to matrix format
data <- as.data.frame(averaged_data)
data <- data.matrix(data)
str(data)        # Structure check
is.numeric(data) # Confirm numeric matrix


# Differential Gene Expression Analysis----

# Define sample groups based on phenotype data
# Adjust group labels according to dataset annotation
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("Gastric tumor tissue", "Adjacent gastric normal tissue"),
                 label = c( "cancer", "normal"))

class(groups)
levels(groups)

# Create design matrix for linear modeling----

# Using no intercept (~0 + groups) allows each group to have its own coefficient
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Fit linear model to expression data
fit_1 <- lmFit(data, design)

# Define contrast to compare cancer vs normal samples
contrast_matrix <- makeContrasts(cancer_vs_normal = cancer - normal,
                                 levels = design)

# Apply contrasts and compute moderated statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)

fit_2 <- eBayes(fit_contrast)

# Extract list of differentially expressed genes (DEGs)----

deg_results <- topTable(fit_2,
                        coef = "cancer_vs_normal",  # Specify contrast of interest
                        number = Inf,               # Return all genes
                        adjust.method = "BH")       # Benjamini-Hochberg correction

# Classify DEGs into Upregulated, Downregulated, or Not Significant----

deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "No")
))

# Subset genes by regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")

# Combine both sets of DEGs
deg_updown <- rbind(upregulated, downregulated)

write.csv(deg_results, file = "DEGs_Results.csv")
write.csv(upregulated, file = "Upregulated_DEGs.csv")
write.csv(downregulated, file = "Downregulated_DEGs.csv")
write.csv(deg_updown, file = "Updown_DEGs.csv")

#### Data Visualization ----
# Volcano Plot: visualizes DEGs by logFC and adjusted p-values

# Note: x-axis = log2 fold change, y-axis = -log10 adjusted p-value

# Save volcano plot as PNG

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(P-value)",
       color = "Regulation")

dev.off()

# Heatmap of Top Differentially Expressed Genes----

# Select top genes with smallest adjusted p-values
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 25)

# Subset averaged expression matrix for selected genes
heatmap_data <- data[top_genes, ]

# Generate unique column names per sample group for display
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))

# Assign formatted names to heatmap columns
colnames(heatmap_data) <- heatmap_names

# Save heatmap as PNG
png("heatmap_top25_DEGs.png", width = 2000, height = 1500, res = 300)

# Generate heatmap without additional scaling
pheatmap(
  heatmap_data,
  scale = "none", # for already normalized data
  cluster_rows = FALSE,              # Disable row clustering
  cluster_cols = TRUE,               # Cluster samples
  show_rownames = TRUE,              # Display gene names
  show_colnames = TRUE,              # Display sample labels
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 25 Differentially Expressed Genes"
)

dev.off()


#RESULT SUMMARY----

# HTA 2.0’s design includes probes for different transcripts or isoforms of the same gene, 
# leading to multiple probe IDs (e.g., TC01000001.hg.1, TC01000002.hg.1) mapping to the same gene symbol (e.g., BRCA1). 
# This is evident from my duplicate_gene result that have 2634 duplicated genes with NA inclusive. 

# The 'cancer_vs_normal = cancer - normal' contrast was performed

# Overall. we have a total of Upregulated = 158 genes , Downregulated = 87 genes

