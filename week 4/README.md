# Week 4: HTA 2.0 Microarray Differential Gene Expression Analysis  

## Overview  
This folder contains my R script (`HTA2.0_analysis.R`) for analyzing data from the Affymetrix Human Transcriptome Array 2.0 (HTA 2.0). 
The goal was to identify genes that differ in expression between cancer and normal samples using the `limma` package. 
The script maps probe IDs to gene symbols, handles genes with multiple probes, and generates plots like the volcano plot and heatmap.  

## Workflow Summary  
- **Probe Mapping:** Probe IDs were linked to gene symbols using `hta20transcriptcluster.db`, and unannotated probes were removed.  
- **Multiple Probes:** Averaged expression values for genes with multiple probes using `limma::avereps`.  
- **Differential Expression:** Performed analysis with `limma` comparing cancer vs normal samples. Found 158 upregulated and 87 downregulated genes (adjusted p < 0.05, |logFC| > 1).  
- **Visualization:**  
  - **Volcano Plot:** `Week_4/Gene Expression Analysis/volcano_plot.png` shows significant upregulated (red) and downregulated (blue) genes.  
  - **Heatmap:** `Week_4/heatmap_top25_DEGs.png` displays the top 25 most significantly changed genes.  

## Result Summary  
The analysis revealed clear differences in gene expression between cancer and normal samples. 
After handling multiple probes per gene, I identified key genes that were significantly up- or downregulated, reflecting major transcriptional changes linked to cancer.  
