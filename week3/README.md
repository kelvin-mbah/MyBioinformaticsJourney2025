# Week 3: Downstream Analysis with R

##  Overview
This week focused on downstream analysis of microarray data using R and Bioconductor packages. The analysis includes quality control, data filtering, and exploratory analysis of gene expression data.

##  Analysis Pipeline

### 1. Data Loading
- Loaded GEO series matrix file using GEOquery
- Extracted expression data, phenotypic data, and feature annotations

### 2. Quality Control
- Principal Component Analysis (PCA), boxplot for outlier detection
- Expression distribution visualization

### 3. Data Filtering
- *Intensity-based filtering*: Removed lowly expressed transcripts (threshold = 3.5)

### 4. Key Insights
- Discovered biological outliers representing real sample heterogeneity
- Balanced outlier distribution (3 normal vs 2 diseased samples)
- Decision to retain outliers for machine learning to capture biological spectrum

##  Tools & Packages Used

### Bioconductor Packages
- GEOquery - GEO data import
- oligo - Microarray data processing
- limma - Differential expression analysis
- affy - Alternative array processing

### Key Findings
- Successful identification of sample outliers
- Effective data cleaning while preserving biological signal
- Dataset prepared for machine learning applications

##  Next Steps
- Differential expression analysis
- Machine learning model development

##  Files
-  Complete R analysis script
- images/ - Analysis visualizations
- data/ - Processed data outputs

##  Lessons Learned
- Importance of biological context in outlier interpretation
- Balance between data cleaning and information preservation
- R/Bioconductor ecosystem power for genomic analysis
