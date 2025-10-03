#Topics Covered in this script -----
# 1. Quality Control (QC) 
# 2. RMA Normalization
# 3. Pre-processing and Filtering

#INSTALLATION OF REQUIRED PACKAGES -----
#Bioconductor repository contains R packages that are used for Omics-Analyis (Genomics, Transcriptomics, and proteomics)
#Install and load Biocmanager
if(!requireNamespace("BiocManager", quietly = TRUE))
install.packages('Biocmanager')
library(BiocManager)

#Install Required Bioconductor packages required to analyze Array dataset
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))
# Load Required Libraries
library(GEOquery)             # Download GEO datasets (series matrix, raw CEL files)
library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(dplyr)                # Data manipulation
library(limma)
library(AnnotationDbi)        #Interface for Annotation database

#DOWNLOAD OR LOAD SERIES MATRIX DATA INTO R ----

gse <- getGEO(filename = 'GSE146996_series_matrix.txt')

#Extract the Expression data
expression_data <- exprs(gse)

#Extract Feature data
feature_data <- fData(gse)

#Extract Phenotype data
phenotype_data <- pData(gse)

sum(is.na(phenotype_data$source_name_ch1))

#Read Cel Files into R
BiocManager::install('oligo')        #Affy could not readin my data type, it asked me to use oligo
library(oligo)                       #Load oligo library
cel_files <- list.celfiles("data/", full.names = TRUE)      
raw_data <- read.celfiles(cel_files)
raw_data


#QUALITY CONTROL BEFORE PRE-PROCESSING -----
#QC identifies outlier arrays, hybridization problems, or technical biases
#It uses different method to examine the data : boxplot, PCA, MA-plot, and Heatmap

arrayQualityMetrics(expressionset = raw_data,
                    outdir = 'Pre-Processing & Normalization/qc_raw_data',
                    force = TRUE,
                    do.logtransform = TRUE)


#RMA(Robust Multi-array Average) NORMALIZATION -----
#A well known method used to study affymetrix microarray data by:
# 1. Background correcting, 
# 2. normalizing probe intensities using quantile normalization and 
# 3. summarizing them into gene-level expression values using a robust median polish algorithm.

# This method reduces experimental variation across multiple arrays, 
# producing more symmetrical and reliable normalized expression data 
# compared to other approaches

normalized_data <- rma(raw_data)

#QC after Data Normalization
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = 'Pre-processing & Normalization/qc_normalized_data',
                    force = TRUE)

#Extract normalized expression 
processed_data <- as.data.frame(exprs(normalized_data))
dim(processed_data)


#FILTER LOW VARIANCE TRANSCRIPT ------
#Filtering removes probes with low signals

#Calculate the Median Intensity of the Probes
row_median <- rowMedians(as.matrix(processed_data))
row_median

#Plot the distribution of median intensities  of probes
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = 'Median Intensity Distribution')

threshold <- 3.5 

abline(v = threshold, col = 'black', lwd = 2)

#Select Probes above threshold
indx <- row_median > threshold
filtered_data <- processed_data[indx, ]

# Rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

#Override processsed_data with filtered data 
processed_data <- filtered_data

#PHENOTYPE DATA PREPARATION----
class(phenotype_data$source_name_ch1)

groups <- factor(phenotype_data$source_name_ch1,
                 levels = c('Adjacent gastric normal tissue', 'Gastric tumor tissue'),
                 labels = c('normal', 'cancer'))

class(groups)
levels(groups)
