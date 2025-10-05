# Bioconductor provides R packages for analyzing omics data (genomics, transcriptomics, proteomics etc).

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))

# Install CRAN packages for data manipulation
install.packages("dplyr")

# Load Required Libraries
library(GEOquery)             # Download GEO datasets (series matrix, raw CEL files)
library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics) # QC reports for microarray data
library(dplyr)                # Data manipulation

# -------------------------------------
#### Download Series Matrix Files ####
# -------------------------------------

# Series matrix files are preprocessed text files containing 
# expression values, sample annotations, and probe information.
# Reason: Useful for quick exploratory analysis when raw CEL files are not needed.
gse <- getGEO("GSE15852", GSEMatrix = TRUE)
gse <- gse[[1]]
exprSet <- exprs(gse)
pheno <- pData(gse)
featureData<-fData(gse)
sum(is.na(pheno$source_name_ch1))

library(limma)

# Normalize and filter for practice
exprSet_norm <- normalizeBetweenArrays(exprSet, method = "quantile")

# Filter low intensity
keep <- rowMeans(exprSet_norm) > 5
exprSet_filtered <- exprSet_norm[keep, ]

cat("Before filtering:", nrow(exprSet), " | After filtering:", nrow(exprSet_filtered), "\n")

# Relabel groups
pheno$Condition <- ifelse(grepl("normal", pheno$source_name_ch1, ignore.case = TRUE),
                          "Normal", "Cancer")

table(pheno$Condition)

library(arrayQualityMetrics)
getGEOSuppFiles("GSE15852", baseDir = "Raw_Data", makeDirectory = TRUE)

untar("Raw_Data/GSE15852/GSE15852_RAW.tar", exdir = "Raw_Data/CEL_Files")

row_median <- rowMedians(as.matrix(processed_data))

hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")