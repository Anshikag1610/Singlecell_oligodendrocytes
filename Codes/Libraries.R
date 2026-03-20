library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(clustree)
library(SingleR)
library(celldex)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(EnhancedVolcano)
library(DOSE)
library(ggalluvial)
library(viridis)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(cowplot)
library(scran)

install.packages(" ")

install.packages(c("remotes", "devtools"))
BiocManager::install(c(
  "SingleCellExperiment",
  "SummarizedExperiment",
  "DelayedArray",
  "batchelor",
  "BiocParallel",
  "HDF5Array"
))
#type n for updaten
# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("slingshot")

# Load required packages
library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(RColorBrewer)

BiocManager::install("ComplexHeatmap")
BiocManager::install("circlize")

library(ComplexHeatmap)
library(circlize)