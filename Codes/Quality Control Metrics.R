# ============================================================
# LOAD REQUIRED LIBRARIES
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)   # REQUIRED for plot_annotation

# ============================================================
# LOAD SEURAT OBJECTS
# ============================================================

GSE115011 <- readRDS("GSE115011_seurat.rds")
GSE146373 <- readRDS("GSE146373_seurat.rds")

# ============================================================
# CALCULATE percent.mt (CRITICAL FIX)
# ============================================================

GSE115011[["percent.mt"]] <- PercentageFeatureSet(
  GSE115011,
  pattern = "^[Mm][Tt]-"
)

GSE146373[["percent.mt"]] <- PercentageFeatureSet(
  GSE146373,
  pattern = "^[Mm][Tt]-"
)

# ============================================================
# ADD DATASET LABELS
# ============================================================

GSE115011$dataset <- "GSE115011"
GSE146373$dataset <- "GSE146373"

# ============================================================
# MERGE OBJECTS
# ============================================================

combined <- merge(
  x = GSE115011,
  y = GSE146373,
  add.cell.ids = c("GSE115011", "GSE146373")
)

# Set identity for plotting
Idents(combined) <- combined$dataset

# ============================================================
# CREATE MERGED VIOLIN PLOT
# ============================================================

p <- VlnPlot(
  combined,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  cols = c("GSE115011" = "#008080", "GSE146373" = "#FF7F50"),
  pt.size = 0,
  ncol = 3
) +
  plot_annotation(
    title = "Quality Control Metrics "
  ) &
  theme_classic() &
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# Display plot
print(p)

# ============================================================
# SAVE FIGURE (PUBLICATION QUALITY)
# ============================================================

ggsave(
  filename = "Fig1B_QC_metrics.pdf",
  plot = p,
  width = 12,
  height = 5,
  dpi = 300,
  device = cairo_pdf
)
