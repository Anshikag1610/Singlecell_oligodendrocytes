library(Seurat)
library(pheatmap)
library(RColorBrewer)

# Load data
combined <- readRDS("outputs/combined_integrated_final.rds")
de_results <- read.csv("TableS4_dataset_differential_expression.csv", row.names = 1)

# ============================================================
# STEP 1: Select top DE genes
# ============================================================

# Top 25 upregulated
de_up <- de_results[order(de_results$avg_log2FC, decreasing = TRUE), ]
top25_up <- head(rownames(de_up), 25)

# Top 25 downregulated
de_down <- de_results[order(de_results$avg_log2FC, decreasing = FALSE), ]
top25_down <- head(rownames(de_down), 25)

# Combine and ensure uniqueness
top50_genes <- unique(c(top25_up, top25_down))

# ============================================================
# STEP 2: Use RNA assay (IMPORTANT FIX)
# ============================================================

DefaultAssay(combined) <- "RNA"

# Scale data if not already scaled
combined <- ScaleData(combined, features = rownames(combined), verbose = FALSE)

expr_data <- GetAssayData(combined, layer = "scale.data")

# Keep only genes present in object
top50_genes <- top50_genes[top50_genes %in% rownames(expr_data)]

expr_subset <- expr_data[top50_genes, ]

# ============================================================
# STEP 3: Order cells by dataset
# ============================================================

cell_order <- order(combined$dataset)
expr_ordered <- expr_subset[, cell_order]

# ============================================================
# STEP 4: Create annotation
# ============================================================

# Fix cluster labels
annotation_col$Cluster <- paste0("C", as.character(annotation_col$Cluster))
annotation_col$Cluster <- factor(annotation_col$Cluster, levels = paste0("C", 0:5))

# Colors
ann_colors <- list(
  Dataset = c(
    "GSE115011" = "#008080",
    "GSE146373" = "#FF7F50"
  ),
  Cluster = setNames(
    RColorBrewer::brewer.pal(6, "Set2"),
    paste0("C", 0:5)
  )
)

# Plot
pdf("Fig3C_heatmap_top50_DE.pdf", width = 12, height = 10)

pheatmap(
  expr_ordered,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks = seq(-2.5, 2.5, length.out = 101),
  fontsize_row = 9,
  main = "Top 50 Differentially Expressed Genes Across Datasets",
  border_color = NA
)

dev.off()
