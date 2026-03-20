# ============================================================================
# COMPREHENSIVE SINGLE-CELL RNA-SEQ INTEGRATION & ANALYSIS PIPELINE
# ============================================================================

# ----------------------------
# 📦 LOAD ALL REQUIRED PACKAGES
# ----------------------------
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

# Create output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)
dir.create("outputs/qc", showWarnings = FALSE)

# ----------------------------
# 📊 1. LOAD DATASETS
# ----------------------------
cat("\n=== LOADING DATASETS ===\n")
GSE115011 <- readRDS("GSE115011_seurat.rds")
GSE146373 <- readRDS("GSE146373_seurat.rds")

cat("GSE115011 cells:", ncol(GSE115011), "\n")
cat("GSE146373 cells:", ncol(GSE146373), "\n")

# ----------------------------
# 🔍 2. QUALITY CONTROL ANALYSIS
# ----------------------------
cat("\n=== QUALITY CONTROL ===\n")

# Calculate mitochondrial percentage if not present
GSE115011[["percent.mt"]] <- PercentageFeatureSet(GSE115011, pattern = "^MT-")
GSE146373[["percent.mt"]] <- PercentageFeatureSet(GSE146373, pattern = "^MT-")

# QC metrics before filtering
qc_before <- data.frame(
  Dataset = c("GSE115011", "GSE146373"),
  N_Cells = c(ncol(GSE115011), ncol(GSE146373)),
  Median_Genes = c(median(GSE115011$nFeature_RNA), median(GSE146373$nFeature_RNA)),
  Median_UMIs = c(median(GSE115011$nCount_RNA), median(GSE146373$nCount_RNA)),
  Median_MT_pct = c(median(GSE115011$percent.mt), median(GSE146373$percent.mt))
)

write.csv(qc_before, "outputs/tables/QC_metrics_before_filtering.csv", row.names = FALSE)

# QC Violin plots
p_qc1 <- VlnPlot(GSE115011, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3, pt.size = 0) + 
  plot_annotation(title = "GSE115011 - QC Metrics")

p_qc2 <- VlnPlot(GSE146373, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3, pt.size = 0) + 
  plot_annotation(title = "GSE146373 - QC Metrics")

pdf("outputs/qc/QC_violin_plots.pdf", width = 12, height = 8)
print(p_qc1)
print(p_qc2)
dev.off()

# Feature scatter plots
p_scatter1 <- FeatureScatter(GSE115011, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("GSE115011")
p_scatter2 <- FeatureScatter(GSE146373, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("GSE146373")

pdf("outputs/qc/QC_scatter_plots.pdf", width = 12, height = 5)
print(p_scatter1 | p_scatter2)
dev.off()

# ----------------------------
# ⚙️ 3. PREPROCESSING & NORMALIZATION
# ----------------------------
cat("\n=== PREPROCESSING INDIVIDUAL DATASETS ===\n")

process_dataset <- function(obj, dataset_name) {
  cat("Processing", dataset_name, "...\n")
  
  DefaultAssay(obj) <- "RNA"
  
  # Normalize
  if (!"data" %in% Layers(obj)) {
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  
  # Find variable features
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  obj <- ScaleData(obj, features = rownames(obj))
  
  # PCA
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 50, verbose = FALSE)
  
  # Clustering
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5)
  
  # UMAP
  obj <- RunUMAP(obj, dims = 1:30)
  
  # Add dataset identifier
  obj$dataset <- dataset_name
  
  return(obj)
}

GSE115011 <- process_dataset(GSE115011, "GSE115011")
GSE146373 <- process_dataset(GSE146373, "GSE146373")

# Individual dataset UMAPs
p_umap_individual <- DimPlot(GSE115011, reduction = "umap", label = TRUE) + 
  ggtitle("GSE115011 - Individual Clustering") |
  DimPlot(GSE146373, reduction = "umap", label = TRUE) + 
  ggtitle("GSE146373 - Individual Clustering")

ggsave("outputs/figures/Fig1_individual_UMAPs.pdf", p_umap_individual, width = 14, height = 6)

# ----------------------------
# 🔗 4. INTEGRATION
# ----------------------------
cat("\n=== INTEGRATING DATASETS ===\n")

# Select integration features
features <- SelectIntegrationFeatures(object.list = list(GSE115011, GSE146373), 
                                      nfeatures = 3000)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = list(GSE115011, GSE146373), 
                                  anchor.features = features,
                                  dims = 1:30)

# Integrate data
combined <- IntegrateData(anchorset = anchors, dims = 1:30)

# ----------------------------
# 📈 5. INTEGRATED ANALYSIS
# ----------------------------
cat("\n=== ANALYZING INTEGRATED DATA ===\n")

DefaultAssay(combined) <- "integrated"

# Scale integrated data
combined <- ScaleData(combined, verbose = FALSE)

# Run PCA
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)

# Elbow plot for PC selection
pdf("outputs/qc/ElbowPlot_integrated.pdf", width = 8, height = 6)
ElbowPlot(combined, ndims = 50)
dev.off()

# UMAP
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

# Clustering at multiple resolutions
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
for (res in c(0.3, 0.5, 0.8, 1.0)) {
  combined <- FindClusters(combined, resolution = res)
}

# Clustree for resolution selection
pdf("outputs/figures/Fig2_clustree.pdf", width = 10, height = 12)
clustree(combined, prefix = "integrated_snn_res.")
dev.off()

# Set optimal resolution
Idents(combined) <- "integrated_snn_res.0.5"
combined$seurat_clusters <- Idents(combined)

# ----------------------------
# 🎨 6. MAIN VISUALIZATIONS
# ----------------------------
cat("\n=== GENERATING MAIN FIGURES ===\n")

# Figure 3: Integration overview
p1 <- DimPlot(combined, reduction = "umap", group.by = "dataset", pt.size = 0.5) + 
  ggtitle("By Dataset") +
  scale_color_manual(values = c("#E64B35", "#4DBBD5"))

p2 <- DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("By Cluster") + 
  NoLegend()

p3 <- DimPlot(combined, reduction = "umap", split.by = "dataset", pt.size = 0.5) + 
  ggtitle("Split by Dataset")

fig3 <- (p1 | p2) / p3
ggsave("outputs/figures/Fig3_integrated_UMAP.pdf", fig3, width = 16, height = 12)


# ----------------------------
# 🧬 7. CELL TYPE ANNOTATION
# ----------------------------
cat("\n=== CELL TYPE ANNOTATION ===\n")

# Automated annotation with SingleR
ref <- celldex::HumanPrimaryCellAtlasData()
#join RNA layers
combined <- JoinLayers(combined, assay = "RNA")
Layers(combined[["RNA"]])
DefaultAssay(combined) <- "RNA"
test_data <- combined[["RNA"]]$data

predictions <- SingleR(test = test_data, 
                       ref = ref, 
                       labels = ref$label.main,
                       de.method = "wilcox")

combined$celltype_auto <- predictions$labels
combined$celltype_scores <- predictions$scores

# Cell type UMAP
p_celltype <- DimPlot(combined, reduction = "umap", group.by = "celltype_auto", 
                      label = TRUE, repel = TRUE, pt.size = 0.5) +
  ggtitle("Automated Cell Type Annotation") +
  theme(legend.position = "bottom")

ggsave("outputs/figures/Fig4_celltype_annotation.pdf", p_celltype, width = 12, height = 10)

# Cell type composition by dataset
celltype_comp <- combined@meta.data %>%
  group_by(dataset, celltype_auto) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(dataset) %>%
  mutate(proportion = n / sum(n))

p_comp <- ggplot(celltype_comp, aes(x = celltype_auto, y = proportion, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type", y = "Proportion", fill = "Dataset") +
  scale_fill_manual(values = c("#E64B35", "#4DBBD5"))

ggsave("outputs/figures/Fig5_celltype_composition.pdf", p_comp, width = 12, height = 6)

# Statistical test for composition
celltype_table <- table(combined$dataset, combined$celltype_auto)
chisq_result <- chisq.test(celltype_table)

write.csv(celltype_table, "outputs/tables/celltype_counts_by_dataset.csv")
write.csv(data.frame(
  Test = "Chi-square test",
  Statistic = chisq_result$statistic,
  P_value = chisq_result$p.value
), "outputs/tables/celltype_composition_stats.csv", row.names = FALSE)

# ----------------------------
# 🔬 8. MARKER GENE IDENTIFICATION
# ----------------------------
cat("\n=== FINDING MARKER GENES ===\n")

# Integrated markers
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "seurat_clusters"

integrated_markers <- FindAllMarkers(combined, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25,
                                     test.use = "wilcox")

write.csv(integrated_markers, "outputs/tables/TableS1_integrated_markers.csv", row.names = FALSE)

# Individual dataset markers
DefaultAssay(GSE115011) <- "RNA"
DefaultAssay(GSE146373) <- "RNA"

markers_115011 <- FindAllMarkers(GSE115011, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

markers_146373 <- FindAllMarkers(GSE146373, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

write.csv(markers_115011, "outputs/tables/TableS2_markers_GSE115011.csv", row.names = FALSE)
write.csv(markers_146373, "outputs/tables/TableS3_markers_GSE146373.csv", row.names = FALSE)

# Top markers per cluster
top_integrated <- integrated_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

top115011 <- markers_115011 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

top146373 <- markers_146373 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

# ----------------------------
# 🧪 9. COMMON VS UNIQUE GENES ANALYSIS
# ----------------------------
cat("\n=== ANALYZING COMMON AND UNIQUE GENES ===\n")

common_genes <- intersect(unique(top115011$gene), unique(top146373$gene))
unique_115011 <- setdiff(unique(top115011$gene), unique(top146373$gene))
unique_146373 <- setdiff(unique(top146373$gene), unique(top115011$gene))

cat("✅ Number of common genes:", length(common_genes), "\n")
cat("✅ Unique to GSE115011:", length(unique_115011), "\n")
cat("✅ Unique to GSE146373:", length(unique_146373), "\n")

# Save gene lists
write.csv(data.frame(gene = common_genes), "outputs/tables/common_marker_genes.csv", row.names = FALSE)
write.csv(data.frame(gene = unique_115011), "outputs/tables/unique_genes_GSE115011.csv", row.names = FALSE)
write.csv(data.frame(gene = unique_146373), "outputs/tables/unique_genes_GSE146373.csv", row.names = FALSE)

# Venn diagram data
venn_data <- data.frame(
  Category = c("Common", "GSE115011 Unique", "GSE146373 Unique"),
  Count = c(length(common_genes), length(unique_115011), length(unique_146373))
)
write.csv(venn_data, "outputs/tables/gene_overlap_summary.csv", row.names = FALSE)

# ----------------------------
# 📊 10. DOT PLOTS
# ----------------------------
cat("\n=== GENERATING DOT PLOTS ===\n")

# Common genes dotplot
if (length(common_genes) > 0) {
  genes_to_plot <- head(common_genes, min(20, length(common_genes)))
  
  p_common <- DotPlot(combined, features = genes_to_plot) + 
    RotatedAxis() +
    ggtitle("Top Common Marker Genes") +
    theme(axis.text.x = element_text(size = 8))
  
  ggsave("outputs/figures/Fig6_dotplot_common_genes.pdf", p_common, width = 14, height = 8)
}

# Unique genes - GSE115011
if (length(unique_115011) > 0) {
  genes_to_plot <- head(unique_115011, min(20, length(unique_115011)))
  genes_to_plot <- intersect(genes_to_plot, rownames(combined))
  
  if (length(genes_to_plot) > 0) {
    p_unique1 <- DotPlot(combined, features = genes_to_plot) + 
      RotatedAxis() +
      ggtitle("Top Unique Genes: GSE115011") +
      theme(axis.text.x = element_text(size = 8))
    
    ggsave("outputs/figures/Fig7_dotplot_unique_GSE115011.pdf", p_unique1, width = 14, height = 8)
  }
}

# Unique genes - GSE146373
if (length(unique_146373) > 0) {
  genes_to_plot <- head(unique_146373, min(20, length(unique_146373)))
  genes_to_plot <- intersect(genes_to_plot, rownames(combined))
  
  if (length(genes_to_plot) > 0) {
    p_unique2 <- DotPlot(combined, features = genes_to_plot) + 
      RotatedAxis() +
      ggtitle("Top Unique Genes: GSE146373") +
      theme(axis.text.x = element_text(size = 8))
    
    ggsave("outputs/figures/Fig8_dotplot_unique_GSE146373.pdf", p_unique2, width = 14, height = 8)
  }
}

# Combined dotplot
all_genes_plot <- c(
  head(common_genes, min(10, length(common_genes))),
  head(intersect(unique_115011, rownames(combined)), 5),
  head(intersect(unique_146373, rownames(combined)), 5)
)

if (length(all_genes_plot) > 0) {
  p_combined_dot <- DotPlot(combined, features = all_genes_plot) + 
    RotatedAxis() +
    ggtitle("Combined: Common + Unique Genes") +
    theme(axis.text.x = element_text(size = 8))
  
  ggsave("outputs/figures/Fig9_dotplot_combined.pdf", p_combined_dot, width = 14, height = 8)
}

# ----------------------------
# 🌡️ 11. HEATMAPS
# ----------------------------
cat("\n=== GENERATING HEATMAPS ===\n")

DefaultAssay(combined) <- "integrated"
scaled_data <- GetAssayData(combined, layer = "scale.data")

# Order cells by cluster
cell_order <- sort(combined$seurat_clusters)

# Heatmap: Common genes
if (length(common_genes) > 0) {
  common_genes_in_data <- intersect(head(common_genes, 20), rownames(scaled_data))
  
  if (length(common_genes_in_data) > 0) {
    common_expr <- scaled_data[common_genes_in_data, names(cell_order)]
    
    pdf("outputs/figures/Fig10_heatmap_common_genes.pdf", width = 12, height = 10)
    pheatmap(common_expr, 
             cluster_rows = TRUE, 
             cluster_cols = FALSE,
             show_colnames = FALSE,
             annotation_col = data.frame(
               Cluster = cell_order,
               Dataset = combined$dataset[names(cell_order)],
               row.names = names(cell_order)
             ),
             main = "Heatmap: Top Common Genes",
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
    dev.off()
  }
}

# Heatmap: Unique GSE115011
unique_115011_in_data <- intersect(head(unique_115011, 20), rownames(scaled_data))
if (length(unique_115011_in_data) > 0) {
  unique115011_expr <- scaled_data[unique_115011_in_data, names(cell_order)]
  
  pdf("outputs/figures/Fig11_heatmap_unique_GSE115011.pdf", width = 12, height = 10)
  pheatmap(unique115011_expr, 
           cluster_rows = TRUE, 
           cluster_cols = FALSE,
           show_colnames = FALSE,
           annotation_col = data.frame(
             Cluster = cell_order,
             Dataset = combined$dataset[names(cell_order)],
             row.names = names(cell_order)
           ),
           main = "Heatmap: Unique Genes - GSE115011",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  dev.off()
}

# Heatmap: Unique GSE146373
unique_146373_in_data <- intersect(head(unique_146373, 20), rownames(scaled_data))
if (length(unique_146373_in_data) > 0) {
  unique146373_expr <- scaled_data[unique_146373_in_data, names(cell_order)]
  
  pdf("outputs/figures/Fig12_heatmap_unique_GSE146373.pdf", width = 12, height = 10)
  pheatmap(unique146373_expr, 
           cluster_rows = TRUE, 
           cluster_cols = FALSE,
           show_colnames = FALSE,
           annotation_col = data.frame(
             Cluster = cell_order,
             Dataset = combined$dataset[names(cell_order)],
             row.names = names(cell_order)
           ),
           main = "Heatmap: Unique Genes - GSE146373",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  dev.off()
}

# ----------------------------
# 🔥 12. VOLCANO PLOTS - DATASET COMPARISON
# ----------------------------
cat("\n=== DIFFERENTIAL EXPRESSION BETWEEN DATASETS ===\n")

DefaultAssay(combined) <- "RNA"
Idents(combined) <- "dataset"

dataset_de <- FindMarkers(combined, 
                          ident.1 = "GSE115011", 
                          ident.2 = "GSE146373",
                          logfc.threshold = 0.25,
                          min.pct = 0.1)

write.csv(dataset_de, "outputs/tables/TableS4_dataset_differential_expression.csv")

# Volcano plot
dataset_de$gene <- rownames(dataset_de)

pdf("outputs/figures/Fig13_volcano_dataset_comparison.pdf", width = 12, height = 10)
EnhancedVolcano(dataset_de,
                lab = dataset_de$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'GSE115011 vs GSE146373',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 0.5,
                drawConnectors = TRUE,
                widthConnectors = 0.5)
dev.off()

# ----------------------------
# 🎯 13. PATHWAY ENRICHMENT ANALYSIS
# ----------------------------
cat("\n=== PATHWAY ENRICHMENT ANALYSIS ===\n")

# GO Enrichment for common genes
if (length(common_genes) >= 10) {
  ego_common <- enrichGO(gene = common_genes,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'SYMBOL',
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)
  
  if (nrow(ego_common@result) > 0) {
    write.csv(ego_common@result, "outputs/tables/TableS5_GO_enrichment_common_genes.csv")
    
    pdf("outputs/figures/Fig14_GO_enrichment_common.pdf", width = 12, height = 10)
    print(dotplot(ego_common, showCategory = 20) + ggtitle("GO Enrichment: Common Genes"))
    print(barplot(ego_common, showCategory = 20) + ggtitle("GO Enrichment: Common Genes"))
    dev.off()
  }
}

# GO Enrichment for unique GSE115011 genes
if (length(unique_115011) >= 10) {
  ego_unique1 <- enrichGO(gene = unique_115011,
                          OrgDb = org.Hs.eg.db,
                          keyType = 'SYMBOL',
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
  
  if (nrow(ego_unique1@result) > 0) {
    write.csv(ego_unique1@result, "outputs/tables/TableS6_GO_enrichment_unique_GSE115011.csv")
    
    pdf("outputs/figures/Fig15_GO_enrichment_unique_GSE115011.pdf", width = 12, height = 10)
    print(dotplot(ego_unique1, showCategory = 20) + ggtitle("GO Enrichment: Unique GSE115011"))
    dev.off()
  }
}

# GO Enrichment for unique GSE146373 genes
if (length(unique_146373) >= 10) {
  ego_unique2 <- enrichGO(gene = unique_146373,
                          OrgDb = org.Hs.eg.db,
                          keyType = 'SYMBOL',
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
  
  if (nrow(ego_unique2@result) > 0) {
    write.csv(ego_unique2@result, "outputs/tables/TableS7_GO_enrichment_unique_GSE146373.csv")
    
    pdf("outputs/figures/Fig16_GO_enrichment_unique_GSE146373.pdf", width = 12, height = 10)
    print(dotplot(ego_unique2, showCategory = 20) + ggtitle("GO Enrichment: Unique GSE146373"))
    dev.off()
  }
}

# KEGG Pathway Analysis for common genes
if (length(common_genes) >= 10) {
  # Convert to Entrez IDs
  gene_ids <- bitr(common_genes, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) >= 10) {
    kegg_common <- enrichKEGG(gene = gene_ids$ENTREZID,
                              organism = 'hsa',
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)
    
    if (nrow(kegg_common@result) > 0) {
      write.csv(kegg_common@result, "outputs/tables/TableS8_KEGG_enrichment_common_genes.csv")
      
      pdf("outputs/figures/Fig17_KEGG_enrichment_common.pdf", width = 12, height = 10)
      print(dotplot(kegg_common, showCategory = 20) + ggtitle("KEGG Pathways: Common Genes"))
      dev.off()
    }
  }
}

# ----------------------------
# 🔬 14. CELL CYCLE ANALYSIS
# ----------------------------
cat("\n=== CELL CYCLE ANALYSIS ===\n")

# Score cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

combined <- CellCycleScoring(combined, 
                             s.features = s.genes, 
                             g2m.features = g2m.genes, 
                             set.ident = TRUE)

# Cell cycle UMAP
p_cc1 <- DimPlot(combined, reduction = "umap", group.by = "Phase", pt.size = 0.5) +
  ggtitle("Cell Cycle Phase")

p_cc2 <- DimPlot(combined, reduction = "umap", group.by = "Phase", 
                 split.by = "dataset", pt.size = 0.5) +
  ggtitle("Cell Cycle by Dataset")

pdf("outputs/figures/Fig18_cell_cycle.pdf", width = 16, height = 6)
print(p_cc1 | p_cc2)
dev.off()

# Cell cycle distribution
cc_table <- table(combined$Phase, combined$dataset)
cc_prop <- prop.table(cc_table, margin = 2)

write.csv(cc_table, "outputs/tables/cell_cycle_counts.csv")
write.csv(cc_prop, "outputs/tables/cell_cycle_proportions.csv")

# Chi-square test
cc_chisq <- chisq.test(cc_table)
write.csv(data.frame(
  Test = "Cell Cycle Chi-square",
  Statistic = cc_chisq$statistic,
  P_value = cc_chisq$p.value
), "outputs/tables/cell_cycle_stats.csv", row.names = FALSE)

# ----------------------------
# 📈 15. FEATURE PLOTS FOR KEY MARKERS
# ----------------------------
cat("\n=== GENERATING FEATURE PLOTS ===\n")

# Top common genes feature plots
if (length(common_genes) >= 9) {
  top9_common <- head(common_genes, 9)
  
  pdf("outputs/figures/Fig19_featureplots_common_genes.pdf", width = 15, height = 15)
  print(FeaturePlot(combined, features = top9_common, ncol = 3, pt.size = 0.5))
  dev.off()
}

# Top unique genes - split by dataset
if (length(unique_115011) >= 4 && length(unique_146373) >= 4) {
  top4_each <- c(head(intersect(unique_115011, rownames(combined)), 4),
                 head(intersect(unique_146373, rownames(combined)), 4))
  
  if (length(top4_each) > 0) {
    pdf("outputs/figures/Fig20_featureplots_unique_genes.pdf", width = 15, height = 15)
    print(FeaturePlot(combined, features = top4_each, ncol = 4, pt.size = 0.5))
    dev.off()
  }
}

# ----------------------------
# 🎨 16. SPLIT VIOLIN PLOTS
# ----------------------------
cat("\n=== GENERATING VIOLIN PLOTS ===\n")

if (length(common_genes) >= 4) {
  top4_violin <- head(common_genes, 4)
  
  pdf("outputs/figures/Fig21_violin_common_genes.pdf", width = 14, height = 10)
  print(VlnPlot(combined, features = top4_violin, 
                split.by = "dataset", ncol = 2, pt.size = 0))
  dev.off()
}

# ----------------------------
# 🌊 17. ALLUVIAL DIAGRAM - CLUSTER CORRESPONDENCE
# ----------------------------
cat("\n=== CLUSTER CORRESPONDENCE ANALYSIS ===\n")

# Prepare data for alluvial plot
cluster_flow <- combined@meta.data %>%
  group_by(dataset, seurat_clusters) %>%
  summarise(count = n(), .groups = 'drop')

p_alluvial <- ggplot(cluster_flow,
                     aes(axis1 = dataset, axis2 = seurat_clusters, y = count)) +
  geom_alluvium(aes(fill = dataset), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/12, fill = "white", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Dataset", "Cluster"), expand = c(.05, .05)) +
  scale_fill_manual(values = c("#E64B35", "#4DBBD5")) +
  theme_minimal() +
  ggtitle("Cell Flow: Dataset to Cluster") +
  theme(legend.position = "bottom")

ggsave("outputs/figures/Fig22_alluvial_cluster_flow.pdf", p_alluvial, 
       width = 10, height = 8)

# ----------------------------
# 📊 18. CLUSTER SIMILARITY HEATMAP
# ----------------------------
cat("\n=== CLUSTER SIMILARITY ANALYSIS ===\n")

# Calculate Jaccard similarity between clusters across datasets
calculate_jaccard_similarity <- function(obj) {
  clusters <- sort(unique(obj$seurat_clusters))
  n_clusters <- length(clusters)
  
  similarity_matrix <- matrix(0, nrow = n_clusters, ncol = n_clusters,
                              dimnames = list(paste0("C", clusters), 
                                              paste0("C", clusters)))
  
  for (i in seq_along(clusters)) {
    for (j in seq_along(clusters)) {
      cells_i <- WhichCells(obj, expression = seurat_clusters == clusters[i])
      cells_j <- WhichCells(obj, expression = seurat_clusters == clusters[j])
      
      intersection <- length(intersect(cells_i, cells_j))
      union <- length(union(cells_i, cells_j))
      
      similarity_matrix[i, j] <- ifelse(union > 0, intersection / union, 0)
    }
  }
  
  return(similarity_matrix)
}

similarity_mat <- calculate_jaccard_similarity(combined)

pdf("outputs/figures/Fig23_cluster_similarity.pdf", width = 10, height = 10)
pheatmap(similarity_mat,
         main = "Cluster Similarity (Jaccard Index)",
         color = colorRampPalette(c("white", "red"))(100),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 8)
dev.off()

# ----------------------------
# 🧬 19. GENE SIGNATURE SCORING
# ----------------------------
cat("\n=== GENE SIGNATURE ANALYSIS ===\n")

# Define gene signatures (customize based on your biology)
stress_genes <- c("HSP90AA1", "HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "HSP90AB1")
hypoxia_genes <- c("HIF1A", "VEGFA", "LDHA", "SLC2A1", "PDK1", "BNIP3")
inflammation_genes <- c("IL6", "TNF", "IL1B", "CXCL8", "CCL2", "PTGS2")
proliferation_genes <- c("MKI67", "TOP2A", "PCNA", "MCM2", "CCNB1", "CCNA2")

# Score signatures
combined <- AddModuleScore(combined, 
                           features = list(stress_genes), 
                           name = "Stress_Score")

combined <- AddModuleScore(combined, 
                           features = list(hypoxia_genes), 
                           name = "Hypoxia_Score")

combined <- AddModuleScore(combined, 
                           features = list(inflammation_genes), 
                           name = "Inflammation_Score")

combined <- AddModuleScore(combined, 
                           features = list(proliferation_genes), 
                           name = "Proliferation_Score")

# Feature plots for signatures
pdf("outputs/figures/Fig24_gene_signatures.pdf", width = 16, height = 12)
print(FeaturePlot(combined, 
                  features = c("Stress_Score1", "Hypoxia_Score1", 
                               "Inflammation_Score1", "Proliferation_Score1"),
                  ncol = 2, pt.size = 0.5))
dev.off()

# Violin plots by dataset
pdf("outputs/figures/Fig25_gene_signatures_violin.pdf", width = 14, height = 10)
print(VlnPlot(combined, 
              features = c("Stress_Score1", "Hypoxia_Score1", 
                           "Inflammation_Score1", "Proliferation_Score1"),
              split.by = "dataset", ncol = 2, pt.size = 0))
dev.off()

# Statistical comparison of signatures between datasets
signature_stats <- data.frame()
for (sig in c("Stress_Score1", "Hypoxia_Score1", "Inflammation_Score1", "Proliferation_Score1")) {
  gse1_vals <- combined@meta.data[combined$dataset == "GSE115011", sig]
  gse2_vals <- combined@meta.data[combined$dataset == "GSE146373", sig]
  
  wilcox_test <- wilcox.test(gse1_vals, gse2_vals)
  
  signature_stats <- rbind(signature_stats, data.frame(
    Signature = sig,
    GSE115011_mean = mean(gse1_vals),
    GSE146373_mean = mean(gse2_vals),
    P_value = wilcox_test$p.value
  ))
}

write.csv(signature_stats, "outputs/tables/TableS9_signature_comparison.csv", row.names = FALSE)

# ----------------------------
# 📋 20. COMPREHENSIVE SUMMARY TABLES
# ----------------------------
cat("\n=== GENERATING SUMMARY TABLES ===\n")

# Table 1: Dataset summary
summary_table <- data.frame(
  Metric = c("Total Cells", "Total Genes Detected", "Median Genes/Cell", 
             "Median UMIs/Cell", "Number of Clusters", "Number of Cell Types"),
  GSE115011 = c(
    ncol(GSE115011),
    length(rownames(GSE115011)[rowSums(GetAssayData(GSE115011, layer = "counts")) > 0]),
    median(GSE115011$nFeature_RNA),
    median(GSE115011$nCount_RNA),
    length(unique(GSE115011$seurat_clusters)),
    length(unique(combined$celltype_auto[combined$dataset == "GSE115011"]))
  ),
  GSE146373 = c(
    ncol(GSE146373),
    length(rownames(GSE146373)[rowSums(GetAssayData(GSE146373, layer = "counts")) > 0]),
    median(GSE146373$nFeature_RNA),
    median(GSE146373$nCount_RNA),
    length(unique(GSE146373$seurat_clusters)),
    length(unique(combined$celltype_auto[combined$dataset == "GSE146373"]))
  ),
  Integrated = c(
    ncol(combined),
    length(rownames(combined)[rowSums(GetAssayData(combined, layer = "counts")) > 0]),
    median(combined$nFeature_RNA),
    median(combined$nCount_RNA),
    length(unique(combined$seurat_clusters)),
    length(unique(combined$celltype_auto))
  )
)

write.csv(summary_table, "outputs/tables/Table1_dataset_summary.csv", row.names = FALSE)

# Table 2: Cluster cell counts
cluster_summary <- combined@meta.data %>%
  group_by(seurat_clusters, dataset) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = dataset, values_from = n_cells, values_fill = 0) %>%
  mutate(Total = GSE115011 + GSE146373)

write.csv(cluster_summary, "outputs/tables/Table2_cluster_cell_counts.csv", row.names = FALSE)

# Table 3: Top markers per cluster (top 5)
top5_per_cluster <- integrated_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 5) %>%
  dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

write.csv(top5_per_cluster, "outputs/tables/Table3_top5_markers_per_cluster.csv", row.names = FALSE)

# ----------------------------
# 💾 21. SAVE FINAL OBJECTS
# ----------------------------
cat("\n=== SAVING FINAL OBJECTS ===\n")

saveRDS(combined, "outputs/combined_integrated_final.rds")
saveRDS(GSE115011, "outputs/GSE115011_processed.rds")
saveRDS(GSE146373, "outputs/GSE146373_processed.rds")

# ----------------------------
# 📄 22. SESSION INFO & REPRODUCIBILITY
# ----------------------------
cat("\n=== SAVING SESSION INFO ===\n")

sink("outputs/sessionInfo.txt")
print(sessionInfo())
sink()

# Create analysis log
analysis_log <- data.frame(
  Step = c("Data Loading", "QC", "Preprocessing", "Integration", "Clustering", 
           "Cell Type Annotation", "Marker Identification", "Pathway Enrichment",
           "Cell Cycle Analysis", "Signature Scoring", "Visualization"),
  Status = rep("Completed", 11),
  Timestamp = Sys.time()
)

write.csv(analysis_log, "outputs/analysis_log.csv", row.names = FALSE)

# ----------------------------
# 📊 23. FINAL MULTI-PANEL FIGURE
# ----------------------------
cat("\n=== GENERATING FINAL MULTI-PANEL FIGURE ===\n")

p_final_1 <- DimPlot(combined, reduction = "umap", group.by = "dataset", pt.size = 0.3) + 
  ggtitle("A) By Dataset") + NoLegend()

p_final_2 <- DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.3) + 
  ggtitle("B) By Cluster") + NoLegend()

p_final_3 <- DimPlot(combined, reduction = "umap", group.by = "celltype_auto", 
                     label = TRUE, repel = TRUE, pt.size = 0.3) + 
  ggtitle("C) Cell Types") + NoLegend()

p_final_4 <- DimPlot(combined, reduction = "umap", group.by = "Phase", pt.size = 0.3) + 
  ggtitle("D) Cell Cycle") + NoLegend()

# Combine into final figure
final_figure <- (p_final_1 | p_final_2) / (p_final_3 | p_final_4)

ggsave("outputs/figures/Figure_MAIN_comprehensive.pdf", final_figure, 
       width = 16, height = 16)

# ----------------------------
# ✅ ANALYSIS COMPLETE
# ----------------------------
cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("              ✅ ANALYSIS PIPELINE COMPLETED!                   \n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("\n📁 All outputs saved to 'outputs/' directory:\n")
cat("   - figures/     : All PDF figures\n")
cat("   - tables/      : All CSV tables\n")
cat("   - qc/          : Quality control plots\n")
cat("   - *.rds files  : Saved Seurat objects\n")
cat("\n📊 Summary Statistics:\n")
cat("   - Total cells (integrated):", ncol(combined), "\n")
cat("   - Number of clusters:", length(unique(combined$seurat_clusters)), "\n")
cat("   - Number of cell types:", length(unique(combined$celltype_auto)), "\n")
cat("   - Common marker genes:", length(common_genes), "\n")
cat("   - Unique to GSE115011:", length(unique_115011), "\n")
cat("   - Unique to GSE146373:", length(unique_146373), "\n")
cat("\n═══════════════════════════════════════════════════════════════\n")


