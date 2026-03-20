combined_integrated_final <- readRDS("outputs/combined_integrated_final.rds")
# ============================================================
# STEP 1: CONVERT SEURAT TO SINGLECELLEXPERIMENT
# ============================================================

# Slingshot works with SingleCellExperiment objects
library(SingleCellExperiment)

# Convert your integrated Seurat object
sce <- as.SingleCellExperiment(combined_integrated_final, assay = "integrated")

# Add UMAP coordinates to reducedDims
reducedDim(sce, "UMAP") <- Embeddings(combined_integrated_final, "umap")

# Verify conversion
print(sce)
print(dim(reducedDim(sce, "UMAP")))

# ============================================================
# STEP 2: RUN SLINGSHOT TRAJECTORY INFERENCE
# ============================================================

# Basic slingshot run
# - Uses UMAP coordinates for trajectory
# - Uses cluster assignments to guide lineage inference
# - Automatically infers lineage structure

sds <- slingshot(sce, 
                 clusterLabels = 'seurat_clusters',  # Your cluster column name
                 reducedDim = 'UMAP',
                 start.clus = "2")  # Start from proliferative cluster (TOP2A/KIFC1)
# Remove start.clus to let slingshot auto-detect

# Check results
print(sds)
summary(sds@metadata$slingshot)

# ============================================================
# STEP 3: EXTRACT PSEUDOTIME VALUES
# ============================================================

# Get pseudotime for all lineages
pseudotime_values <- slingPseudotime(sds)

# Check how many lineages were inferred
n_lineages <- ncol(pseudotime_values)
print(paste("Number of lineages detected:", n_lineages))

# Add pseudotime back to Seurat object for easy plotting
for(i in 1:n_lineages) {
  combined_integrated_final@meta.data[, paste0("slingPseudotime_", i)] <- pseudotime_values[, i]
  
}

# Add lineage assignment (which lineage each cell belongs to)
lineage_weights <- slingCurveWeights(sds)
combined_integrated_final$sling_lineage <- apply(lineage_weights, 1, which.max)

# ============================================================
# STEP 4: VISUALIZATION - TRAJECTORY ON UMAP
# ============================================================

# Plot 1: UMAP with trajectory curves overlaid
library(RColorBrewer)

# Extract curve data for plotting
curves <- slingCurves(sds)

# Base UMAP plot
p1 <- DimPlot(combined_integrated_final, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Slingshot Trajectory - Clusters")

# Add slingshot curves
for(i in seq_along(curves)) {
  curve_data <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
  colnames(curve_data) <- c("UMAP_1", "UMAP_2")
  
  p1 <- p1 + geom_path(data = curve_data, 
                       aes(x = UMAP_1, y = UMAP_2),
                       size = 2, 
                       color = "black",
                       arrow = arrow(length = unit(0.3, "cm"), type = "closed"))
}

print(p1)
ggsave("outputs/figures/slingshot_trajectory_clusters.pdf",
       p1,
       width = 10,
       height = 8)


# ============================================================
# STEP 5: PSEUDOTIME VISUALIZATION
# ============================================================

# Plot 2: UMAP colored by pseudotime (Lineage 1)
p2 <- FeaturePlot(combined_integrated_final, 
                  features = "slingPseudotime_1",
                  reduction = "umap") +
  scale_color_viridis_c(option = "plasma") +
  ggtitle("Pseudotime - Lineage 1")

# Add trajectory curves
for(i in seq_along(curves)) {
  curve_data <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
  colnames(curve_data) <- c("UMAP_1", "UMAP_2")
  
  p2 <- p2 + geom_path(data = curve_data, 
                       aes(x = UMAP_1, y = UMAP_2),
                       size = 1.5, 
                       color = "black",
                       linetype = "dashed")
}

print(p2)
ggsave("outputs/figures/slingshot_pseudotime_lineage1.pdf", 
       p2, width = 10, height = 8)

# If multiple lineages, plot all
if(n_lineages > 1) {
  for(i in 2:n_lineages) {
    p_temp <- FeaturePlot(combined_integrated_final, 
                          features = paste0("slingPseudotime_", i),
                          reduction = "umap") +
      scale_color_viridis_c(option = "plasma") +
      ggtitle(paste("Pseudotime - Lineage", i))
    
    print(p_temp)
    ggsave(paste0("outputs/figures/slingshot_pseudotime_lineage", i, ".pdf"), 
           p_temp, width = 10, height = 8)
  }
}

# ============================================================
# STEP 6: DATASET COMPARISON ALONG PSEUDOTIME
# ============================================================

# Plot 3: Split by dataset to compare trajectories
p3 <- DimPlot(combined_integrated_final, 
              reduction = "umap", 
              group.by = "dataset",
              split.by = "dataset") +
  ggtitle("Trajectory by Dataset")

# Add curves to both panels
for(i in seq_along(curves)) {
  curve_data <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
  colnames(curve_data) <- c("UMAP_1", "UMAP_2")
  
  p3 <- p3 + geom_path(data = curve_data, 
                       aes(x = UMAP_1, y = UMAP_2),
                       size = 1.5, 
                       color = "black")
}

print(p3)
ggsave("outputs/figures/slingshot_trajectory_by_dataset.pdf", 
       p3, width = 14, height = 6)

# ============================================================
# STEP 7: GENE EXPRESSION ALONG PSEUDOTIME
# ============================================================

# Plot genes of interest along pseudotime
genes_of_interest <- c("TOP2A", "KIFC1",  # Proliferation
                       "SOX9", "MOG", "MAG", "ENPP6",  # 3D/Myelination
                       "OGN", "CYP1B1", "ASPN")  # 2D/Metabolic

# For primary lineage
lineage_1_data <- data.frame(
  pseudotime = combined_integrated_final$slingPseudotime_1,
  dataset = combined_integrated_final$dataset,
  cluster = combined_integrated_final$seurat_clusters
)

# Add gene expression
for(gene in genes_of_interest) {
  if(gene %in% rownames(combined_integrated_final)) {
    lineage_1_data[[gene]] <- GetAssayData(combined_integrated_final, 
                                           assay = "RNA", 
                                           layer = "data")[gene, ]
  }
}

# Remove cells not in this lineage (NA pseudotime)
lineage_1_data <- lineage_1_data[!is.na(lineage_1_data$pseudotime), ]

# Plot each gene
library(ggplot2)
library(cowplot)

gene_plots <- list()
for(gene in genes_of_interest) {
  if(gene %in% colnames(lineage_1_data)) {
    p <- ggplot(lineage_1_data, aes(x = pseudotime, y = .data[[gene]])) +
      geom_point(aes(color = dataset), alpha = 0.3, size = 0.5) +
      geom_smooth(aes(color = dataset), method = "loess", se = TRUE, span = 0.5) +
      scale_color_manual(values = c("GSE115011" = "#008080", "GSE146373" = "#FF7F50")) +
      theme_classic() +
      labs(title = gene, x = "Pseudotime", y = "Expression") +
      theme(legend.position = "bottom")
    
    gene_plots[[gene]] <- p
  }
}

# Combine plots
combined_gene_plot <- plot_grid(plotlist = gene_plots, ncol = 3)
print(combined_gene_plot)

ggsave("outputs/figures/slingshot_gene_expression_pseudotime.pdf", 
       combined_gene_plot, width = 15, height = 12)

# ============================================================
# STEP 8: STATISTICAL TESTING - DO DATASETS DIFFER IN PSEUDOTIME?
# ============================================================

# Test if datasets occupy different pseudotime positions
wilcox_test <- wilcox.test(
  lineage_1_data$pseudotime[lineage_1_data$dataset == "GSE115011"],
  lineage_1_data$pseudotime[lineage_1_data$dataset == "GSE146373"]
)

print("Wilcoxon test - Dataset pseudotime difference:")
print(wilcox_test)

# Calculate median pseudotime per dataset
median_pseudotime <- tapply(lineage_1_data$pseudotime, 
                            lineage_1_data$dataset, 
                            median, na.rm = TRUE)
print("Median pseudotime by dataset:")
print(median_pseudotime)

# Violin plot of pseudotime distribution
p_violin <- ggplot(lineage_1_data, aes(x = dataset, y = pseudotime, fill = dataset)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.alpha = 0.3) +
  scale_fill_manual(values = c("GSE115011" = "#008080", "GSE146373" = "#FF7F50")) +
  theme_classic() +
  labs(title = "Pseudotime Distribution by Dataset",
       subtitle = paste("Wilcoxon p =", signif(wilcox_test$p.value, 3)),
       x = "Dataset", y = "Pseudotime") +
  theme(legend.position = "none")

print(p_violin)
ggsave("outputs/figures/slingshot_pseudotime_distribution.pdf", 
       p_violin, width = 6, height = 6)

# ============================================================
# STEP 9: HEATMAP OF GENE EXPRESSION ALONG PSEUDOTIME
# ============================================================
BiocManager::install("ComplexHeatmap")


library(ComplexHeatmap)
library(circlize)

# Bin cells by pseudotime
lineage_1_data$pseudotime_bin <- cut(
  lineage_1_data$pseudotime,
  breaks = quantile(lineage_1_data$pseudotime,
                    probs = seq(0, 1, length.out = 21),
                    na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

# Calculate mean expression per bin per dataset
heatmap_data <- list()
for(gene in genes_of_interest) {
  if(gene %in% colnames(lineage_1_data)) {
    for(ds in c("GSE115011", "GSE146373")) {
      subset_data <- lineage_1_data[lineage_1_data$dataset == ds, ]
      mean_expr <- tapply(subset_data[[gene]], 
                          subset_data$pseudotime_bin, 
                          mean, na.rm = TRUE)
      
      heatmap_data[[paste(gene, ds, sep = "_")]] <- mean_expr
    }
  }
}

# Convert to matrix
heatmap_matrix <- do.call(rbind, heatmap_data)
colnames(heatmap_matrix) <- paste("Bin", 1:ncol(heatmap_matrix))


# Create heatmap
pdf("outputs/figures/slingshot_gene_heatmap_pseudotime.pdf", 
    width = 12, height = 8)

Heatmap(heatmap_matrix,
        name = "Expression",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        col = colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 
                           median(heatmap_matrix, na.rm = TRUE), 
                           max(heatmap_matrix, na.rm = TRUE)),
                         c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "Mean Expression"))

dev.off()

# ============================================================
# STEP 10: SAVE RESULTS
# ============================================================

# Save updated Seurat object with pseudotime
saveRDS(combined_integrated_final, 
        "outputs/seurat_integrated_with_slingshot.rds")

# Save slingshot object
saveRDS(sds, 
        "outputs/slingshot_trajectory_object.rds")

# Create summary statistics table
trajectory_summary <- data.frame(
  Dataset = c("GSE115011", "GSE146373"),
  N_cells = c(sum(lineage_1_data$dataset == "GSE115011"),
              sum(lineage_1_data$dataset == "GSE146373")),
  Median_Pseudotime = as.numeric(median_pseudotime),
  Mean_Pseudotime = tapply(lineage_1_data$pseudotime, 
                           lineage_1_data$dataset, 
                           mean, na.rm = TRUE),
  SD_Pseudotime = tapply(lineage_1_data$pseudotime, 
                         lineage_1_data$dataset, 
                         sd, na.rm = TRUE)
)

write.csv(trajectory_summary, 
          "outputs/slingshot_trajectory_summary.csv",
          row.names = FALSE)

print("Slingshot analysis complete!")
print(trajectory_summary)