Single-Cell RNA-seq Analysis of iPSC-Derived Oligodendrocyte Differentiation

📌 Overview

This project investigates how culture dimensionality (2D vs 3D) influences the differentiation of human induced pluripotent stem cells (iPSCs) into oligodendrocytes using single-cell RNA sequencing (scRNA-seq).

We performed an integrative analysis of 1,178 high-quality cells from:

3D spheroid culture (GSE115011)

2D monolayer culture (GSE146373)

The workflow includes data integration, clustering, annotation, differential expression, pathway analysis, and trajectory inference.

⚙️ Workflow Summary

Quality control & filtering

Data normalization & scaling

Dataset integration (CCA-based, Seurat)

Clustering (Louvain algorithm)

Cell type annotation (SingleR)

Differential expression analysis

Gene Ontology (GO) enrichment

Pseudotime trajectory analysis (Slingshot)

Visualization: UMAP, heatmaps, volcano plots, Venn diagrams

📂 Repository Structure
🧪 Core Analysis Scripts

SINGLE-CELL RNA-SEQ INTEGRATION & ANALYSIS PIPELINE.R → Main pipeline

Libraries.R → Required packages

Quality Control Metrics.R → QC and filtering

SINGLECELLEXPERIMENT.R → Data structure handling

🔍 Downstream Analysis

Pseudotime.R / Trajectory_Analysis.R → Lineage inference

HEATMAP_50.R → Top marker visualization

volcano clean.R → Differential expression plots

venn.R → Marker gene overlap

💾 Data

GSE115011_seurat.rds → Processed Seurat object
GSEGSE146373_seurat.rds → Processed Seurat object

📊 Outputs

Figures: UMAPs, heatmaps, volcano plots

Tables: Differentially expressed genes, annotations

🔬 Key Findings

3D cultures:

Neural/glial dominance (~65% neurons, ~32% astrocytes)

Strong oligodendrocyte marker expression (MBP, MOG, ENPP6)

Enriched in myelination pathways

2D cultures:

Predominantly mesenchymal/fibroblast-like cells (~84.7%)

No oligodendrocyte lineage commitment

Enriched in lipid metabolism pathways

Only 2/78 marker genes (TOP2A, KIFC1) shared

Significant transcriptomic divergence (>1000-fold changes)

Distinct developmental trajectories confirmed by pseudotime

📊 Conclusion

This study demonstrates that culture architecture is a key determinant of lineage fidelity:

3D systems → neurogenic & oligodendrocyte differentiation

2D systems → mesenchymal drift
