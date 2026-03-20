library(EnhancedVolcano)
library(ggplot2)

# Load DE results
de_results <- read.csv("TableS4_dataset_differential_expression.csv", row.names = 1)

# Define genes to highlight
highlight_genes <- c(
  # GSE115011-enriched
  "MBP", "MOG", "ENPP6", "ATP6", "COX1", "COX2", "COX3", "CYTB",
  # GSE146373-enriched
  "OGN", "ASPN", "IL8", "CXCL1", "CYP1B1"
)

# Create volcano plot
p <- EnhancedVolcano(de_results,
                     lab = rownames(de_results),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     selectLab = highlight_genes,
                     title = 'GSE115011 vs GSE146373',
                     subtitle = 'Genome-wide Differential Expression',
                     pCutoff = 0.05,
                     FCcutoff = 1.0,  # log2FC > 1 = 2-fold
                     pointSize = 2.0,
                     labSize = 5.0,
                     labCol = 'black',
                     labFace = 'bold',
                     boxedLabels = TRUE,
                     col = c('gray30', 'forestgreen', 'royalblue', 'red2'),
                     colAlpha = 0.5,
                     drawConnectors = TRUE,
                     widthConnectors = 0.5,
                     max.overlaps = 20) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

ggsave("Fig3B_volcano_clean.pdf", p, width = 12, height = 10, dpi = 300)