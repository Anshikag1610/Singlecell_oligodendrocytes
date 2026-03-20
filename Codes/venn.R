# Install from CRAN
install.packages("ggVennDiagram")

# Then load it
library(ggVennDiagram)
library(ggVennDiagram)
library(ggplot2)

# Load gene lists
common <- read.csv("common_marker_genes.csv")$gene
unique_115011 <- read.csv("unique_genes_GSE115011.csv")$gene
unique_146373 <- read.csv("unique_genes_GSE146373.csv")$gene

# Create gene list object
gene_lists <- list(
  "GSE115011\n3D Spheroid" = c(unique_115011, common),
  "GSE146373\n2D Monolayer" = c(unique_146373, common)
)

# Create Venn diagram
p <- ggVennDiagram(gene_lists,
                   label_alpha = 0,
                   category.names = names(gene_lists),
                   set_color = c("#008080", "#FF7F50")) +
  scale_fill_gradient(low = "white", high = "gray90") +
  scale_color_manual(values = c("#006666", "#CC5533")) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Marker Gene Overlap",
       subtitle = "Only 2.6% shared (2 of 78 genes: TOP2A, KIFC1)")

ggsave("Fig3A_Venn_enhanced.pdf", p, width = 8, height = 8, dpi = 300)