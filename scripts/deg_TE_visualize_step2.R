# DESeq2 for 72 transcriptomes
# Apoly
# visualize TE related Pogo 
# May 13 2024
# Author KD

library(DESeq2)
library(ggrepel)
library(pheatmap)
library(ggsignif)

library(ggplot2)
library(tools)
library(ggtext)  # for rich text formatting
#library(EnhancedVolcano)
#library(limma)
#library(qvalue)
library(dplyr)
library(tidyverse)

options(width=100)

#rm(list=ls())
getwd()
setwd("data/")



# Reload the DESeqDataSet object from the RDS file
dds <- readRDS("dds_object_72transcripts_new.rds")


# Normalize counts
normalized_counts <- counts(dds, normalized=TRUE)

# Let's assume you know the row name or ID of the significant gene.
# For example, "POGO element"
# Apoly008489
gene_id <- "Apoly008489" 

# Extract the normalized counts for this specific gene
gene_data <- data.frame(
  gene_expression = normalized_counts[gene_id, ],
  sample_info = colData(dds)
)

# Since colData might have more info, let's assume you have 'parental_phenotype' and 'treatment' as part of the experimental design
gene_data$parental_phenotype <- gene_data$sample_info$parental_phenotype
gene_data$treatment <- gene_data$sample_info$treatment

# Create the box plot with an annotation for log scale
pogo_plot <- ggplot(gene_data, aes(x = sample_info.treatment, y = gene_expression, fill = sample_info.parental_phenotype)) +
  geom_boxplot() +
  labs(title = "Expression of Pogo transposable element with KRAB domain\nAcross Phenotypes and Treatments",
       x = "Treatment Group",
       y = "Normalized Gene Expression",
       fill = "Parental Phenotype") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold", margin = margin(10, 0, 10, 0)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  scale_y_log10() +  # Apply log transformation to y-axis
  annotate("text", x = Inf, y = Inf, label = "Scale: Log10", hjust = 1.1, vjust = 2, size = 3, color = "red")

print(pogo_plot)

# Save the plot with appropriate size and resolution for Nature
ggsave("pogo_levels_72transcriptomes.png", plot = pogo_plot, width = 8, height = 6, dpi = 300, device = 'png')

# make final plots and interpret the results


rm(list = ls())

