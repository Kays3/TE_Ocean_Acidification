# DESeq2 for 72 transcriptomes
# Apoly
# Z SCORE expression across all 12 pairwise comparison
# May 13, 2024
# 12 highly DEG-TE visualized
# Author KD

library(DESeq2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggtext)  # for rich text formatting
library(viridis)
library(dplyr)
library(tidyverse)

options(width=100)

#rm(list=ls())
getwd()
setwd("../data/")



# Reload the DESeqDataSet object from the RDS file
#dds <- readRDS("reference_files/dds_object_72transcripts_new.rds")
te_deg_all<-read.csv("reference_files/TE_DEG_heatmap.csv")

te_deg_all2<- te_deg_all %>% distinct(Apoly_gene, function_protein, .keep_all = TRUE)

te_deg_all3<- te_deg_all2 %>% 
  mutate(id = make.unique(as.character(function_protein)))

# Assuming 'parental_phenotype' and 'treatment' are the names of your two conditions in the colData
dds$combinedCondition <- paste(dds$parental_phenotype, dds$treatment, sep = "_")

# Assume the metadata file contains columns "sample_name" and "condition"
metadata <- data.frame(
  sample_name = colnames(dds),
  condition = colData(dds)$combinedCondition  # Replace "condition" with your specific column name
)
# Specify the list of desired transcripts
desired_transcripts <- te_deg_all2$Apoly_gene
# Replace with your actual list

# Extract normalized counts
norm_counts <- counts(dds, normalized = TRUE)


#colnames(norm_counts) <- dds$combinedCondition

# Filter for specific transcripts
filtered_counts <- norm_counts[rownames(norm_counts) %in% desired_transcripts, ]



# Set the column names based on sample names
#colnames(filtered_counts) <- metadata$sample_name

# zscore-transform the data
log10_counts <- t(apply(filtered_counts,1,scale))

colnames(log10_counts) <- metadata$sample_name


# Replace rownames using a named vector
new_rownames <- setNames(te_deg_all3$id, te_deg_all2$Apoly_gene)
# Match and replace rownames
rownames(log10_counts) <- new_rownames[rownames(log10_counts)]

# Create a data frame with zscore-transformed expression values grouped by condition
log10_counts_grouped <- as.data.frame(log10_counts) %>%
  rownames_to_column("transcript") %>%
  pivot_longer(-transcript, names_to = "sample_name", values_to = "expression") %>%
  left_join(metadata, by = "sample_name") %>%
  group_by(transcript, condition) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_expression) %>%
  column_to_rownames("transcript")


# Convert the data frame to long format for ggplot2
log10_counts_long <- as.data.frame(log10_counts_grouped) %>%
  rownames_to_column("transcript") %>%
  pivot_longer(-transcript, names_to = "condition", values_to = "mean_expression")


# Extract the order of the clustered columns
ordered_sample_names <- colnames(log2_counts)[col_clust$order]

# Factor the sample names based on the clustering
log2_counts_long$sample_name <- factor(log2_counts_long$sample_name, levels = ordered_sample_names)


# Compute row and column clustering for VST data
row_clust_z <- hclust(dist(log10_counts))
col_clust_z <- hclust(dist(t(log10_counts)))


# Create the heatmap with ggplot2
gg_heatmap <- ggplot(log10_counts_long, aes(x = condition, y = transcript, fill = mean_expression)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "Expression Patterns of TE-related Transcripts", x = "Condition", y = "Transcript", fill = "Z-score Expression")+
  facet_wrap(~ factor(transcript, levels = rownames(log10_counts)[row_clust_z$order]))


# Display the plot
print(gg_heatmap)

# Save the ggplot2 heatmap as a high-quality PNG
#ggsave("z-score_ALL_DEG_TE_heatmap.png", plot = gg_heatmap, width = 10, height = 8, dpi = 300)

####
#trial

# Convert the data frame to long format for ggplot2
log10_counts_long <- as.data.frame(log10_counts) %>%
  rownames_to_column("transcript") %>%
  pivot_longer(-transcript, names_to = "sample_name", values_to = "expression") %>%
  left_join(as.data.frame(colData(dds)) %>%
              rownames_to_column("sample_name"),
            by = "sample_name")

# Aggregate data by condition and transcript
log10_counts_grouped <- log10_counts_long %>%
  group_by(transcript, combinedCondition) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = combinedCondition, values_from = mean_expression) %>%
  column_to_rownames("transcript")

# Perform hierarchical clustering
row_clust <- hclust(dist(log10_counts_grouped))
col_clust <- hclust(dist(t(log10_counts_grouped)))

# Convert to long format for ggplot2
log10_counts_long <- as.data.frame(log10_counts_grouped) %>%
  rownames_to_column("transcript") %>%
  pivot_longer(-transcript, names_to = "condition", values_to = "mean_expression")

# Extract ordered names based on clustering
ordered_transcripts <- rownames(log10_counts_grouped)[row_clust$order]
ordered_conditions = colnames(log10_counts_grouped)[col_clust$order]

# Factor the names based on the clustering
log10_counts_long$transcript <- factor(log10_counts_long$transcript, levels = ordered_transcripts)
log10_counts_long$condition <- factor(log10_counts_long$condition, levels = ordered_conditions)

# Create the heatmap with ggplot2
gg_log10_heatmap <- ggplot(log10_counts_long, aes(x = condition, y = transcript, fill = mean_expression)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "navy", mid = "grey", high = "firebrick", midpoint = median(log10_counts_long$mean_expression, na.rm = TRUE)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "Differential Expression of TE-related transcripts in the brain", x = "Condition", y = "Transcript", fill = "Z-score Expression")

# Display the plot
print(gg_log10_heatmap)

# Save the heatmap
ggsave("Zscore_expression_heatmap_grouped_clustered_bw.png", plot = gg_log10_heatmap, width = 12, height = 10, dpi = 300)



#rm(list = ls())

