# TE-DEG collection visualize
# Author KD
# DESeq2 for 72 transcriptomes
# Apoly
# May 13, 2024
# Author KD



getwd()


#gdrive beholder
setwd("data/")

library(ggplot2)
library(reshape2)  # For melting your data frame if it's not already in long format
# Load the libraries
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)


# Example data in wide format
# manually extracted from excell files of all DEG contrasts
# using filter funtion search for "transpos"

data <- read.csv("reference_files/TE_DEG_heatmap.csv")
data1<-as.data.frame(cbind(data$function_protein,data$# Assuming 'ComparisonGroup' is the column indicating the group each element belongs to


element_counts2 <- data %>% count(function_protein,Contrast) %>% spread(key = function_protein, value = n, fill = 0) 
# Spread to wide format suitable for a heatmap
print(element_counts2)

# Ensure the data is in matrix format for pheatmap
count_matrix <- as.matrix(element_counts2[-1])  # Exclude any non-numeric columns if present
rownames(count_matrix) <- element_counts2[[1]]  # Set the first column as row names if it represents the element types

# Create the heatmap
pheatmap(count_matrix,
         color = colorRampPalette(c("white", "blue"))(50),  # Use a sequential color palette
         scale = "row",  # Normalize across rows
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Heatmap of Transposable Elements Counts")

#need to reverse

# Assuming 'data' is already loaded and transformed
library(reshape2)  # For melting the data matrix
library(ggplot2)

# Assuming 'count_matrix' is your matrix
# Convert the matrix to a long format data frame
data_long <- melt(as.matrix(count_matrix))

# Rename columns appropriately
colnames(data_long) <- c("Condition", "TransposableElement", "Count")

# Define the color scale
color_scale <- scale_fill_gradientn(colors = rev(colorRampPalette(c("white", "blue"))(50)))

# Create the heatmap
heatmap_plot <- ggplot(data_long, aes(x = Condition , y = TransposableElement , fill = Count)) +
  geom_tile() +  # Use tiles for the heatmap
  color_scale +  # Apply the color scale
  labs(
    title = "Heatmap of Transposable Elements Counts",
    x = "Transposable Element",
    y = "Condition",
    fill = "Count"
  ) +
  theme_minimal() +  # Nature-style theme
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  coord_fixed()  # Keep aspect ratio to 1 to make tiles square

# Print the plot
print(heatmap_plot)

#version2

# Define a suitable color scale
color_scale <- scale_fill_gradientn(colors = rev(colorRampPalette(c("black", "white"))(100)), 
                                    name = "Count",
                                    limits = c(min(data_long$Count), max(data_long$Count)),
                                    breaks = pretty(data_long$Count, n = 5))

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(data_long, aes(x = Condition, y = TransposableElement, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.1) +  # Add white borders for clarity
  color_scale +
  labs(title = "Heatmap of Differentially Expressed Transposable Elements \n Across Pairwise Comparisons",
       x = "Pairwise Comparisons",
       y = "Transposable Element") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  coord_fixed(ratio = 1 / 3)  # Adjust the aspect ratio

# Print the plot
print(heatmap_plot)


# Save the heatmap to a file in a high-quality format
#ggsave("heatmap_transposable_elements.png", plot = heatmap_plot, width = 10, height = 6, device = "png")

rm(list=ls())
