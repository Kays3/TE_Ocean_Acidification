# TE-DEG collection visualize
# Author KD
# DESeq2 for 72 transcriptomes
# Apoly
# May 13, 2024
# revised May 27, 2024 - TBD
# clustered along Contrasts, legend with whole numbers

# Author KD



getwd()


#gdrive beholder
setwd("TE_Ocean_Acidification/data/")

library(ggplot2)
library(reshape2)  # For melting your data frame if it's not already in long format
# Load the libraries
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)
library(pheatmap)


# Example data in wide format
# manually extracted from excell files of all DEG contrasts
# using filter funtion search for "transpos"

data <- read.csv("reference_files/TE_DEG_heatmap.csv")
data1<-as.data.frame(cbind(data$function_protein,data$Contrast))
# Assuming 'Contrast' is the column indicating the group each element belongs to


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

#need to try this
# Ensure the data is in matrix format
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# Ensure the data is in matrix format
count_matrix <- as.matrix(element_counts2[-1])  # Exclude any non-numeric columns if present
rownames(count_matrix) <- element_counts2[[1]]  # Set the first column as row names if it represents the element types

# Perform hierarchical clustering
row_clust <- hclust(dist(count_matrix))
col_clust <- hclust(dist(t(count_matrix)))

# Convert to long format for ggplot2
count_matrix_long <- as.data.frame(count_matrix) %>%
  rownames_to_column("TE_family") %>%
  pivot_longer(-TE_family, names_to = "sample_name", values_to = "count")

# Custom ordering based on presence of "tol" or "sens"
custom_order <- c(grep("tol", rownames(count_matrix), value = TRUE), 
                  grep("sens", rownames(count_matrix), value = TRUE), 
                  setdiff(rownames(count_matrix), c(grep("tol", rownames(count_matrix), value = TRUE), 
                                                    grep("sens", rownames(count_matrix), value = TRUE))))

# Extract ordered names based on clustering and custom order
ordered_TE_families <- intersect(custom_order, rownames(count_matrix)[row_clust$order])
ordered_conditions <- colnames(count_matrix)[col_clust$order]

# Factor the names based on the clustering and custom order
count_matrix_long$TE_family <- factor(count_matrix_long$TE_family, levels = ordered_TE_families)
count_matrix_long$sample_name <- factor(count_matrix_long$sample_name, levels = ordered_conditions)

# Create a new column for custom colors based on "tol" or "sens"
count_matrix_long$color_group <- ifelse(grepl("tol", count_matrix_long$TE_family), "Tolerant parents",
                                        ifelse(grepl("sens", count_matrix_long$TE_family), "Sensitive parents", "other"))

# Define custom colors
custom_colors <- c("tol" = "red", "sens" = "blue", "other" = "grey80")

# Create the heatmap with ggplot2
gg_heatmap <- ggplot(count_matrix_long, aes(x = sample_name, y = TE_family)) +
  geom_tile(aes(fill = count), color = "white") +
  geom_text(aes(label = ifelse(count == 0, "", round(count, 0))), size = 3, color = "black") +
  scale_fill_gradientn(colors = c("white", "lightblue", "blue", "red"), 
                       values = scales::rescale(c(0, 1, 2)), 
                       breaks = c(0, 1, 2), 
                       limits = c(0, 2), 
                       guide = "colorbar") +
  facet_wrap(~ color_group, scales = "free_y", ncol = 1, strip.position = "left") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text.y = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Heatmap of Transposable Elements Counts", x = "TE Family", y = "Fish Comparison", fill = "Count")

# Display the plot
print(gg_heatmap)

# Save the heatmap
ggsave("expression_heatmap_clustered.png", plot = gg_heatmap, width = 12, height = 10, dpi = 300)



# black and white further - no clustering
# Assuming 'data' is already loaded and transformed

data<-as.data.frame(count_matrix)

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

#clustered





# Save the heatmap to a file in a high-quality format
#ggsave("heatmap_transposable_elements.png", plot = heatmap_plot, width = 10, height = 6, device = "png")

#try



rm(list=ls())
