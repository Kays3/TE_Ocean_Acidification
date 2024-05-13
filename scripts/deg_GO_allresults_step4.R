# GO collection visualize
# Author KD
# DESeq2 for 72 transcriptomes
# Apoly
# done normalization
# May 13, 2024
# Author KD




library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(reshape2)  # For melting the data matrix
library(dplyr)
library(readr)
library(stringr)  # For string manipulation



getwd()


#gdrive beholder
setwd("data/")



file_path <- "GO_terms/"  # Update this to your directory containing the CSV files
file_list <- list.files(path = file_path, pattern = "*.csv", full.names = TRUE)



# Read CSV files and select/rename necessary columns
data_frames <- lapply(file_list, function(file) {
  data <- read_csv(file)
  data <- data %>% 
    select(ID, p.adjust) %>%  # Adjust 'p_value' if your column name differs
    rename(!!paste0("p_value_", gsub("^.*/|\\.csv$", "", file)) := p.adjust)  # Rename p-value column based on file name
})


# Merge all data frames on 'GO_Term'
merged_data <- Reduce(function(x, y) full_join(x, y, by = "ID"), data_frames)

# Replace NA with zero in P_value columns
merged_data[grepl("p_value", names(merged_data))] <- lapply(merged_data[grepl("p_value", names(merged_data))], function(x) ifelse(is.na(x), 1, x))

# Filter GO terms with P-value less than 0.05 in at least one comparison
merged_data <- merged_data %>%
  filter(if_any(starts_with("P_value"), ~ .x < 0.05))


# Load the data
go_data <- read_csv("reference_files/final_merged_GO_terms.csv")
url_GO_annot<-"reference_files/ultimate_apoly_GO_annot.txt"

GO_annot<-read_delim(url_GO_annot)
GO_annot2<- GO_annot %>% distinct(go, .keep_all=TRUE)

colnames(GO_annot2)<-c("ID","gene_id","function")

go_data2<-merge(go_data,GO_annot2)
go_data3<-go_data2[1:20,]

str(go_data)


# Extracting only the necessary columns for the heatmap
heatmap_data <- select(go_data, starts_with("P_value"))
str(heatmap_data)


# Apply -log10 transformation to numeric columns, adding a small constant to avoid taking log10(0)
heatmap_data[] <- lapply(heatmap_data, function(x) if(is.numeric(x)) -log10(x + 1e-10) else x)

# Now convert the entire dataframe to a matrix for the heatmap, assuming all columns now are appropriate for the plot
#heatmap_data_matrix <- as.matrix(heatmap_data)
rownames(heatmap_data) <- go_data2$`function`  # Set GO terms as row names


str(heatmap_data)  # Adjust 'heatmap_data' to the name of your actual dataset

heatmap_data_matrix<-heatmap_data[1:20,]

rownames(heatmap_data_matrix) <- go_data3$`function`  # Set GO terms as row names

# Plotting the heatmap
pheatmap(heatmap_data, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Heatmap of -log10 Transformed P-values")



# Compute the minimum P-value across each row (GO term)
heatmap_data$min_p_value <- apply(heatmap_data, 1, min, na.rm = TRUE)
rownames(heatmap_data) <- go_data2$`function`  # Set GO terms as row names

# Order the data by min_p_value and select the top 20 GO terms
top_terms <- heatmap_data %>%
  arrange(min_p_value) %>%
  slice(1:20)

# Prepare the matrix for plotting
top_terms_matrix <- as.matrix(top_terms[, -which(names(top_terms) %in% c("GO_Term", "min_p_value"))])
rownames(top_terms_matrix) <- rownames(top_terms)

#plot


# Example of simplifying column names for the heatmap
# Assume original column names are like "Treatment_Group1_Control", "Treatment_Group2_Acute", etc.
short_labels <- c("acute_sensitive_devo_sensitive         ",
                  "acute_tolerant_devo_tolerant           ",
                  "control_sensitive_acute_sensitive      ",
                  "control_sensitive_devo_sensitive       ",
                  "control_tolerant_acute_tolerant        ",
                  "control_tolerant_control_sensitive     ",
                  "control_tolerant_developmental_tolerant",
                  "control_tolerant_trans_tolerant        ",
                  "trans_sensitive_acute_sensitive        ",
                  "trans_sensitive_devo_sensitive         ",
                  "trans_tolerant_acute_tolerant          ",
                  "trans_tolerant_devo_tolerant           "
                  )  # Simplified labels

full_labels <- c("p_value_acute_sensitive_devo_sensitive         ",
                 "p_value_acute_tolerant_devo_tolerant           ",
                 "p_value_control_sensitive_acute_sensitive      ",
                 "p_value_control_sensitive_devo_sensitive       ",
                 "p_value_control_tolerant_acute_tolerant        ",
                 "p_value_control_tolerant_control_sensitive     ",
                 "p_value_control_tolerant_developmental_tolerant",
                 "p_value_control_tolerant_trans_tolerant        ",
                 "p_value_trans_sensitive_acute_sensitive        ",
                 "p_value_trans_sensitive_devo_sensitive         ",
                 "p_value_trans_tolerant_acute_tolerant          ",
                 "p_value_trans_tolerant_devo_tolerant           "
                )  # Full descriptions

# Create a mapping of short labels to full descriptions
label_description <- setNames(full_labels, short_labels)

# Prepare the data with simplified labels (make sure to match these to your actual data structure)
#top_terms_matrix <- as.matrix(top_terms[, -which(names(top_terms) %in% c("GO_Term", "min_p_value"))])
colnames(top_terms_matrix) <- short_labels  # Replace column names with short labels
#rownames(top_terms_matrix) <- top_terms$GO_Term



# Prepare the data
data_melted <- melt(as.matrix(top_terms_matrix))
colnames(data_melted) <- c("GO_Term", "Condition", "Value")

# Define a threshold for significant values, for example:
threshold_value <- 2  # Adjust this based on your specific needs, such as a -log10 p-value threshold


# Create the heatmap using ggplot2
ggplot(data_melted, aes(x = Condition, y = GO_Term, fill = Value)) +
  geom_tile() +  # Create the tiles for the heatmap
  scale_fill_gradient2(low = "blue", high = "black", mid = "white", midpoint = median(data_melted$Value, na.rm = TRUE)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Heatmap of Top 20 GO Terms Significance",
    subtitle = "Visualization of -log10 transformed P-values across conditions",
    x = "Pairwise Comparisons",
    y = "GO Terms",
    fill = "Significance"
  )
ggsave("refined_heatmap_GO_Terms.png", plot = last_plot(), width = 10, height = 8, dpi = 300)


#rm(list=ls())
