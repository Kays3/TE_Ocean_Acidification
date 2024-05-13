# DESeq2 for 72 transcriptomes
# Apoly
# done normalization
# May 13
# Author KD

library(DESeq2)
library(ggrepel)
library(pheatmap)
library(ggplot2)
library(tools)
library(dplyr)
library(tidyverse)

options(width=100)

#rm(list=ls())
getwd()

datadir<-"data/htseq_apoly"

refere<-read.table("data/reference_files/neo_gene_ref_apoly.txt", header = T)

new_refere<-read.delim("data/reference_files/apoly_gene_protein_refere.txt", header = T)

phdata<-data.frame(fname=list.files(path=datadir,pattern="*Aligned.sortedByCoord.out.counts.txt"),stringsAsFactors=FALSE)
head(phdata)
dim(phdata)


phdata <- phdata %>% transmute(sample=substr(fname,1,7),fname)
head(phdata)


# Read in the sample description file
sample_desc <- read.table("data/reference_files/Apoly_CO2_samplesheet.2021.txt", header=TRUE, sep="\t")

sample_desc <- sample_desc[-c(123, 128, 216),]

NAMES = gsub("Aligned.sortedByCoord.out.counts.txt","", grep(".txt", dir(path = datadir ), value=TRUE))
DATA=data.frame("file"=NAMES,sample=substr(NAMES,1,6),stringsAsFactors=FALSE)
DATA$NUMERIC <-c(1:72)

colnames(DATA)[1]<-"RNAseq.fastq"
colnames(DATA)[2]<-"sample_Mirko"


DATA = merge(DATA,sample_desc, merge="RNAseq.fastq", all.x=TRUE)
DATA$"rep" = substr(DATA$"RNAseq.fastq",7,7)

DATA$"Treat_pheno" =paste(DATA$"treatment",DATA$"parental_phenotype",sep="_")

minidata<-as.data.frame(cbind(DATA$RNAseq.fastq,DATA$treatment,DATA$parental_phenotype))
colnames(minidata)<-c("sample","treatment","parental_phenotype")

phdata1<-merge(phdata,minidata, merge="sample", all.x=TRUE)

phdata1 <- phdata1 %>% mutate(treatment=as.factor(treatment),
                              parental_phenotype=as.factor(parental_phenotype),
                              md5=tools::md5sum(file.path(datadir,fname)))
head(phdata1)



phdata1$treatment <- relevel(phdata1$treatment, "control")
levels(phdata1$treatment)

levels(phdata1$parental_phenotype)

#make dds pbject for deseq2

dds<-DESeqDataSetFromHTSeqCount(sampleTable=phdata1,directory=datadir,
                                design= ~ parental_phenotype + treatment + parental_phenotype:treatment)



colData(dds)
dds@design
#resultsNames(dds)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)
dim(counts(dds))
head(counts(dds),3)

gene_means <- rowMeans(norm_counts)
gene_sds <- apply(norm_counts, 1, sd)  # Calculate standard deviation for each gene

# Create a logical vector where TRUE represents genes where SD < Mean
valid_genes <- gene_sds < gene_means

# Subset the DESeqDataSet to keep only valid genes
dds <- dds[valid_genes, ]


# Filter to keep genes with counts greater than 10 in at least 3 samples
n_samples <- 3  # Minimum number of samples
threshold <- 10  # Minimum count threshold

# Apply the filter
dds <- dds[rowSums(counts(dds) >= threshold) >= n_samples, ]


dds<-DESeq(dds)


# see your results

resultsNames(dds)

# get the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
colnames(mod_mat)

# Define coefficient vectors for each condition
control_tolerant <- colMeans(mod_mat[dds$treatment == "control" & dds$parental_phenotype == "tolerant", ])
control_sensitive <- colMeans(mod_mat[dds$treatment == "control" & dds$parental_phenotype == "sensitive", ])

acute_tolerant <- colMeans(mod_mat[dds$treatment == "acute" & dds$parental_phenotype == "tolerant", ])
acute_sensitive <- colMeans(mod_mat[dds$treatment == "acute" & dds$parental_phenotype == "sensitive", ])

devo_tolerant <- colMeans(mod_mat[dds$treatment == "developmental" & dds$parental_phenotype == "tolerant", ])
devo_sensitive <- colMeans(mod_mat[dds$treatment == "developmental" & dds$parental_phenotype == "sensitive", ])

trans_tolerant <- colMeans(mod_mat[dds$treatment == "transgenerational" & dds$parental_phenotype == "tolerant", ])
trans_sensitive <- colMeans(mod_mat[dds$treatment == "transgenerational" & dds$parental_phenotype == "sensitive", ])

#results(dds, contrast=c("group", "controlsensitive", "controltolerant"))
#comparison to do
# 1 control tolerant vs control sensitive
# 2 trans tol vs trans sens
# 3 dev tol vs def sens
# 4 acute tol vs acute sens
# 5 control tol vs acute tol
# 6 control tol vs dev tol
# 7 control tol vs trans tol
# 8 dev tol vs acute tol
# 9 trasn tol vs dev tol
# 10 trans tol vs acute tol
# 11 control sens vs acute sens
# 12 control sens vs dev sens
# 13 control sens vs trans sens
# 14 dev sens vs acute sens
# 15 trasn sens vs dev sens
# 16 trans sens vs acute sens

# no need next since no biological relevance - mixed up 4 groups per Parental Phenotype
# 17 all TOLERANT VS ALL SENSITIVE 


### Set thresholds as in Schunter et al 2018
padj.cutoff <- 0.05
lfc.cutoff <- 0.5


res1 <- results(dds, contrast = control_tolerant - control_sensitive)
#res <- results(dds)



## Defined contrasts, extract results table and shrink log2 fold changes
# control is always is REFERENCE - thus this comparison is control tolerant vs control sensitive
#more stringent test giving that control and sensitive values are used as reference

#1
#contrast_kd=c("group", "controlsensitive", "controltolerant")


res_table1_tb <- res1 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

resOrdered <- res_table1_tb[order(res_table1_tb$pvalue),]



sig1 <- res_table1_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sig1)

colnames(sig1)[1]<-"Apoly_gene"
new_control_tol_v_sense<-merge(sig1,new_refere, by = "Apoly_gene")


#do exploratory GO
# report genes per each comparison


plotMA(res1)

res_tableOE_tb <- new_control_tol_v_sense %>% 
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.5)

write.table(res_tableOE_tb, file = "1_control_tolerant_control_sensitive_464_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

# DEGs list for Control tolerant vs Control Sensitive
# search for TE 

# !!!! found POGO TE

# Identify the point to highlight

# Example criteria for "significant" might be padj < 0.05 and abs(log2FoldChange) >= 1

# Highlight label - assuming 'Protein_X' is the one to be highlighted

# Calculate significance and sort by padj to get the top 10 significant points
res_tableOE_tb$significant <- with(res_tableOE_tb, padj < 0.05 & abs(log2FoldChange) >= 0.5)
top_significant <- res_tableOE_tb[order(res_tableOE_tb$padj), ][1:10, ]

# Highlighted label
highlight_label <- "Pogo transposable element with KRAB domain"

# Plot
volcano_plot <- ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj), colour = significant)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_point(data = subset(res_tableOE_tb, function_protein == highlight_label), 
             aes(x = log2FoldChange, y = -log10(padj)), 
             colour = "gold", size = 4, shape = 17) +
  geom_text_repel(data = top_significant,
                  aes(label = function_protein, x = log2FoldChange, y = -log10(padj)),
                  size = 3.5, 
                  max.overlaps = Inf) +
  geom_text_repel(data = subset(res_tableOE_tb, function_protein == highlight_label),
                  aes(label = function_protein, x = log2FoldChange, y = -log10(padj)),
                  color = "black", 
                  size = 4, 
                  nudge_x = 1,  # Adjust as needed to move the label
                  nudge_y = 1,  # Adjust as needed to move the label
                  segment.color = 'grey50') +
  scale_colour_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(
    title = "Control tolerant vs Control sensitive fish (643 DEG)",
    subtitle = "significant: p-value adjusted < 0.05 and log2 fold change > 0.5",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value"
  )

# Display the plot
print(volcano_plot)

#make changes:
#add TE transcripts to legend
# POGO LFC 0.93 
# add if upregulated of downregulated

# 
# Save the plot with appropriate size and resolution for Nature
ggsave("may2024/04292024_Plot_Control_tol_vs_sens_.png", plot = volcano_plot, width = 8, height = 6, dpi = 300, device = 'png')

length(which(res_tableOE_tb$significant == "TRUE"))




## Define contrasts, extract results table and shrink log2 fold changes Control tolerant vs control sensitive

### 2 transgene tol vs sens

res_tableKD_trans_tol_sens <- results(dds, contrast=trans_tolerant - trans_sensitive)
res_tableKD_trans_tol_sens_tb <- res_tableKD_trans_tol_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_trans_tol_sens <- res_tableKD_trans_tol_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_trans_tol_sens)


colnames(sigKD_trans_tol_sens)[1]<-"Apoly_gene"
new_trans<-merge(sigKD_trans_tol_sens,new_refere, by = "Apoly_gene")

write.table(new_trans, file = "may2024/2_transgenerationaltolerant_transgenerationalsensitive_20_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

#plot volcano
#no exploratory GO
# report genes per each comparison


plotMA(res_tableKD_trans_tol_sens)

res_tableOE_trans_tol_sens <- new_trans %>% 
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.5)


#20 genes (vs Celia Paper 2018 - 152)

# no TEs found
#plot

# Calculate significance and sort by padj to get the top 10 significant points
res_tableOE_trans_tol_sens$significant <- with(res_tableOE_trans_tol_sens, padj < 0.05 & abs(log2FoldChange) >= 0.5)
top_significant_trans_tol_sens <- res_tableOE_trans_tol_sens[order(res_tableOE_trans_tol_sens$padj), ][1:10, ]

# Highlighted label
  highlight_label2 <- "Transcription factor NF-E2 45 kDa subunit"


# Plot
volcano_plot2 <- ggplot(res_tableOE_trans_tol_sens, aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_point(data = subset(res_tableOE_trans_tol_sens, function_protein == highlight_label2), 
             aes(x = log2FoldChange, y = -log10(padj)), 
             colour = "gold", size = 4, shape = 17) +
  geom_text_repel(data = top_significant_trans_tol_sens,
                  aes(label = function_protein, x = log2FoldChange, y = -log10(padj)),
                  size = 3.5, 
                  max.overlaps = Inf) +
  geom_text_repel(data = subset(res_tableOE_trans_tol_sens, function_protein == highlight_label),
                  aes(label = function_protein, x = log2FoldChange, y = -log10(padj)),
                  color = "black", 
                  size = 4, 
                  nudge_x = 1,  # Adjust as needed to move the label
                  nudge_y = 1,  # Adjust as needed to move the label
                  segment.color = 'grey50') +
  scale_colour_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(
    title = "Transgenerational tolerant vs Transgenerational sensitive fish",
    subtitle = "p-value adjusted < 0.05 and minimum log2 fold change > 0.3 (66 DEG)",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value"
  )

# Display the plot
print(volcano_plot2)


# Save the plot with appropriate size and resolution for Nature
#ggsave("Volcano_Plot_Transgenerational_tol_vs_sens_66DEG.png", plot = volcano_plot, width = 8, height = 6, dpi = 300, device = 'png')

#write.table(res_tableOE_trans_tol_sens, file = "2_trans_tolerant_trans_sensitive_66_DEG_April2024.txt", row.names = F, sep = "\t", quote = F)


### 3 developmental tol vs sens

res_tableKD_devo_tol_sens <- results(dds, contrast=devo_tolerant - devo_sensitive)
res_tableKD_devo_tol_sens_tb <- res_tableKD_devo_tol_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_developmental_tol_sens <- res_tableKD_devo_tol_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_developmental_tol_sens)

#23 genes (vs Celia Paper 2018 - 359)
#xplore


colnames(sigKD_developmental_tol_sens)[1]<-"Apoly_gene"
new_developmental<-merge(sigKD_developmental_tol_sens,new_refere, by = "Apoly_gene")



#plot volcano
#do exploratory GO
# report genes per each comparison



res_tableOE_developmental_tol_sens <- new_developmental %>% 
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.5)

# no TEs found
#no plot
write.table(res_tableOE_developmental_tol_sens, file = "3_developmentaltolerant_developmentalsensitive_23_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

### 4 acute tol vs sens


res_tableKD_acute_tol_sens <- results(dds, contrast=acute_tolerant - acute_sensitive)
res_tableKD_acute_tol_sens_tb <- res_tableKD_acute_tol_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_acute_tol_sens <- res_tableKD_acute_tol_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_acute_tol_sens)

#16 genes (vs Celia Paper 2018 - 144)



colnames(sigKD_acute_tol_sens)[1]<-"Apoly_gene"
new_acute<-merge(sigKD_acute_tol_sens,new_refere, by = "Apoly_gene")


# no TEs found


write.table(new_acute, file = "4_acutetolerant_acutesensitive_16_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)



# 5 control tol vs acute tol


res_tableKD_control_tol_acute_tol <- results(dds, contrast=control_tolerant - acute_tolerant)

res_tableKD_control_tol_acute_tol_tb <- res_tableKD_control_tol_acute_tol %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_control_tol_acute_tol <- res_tableKD_control_tol_acute_tol_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_control_tol_acute_tol)



#explore for TEs...

colnames(sigKD_control_tol_acute_tol)[1]<-"Apoly_gene"
new_control_tol_acute_tol<-merge(sigKD_control_tol_acute_tol,new_refere, by = "Apoly_gene")



#6501 genes (vs Celia Paper 2018 - 3669)
#6 transcripts with TE associated terms


write.table(new_control_tol_acute_tol, file = "5_control_tolerant_acute_tolerant_6501_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)


# 6 control tol  vs developmental tol

res_tableKD_control_tol_devo_tol <- results(dds, contrast=control_tolerant - devo_tolerant)

res_tableKD_control_tol_devo_tol_tb <- res_tableKD_control_tol_devo_tol %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_control_tol_devo_tol <- res_tableKD_control_tol_devo_tol_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_control_tol_devo_tol)

#explore for TEs...



colnames(sigKD_control_tol_devo_tol)[1]<-"Apoly_gene"
new_control_tol_devo_tol<-merge(sigKD_control_tol_devo_tol,new_refere, by = "Apoly_gene")


#3829 genes (vs Celia Paper 2018 - 1142)

write.table(new_control_tol_devo_tol, file = "6_control_tolerant_developmental_tolerant_3829_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)



# 7 control tol vs transgene tol

#342 genes (vs Celia Paper 2018 - 173)




res_tableKD_control_tol_trans_tol <- results(dds, contrast=control_tolerant - trans_tolerant)

res_tableKD_control_tol_trans_tol_tb <- res_tableKD_control_tol_trans_tol %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_control_tol_trans_tol <- res_tableKD_control_tol_trans_tol_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_control_tol_trans_tol)

#explore for TEs...



colnames(sigKD_control_tol_trans_tol)[1]<-"Apoly_gene"
new_control_tol_trans_tol<-merge(sigKD_control_tol_trans_tol,new_refere, by = "Apoly_gene")


write.table(new_control_tol_trans_tol, file = "7_control_tolerant_trans_tolerant_151_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)


# 8 acute vs developmental tol



res_tableKD_acute_tol_devo_tol <- results(dds, contrast=acute_tolerant - devo_tolerant)

res_tableKD_acute_tol_devo_tol_tb <- res_tableKD_acute_tol_devo_tol %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_acute_tol_devo_tol <- res_tableKD_acute_tol_devo_tol_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_acute_tol_devo_tol)

#explore for TEs...



colnames(sigKD_acute_tol_devo_tol)[1]<-"Apoly_gene"
new_acute_tol_devo_tol<-merge(sigKD_acute_tol_devo_tol,new_refere, by = "Apoly_gene")




#456 genes (vs Celia Paper 2018 - 1519)
write.table(new_acute_tol_devo_tol, file = "8_acute_tolerant_devo_tolerant_456_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)


# 9 trans vs developmental tol



res_tableKD_trans_tol_devo_tol <- results(dds, contrast=trans_tolerant - devo_tolerant)

res_tableKD_trans_tol_devo_tol_tb <- res_tableKD_trans_tol_devo_tol %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_trans_tol_devo_tol <- res_tableKD_trans_tol_devo_tol_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_trans_tol_devo_tol)

#explore for TEs...



colnames(sigKD_trans_tol_devo_tol)[1]<-"Apoly_gene"
new_trans_tol_devo_tol<-merge(sigKD_trans_tol_devo_tol,new_refere, by = "Apoly_gene")


write.table(new_trans_tol_devo_tol, file = "9_trans_tolerant_devo_tolerant_2320_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

#2320 genes (vs Celia Paper 2018 - 2497)

# 10 trans vs acute tol



res_tableKD_trans_tol_acute_tol <- results(dds, contrast=trans_tolerant - acute_tolerant)

res_tableKD_trans_tol_acute_tol_tb <- res_tableKD_trans_tol_acute_tol %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_trans_tol_acute_tol <- res_tableKD_trans_tol_acute_tol_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_trans_tol_acute_tol)

#explore for TEs...



colnames(sigKD_trans_tol_acute_tol)[1]<-"Apoly_gene"
new_trans_tol_acute_tol<-merge(sigKD_trans_tol_acute_tol,new_refere, by = "Apoly_gene")


write.table(new_trans_tol_acute_tol, file = "10_trans_tolerant_acute_tolerant_5363_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

#5363 genes (vs Celia Paper 2018 - not evaluated)


# 11 acute vs control sensitive


res_tableKD_control_sens_acute_sens <- results(dds, contrast=control_sensitive - acute_sensitive)

res_tableKD_control_sens_acute_sens_tb <- res_tableKD_control_sens_acute_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_control_sens_acute_sens <- res_tableKD_control_sens_acute_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_control_sens_acute_sens)



#explore for TEs...

colnames(sigKD_control_sens_acute_sens)[1]<-"Apoly_gene"
new_control_sens_acute_sens<-merge(sigKD_control_sens_acute_sens,new_refere, by = "Apoly_gene")



#928 genes (vs Celia Paper 2018 - 2010)


write.table(new_control_sens_acute_sens, file = "11_control_sensitive_acute_sensitive_928_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

# 12 developmental vs control sensitive

# 12 control sens vs devo sens


res_tableKD_control_sens_devo_sens <- results(dds, contrast=control_sensitive - devo_sensitive)

res_tableKD_control_sens_devo_sens_tb <- res_tableKD_control_sens_devo_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_control_sens_devo_sens <- res_tableKD_control_sens_devo_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_control_sens_devo_sens)



#explore for TEs...

colnames(sigKD_control_sens_devo_sens)[1]<-"Apoly_gene"
new_control_sens_devo_sens<-merge(sigKD_control_sens_devo_sens,new_refere, by = "Apoly_gene")



write.table(new_control_sens_devo_sens, file = "12_control_sensitive_devo_sensitive_1240_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

#1240 genes (vs Celia Paper 2018 - 2590)


### 13 control sens vs trans sens


res_tableKD_control_sens_trans_sens <- results(dds, contrast=control_sensitive - trans_sensitive)

res_tableKD_control_sens_trans_sens_tb <- res_tableKD_control_sens_trans_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_control_sens_trans_sens <- res_tableKD_control_sens_trans_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_control_sens_trans_sens)



#explore for TEs...

colnames(sigKD_control_sens_trans_sens)[1]<-"Apoly_gene"
new_control_sens_trans_sens<-merge(sigKD_control_sens_trans_sens,new_refere, by = "Apoly_gene")



write.table(new_control_sens_trans_sens, file = "13_control_sensitive_trans_sensitive_0_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)


#0 genes (vs Celia Paper 2018 - 62)



#14 acute vs developmental sensitive



res_tableKD_acute_sens_devo_sens <- results(dds, contrast=acute_sensitive - devo_sensitive)

res_tableKD_acute_sens_devo_sens_tb <- res_tableKD_acute_sens_devo_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_acute_sens_devo_sens <- res_tableKD_acute_sens_devo_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_acute_sens_devo_sens)

#explore for TEs...



colnames(sigKD_acute_sens_devo_sens)[1]<-"Apoly_gene"
new_acute_sens_devo_sens<-merge(sigKD_acute_sens_devo_sens,new_refere, by = "Apoly_gene")




write.table(new_acute_sens_devo_sens, file = "14_acute_sensitive_devo_sensitive_716_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)


#716 genes (vs Celia Paper 2018 - 1500)



#15 trans vs developmental sensitive


# 15 trans vs developmental sens



res_tableKD_trans_sens_devo_sens <- results(dds, contrast=trans_sensitive - devo_sensitive)

res_tableKD_trans_sens_devo_sens_tb <- res_tableKD_trans_sens_devo_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_trans_sens_devo_sens <- res_tableKD_trans_sens_devo_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_trans_sens_devo_sens)

#explore for TEs...



colnames(sigKD_trans_sens_devo_sens)[1]<-"Apoly_gene"
new_trans_sens_devo_sens<-merge(sigKD_trans_sens_devo_sens,new_refere, by = "Apoly_gene")




write.table(new_trans_sens_devo_sens, file = "15_trans_sensitive_devo_sensitive_1650_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

#1650 genes (vs Celia Paper 2018 - 2507)



#16 trans vs acute sensitive


res_tableKD_trans_sens_acute_sens <- results(dds, contrast=trans_sensitive - acute_sensitive)

res_tableKD_trans_sens_acute_sens_tb <- res_tableKD_trans_sens_acute_sens %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()



sigKD_trans_sens_acute_sens <- res_tableKD_trans_sens_acute_sens_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

dim(sigKD_trans_sens_acute_sens)

#explore for TEs...



colnames(sigKD_trans_sens_acute_sens)[1]<-"Apoly_gene"
new_trans_sens_acute_sens<-merge(sigKD_trans_sens_acute_sens,new_refere, by = "Apoly_gene")




write.table(new_trans_sens_acute_sens, file = "16_trans_sensitive_acute_sensitive_1203_DEG_may2024.txt", row.names = F, sep = "\t", quote = F)

#1203 genes (vs Celia Paper 2018 - not evaluated)


#17 all tolerant vs all sensitive doesnt make sense.

# make final plots and interpret the results
# Save the DESeqDataSet object to an RDS file
saveRDS(dds, file = "dds_object_72transcripts_new.rds")



rm(list = ls())
