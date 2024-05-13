# GO for 12 pairwise comparisons
# Author KD
# DESeq2 for 72 transcriptomes
# Apoly
# done normalization
# May 13, 2024
# Author KD





library(gplots)
library(ggvenn)
library(stringr)
library(dplyr)
library(tidyverse)
library(VennDiagram)
library(qvalue)
library(ggplot2)
library(svMisc)
library(wesanderson)
# load the library
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(readr)
library(dplyr)
library(tidyverse)


getwd()


#gdrive beholder
setwd("data/deg_results")

#insr Done
Cont_tol_vs_cont_sens<-read.delim("1_control_tolerant_control_sensitive_464_DEG_may2024.txt",header = T)



#load GO temrs for enrichers
#unzip it first
url_GO_annot<-"~data/ultimate_apoly_GO_annot.txt"

GO_annot<-read_delim(url_GO_annot)

disease2gene=GO_annot[, c("go", "gene_id")]
disease2name=GO_annot[, c("go", "function")]





qval.001 <- c(Cont_tol_vs_cont_sens$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
head(summary(enr.qval.001))


#dotplot(enr.qval.001)
#ggsave("3974_TEs_GO.png",width = 10,height = 10, units = "in")


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "Cont_tol_vs_cont_sens.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)




control_tolerant_acute_tolerant<-read.delim("5_control_tolerant_acute_tolerant_6501_DEG_may2024.txt",header = T)




qval.001 <- c(control_tolerant_acute_tolerant$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "5_control_tolerant_acute_tolerant.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)



control_tolerant_developmental_tolerant<-read.delim("6_control_tolerant_developmental_tolerant_3829_DEG_may2024.txt",header = T)




qval.001 <- c(control_tolerant_developmental_tolerant$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "control_tolerant_developmental_tolerant.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)



#control_tolerant_trans_tolerant

control_tolerant_trans_tolerant<-read.delim("7_control_tolerant_trans_tolerant_151_DEG_may2024.txt",header = T)




qval.001 <- c(control_tolerant_trans_tolerant$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "control_tolerant_trans_tolerant.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)


#8_acute_tolerant_devo_tolerant_456_DEG_may2024

acute_tolerant_devo_tolerant<-read.delim("8_acute_tolerant_devo_tolerant_456_DEG_may2024.txt",header = T)




qval.001 <- c(acute_tolerant_devo_tolerant$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "acute_tolerant_devo_tolerant.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)




#9_trans_tolerant_devo_tolerant_2320_DEG_may2024

trans_tolerant_devo_tolerant<-read.delim("9_trans_tolerant_devo_tolerant_2320_DEG_may2024.txt",header = T)




qval.001 <- c(trans_tolerant_devo_tolerant$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "trans_tolerant_devo_tolerant.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)


#10_trans_tolerant_acute_tolerant_5363_DEG_may2024

trans_tolerant_acute_tolerant<-read.delim("10_trans_tolerant_acute_tolerant_5363_DEG_may2024.txt",header = T)




qval.001 <- c(trans_tolerant_acute_tolerant$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "trans_tolerant_acute_tolerant.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)




#11_control_sensitive_acute_sensitive_928_DEG_may2024

control_sensitive_acute_sensitive<-read.delim("11_control_sensitive_acute_sensitive_928_DEG_may2024.txt",header = T)




qval.001 <- c(control_sensitive_acute_sensitive$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
#head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "control_sensitive_acute_sensitive.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)




#12_control_sensitive_devo_sensitive_1240_DEG_may2024

control_sensitive_devo_sensitive<-read.delim("12_control_sensitive_devo_sensitive_1240_DEG_may2024.txt",header = T)




qval.001 <- c(control_sensitive_devo_sensitive$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
#head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "control_sensitive_devo_sensitive.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)



#14_acute_sensitive_devo_sensitive_716_DEG_may2024

acute_sensitive_devo_sensitive<-read.delim("14_acute_sensitive_devo_sensitive_716_DEG_may2024.txt",header = T)




qval.001 <- c(acute_sensitive_devo_sensitive$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
#head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "acute_sensitive_devo_sensitive.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)



#15_trans_sensitive_devo_sensitive_1650_DEG_may2024

trans_sensitive_devo_sensitive<-read.delim("15_trans_sensitive_devo_sensitive_1650_DEG_may2024.txt",header = T)




qval.001 <- c(trans_sensitive_devo_sensitive$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
#head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "trans_sensitive_devo_sensitive.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)



#16_trans_sensitive_acute_sensitive_1203_DEG_may2024

trans_sensitive_acute_sensitive<-read.delim("16_trans_sensitive_acute_sensitive_1203_DEG_may2024.txt",header = T)




qval.001 <- c(trans_sensitive_acute_sensitive$gene_id)

enr.qval.001 = enricher(gene = qval.001,
                        
                        TERM2GENE=disease2gene, TERM2NAME=disease2name)
#head(summary(enr.qval.001))


up.tab.tol = enr.qval.001@result

write.table(up.tab.tol, file = "trans_sensitive_acute_sensitive.csv", sep = ",", quote = F, 
            row.names = F, col.names = T)




rm(list=ls())


