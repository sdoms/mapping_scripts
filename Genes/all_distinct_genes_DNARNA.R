setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/")
library(tidyverse)
library(stringr)
library(readxl)
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_22022022_SW_pval_corr.Rdata")
input_DNA <- all_dna_rna_SW %>% filter(dna.rna=="DNA")
input_RNA<- all_dna_rna_SW %>% filter(dna.rna=="RNA")
input_SW<- all_dna_rna_SW

all_genes<- unique(as.character(unlist(sapply(input_SW$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
all.peak.genes <- unique(all_dna_rna_SW$closest_gene)
all_genes_RNA<- unique(as.character(unlist(sapply(input_RNA$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
all_genes_DNA<- unique(as.character(unlist(sapply(input_DNA$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
write.table(all_genes_DNA, "../../all_genes_intervals/SW/all_genes_DNA_SW_pval_corr.txt", col.names = F, row.names = F, quote = F)
write.table(all_genes_RNA, "../../all_genes_intervals/SW/all_genes_RNA_SW_pval_corr.txt", col.names = F, row.names = F, quote = F)
write.table(all_genes, "../../all_genes_intervals/SW/all_genes_DNA_RNA_SW_pval_corr.txt", col.names = F, row.names = F, quote = F)
write.table(all.peak.genes, "../../all_genes_intervals/SW/all_peak_genes_DNA_RNA_SW_pval_corr.txt", col.names = F, row.names = F, quote = F)

all_genes <- splus2R::upperCase(all_genes)
# snps ####

library(tidyverse)
library(stringr)
library(readxl)

infile <- read.csv("../../all_genes_snps/all_markers_RNA_DNA-with-genes_SW.csv", sep=";")
snp_genes<- as.character(unlist(sapply(infile$all_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
snp_genes<- unique(snp_genes)
preced_genes <-as.character(unlist(sapply(infile$all_preceding_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
preced_genes<- unique(preced_genes)
follow_genes <- as.character(unlist(sapply(infile$all_following_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
follow_genes<- unique(follow_genes)

all_genes <- c(snp_genes, preced_genes, follow_genes)
all_genes<- unique(na.omit(all_genes))
all_genes<- all_genes[-2] #remove NA

snp_genes <- snp_genes[c(-1,-3)] # remoce NA
write(all_genes, "all_genes_snps_SW_pval_corr.txt")
write(snp_genes, "all_snp_genes_snps_SW_pval_corr.txt")

