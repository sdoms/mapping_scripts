setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/mills/")
library(tidyverse)
library(stringr)
library(readxl)
load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_03032021_SW.Rdata")
input_DNA <- all_dna_rna_SW %>% filter(dna.rna=="DNA")
input_RNA<- all_dna_rna_SW %>% filter(dna.rna=="RNA")
input_SW<- all_dna_rna_SW

all_genes<- unique(as.character(unlist(sapply(input_SW$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
all.peak.genes <- unique(all_dna_rna_SW$closest_gene)
all_genes_RNA<- unique(as.character(unlist(sapply(input_RNA$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
all_genes_DNA<- unique(as.character(unlist(sapply(input_DNA$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
write.table(all_genes_DNA, "../../../all_genes_intervals/SW/all_genes_DNA_SW.txt", col.names = F, row.names = F, quote = F)
write.table(all_genes_RNA, "../../../all_genes_intervals/SW/all_genes_RNA_SW.txt", col.names = F, row.names = F, quote = F)
write.table(all_genes, "../../../all_genes_intervals/SW/all_genes_DNA_RNA_SW.txt", col.names = F, row.names = F, quote = F)
write.table(all.peak.genes, "../../../all_genes_intervals/SW/all_peak_genes_DNA_RNA_SW.txt", col.names = F, row.names = F, quote = F)

all_genes <- splus2R::upperCase(all_genes)
# snps ####

library(tidyverse)
library(stringr)
library(readxl)

infile <- read.csv("all_markers_RNA_DNA-with-genes_SW.csv", sep=";")
snp_genes<- as.character(unlist(sapply(infile$all_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
snp_genes<- unique(snp_genes)
preced_genes <-as.character(unlist(sapply(infile$all_preceding_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
preced_genes<- unique(preced_genes)
follow_genes <- as.character(unlist(sapply(infile$all_following_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
follow_genes<- unique(follow_genes)

all_genes <- c(snp_genes, preced_genes, follow_genes)
all_genes<- unique(na.omit(all_genes))
all_genes<- all_genes[-1] #remove NA

snp_genes <- snp_genes[c(-1,-2)] # remoce NA
write(all_genes, "all_genes_snps_SW.txt")
write(snp_genes, "all_snp_genes_snps_SW.txt")

