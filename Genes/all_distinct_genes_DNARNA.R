setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/mills/")
library(tidyverse)
library(stringr)
library(readxl)
input_DNA <- read_excel("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/all_sig_DNA_core.xlsx")
input_RNA<- read_excel("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/all_sig_RNA_core.xlsx")

input_DNA$DNAorRNA <- "DNA"
input_RNA$DNAorRNA<- "RNA"
SW_threshold <- (0.05/32625)/118.2177
input_DNA<- input_DNA %>% 
  filter(peak.P.type< SW_threshold)
input_RNA<- input_RNA %>% 
  filter(peak.P.type< SW_threshold)# input_SW<- read_csv2("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/results_DNA_RNA_SW.csv")
input_SW <- rbind(input_DNA,input_RNA)
all_genes<- unique(as.character(unlist(sapply(input_SW$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
all_genes_RNA<- unique(as.character(unlist(sapply(input_RNA$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
all_genes_DNA<- unique(as.character(unlist(sapply(input_DNA$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T)))))
write.table(all_genes_DNA, "../../../all_genes_intervals/SW/all_genes_DNA_SW.txt", col.names = F, row.names = F, quote = F)
write.table(all_genes_RNA, "../../../all_genes_intervals/SW/all_genes_RNA_SW.txt", col.names = F, row.names = F, quote = F)
write.table(all_genes, "../../../all_genes_intervals/SW/all_genes_DNA_RNA_SW.txt", col.names = F, row.names = F, quote = F)

all_genes <- splus2R::upperCase(all_genes)
# snps ####
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")
library(tidyverse)
library(stringr)
library(readxl)
# 
# infile <- read.csv("../../Shared/all_markers_RNA_DNA-with-genes_GW.csv", sep=";")
# snp_genes<- as.character(unlist(sapply(infile$all_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
# snp_genes<- unique(snp_genes)
# preced_genes <-as.character(unlist(sapply(infile$all_preceding_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
# preced_genes<- unique(preced_genes)
# follow_genes <- as.character(unlist(sapply(infile$all_following_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
# follow_genes<- unique(follow_genes)
# 
# all_genes <- c(snp_genes, preced_genes, follow_genes)
# all_genes<- unique(na.omit(all_genes))
# 
# 
# snp_genes <- na.omit(snp_genes) # remoce NA
# write(all_genes, "all_genes_snps_GW.txt")
# write(snp_genes, "all_snp_genes_snps_GW.txt")

##### study-wide ######

# ASV<- read_excel("./DNA/summary_tables/study_wide_significant/gwscan/SV_results_DNA_SW.xlsx", skip=2)
# ASV_genes<- as.character(unlist(sapply(ASV$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# ASV_genes<- unique(ASV_genes)
# 
# all_genes<- ASV_genes
# 
# genus<- read_excel("./DNA/summary_tables/study_wide_significant/gwscan/results_genus_DNA_SW.xlsx", sheet=1)
# genus_genes<- as.character(unlist(sapply(genus$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# genus_genes <- unique(genus_genes)
# 
# all_genes<- unique(append(all_genes, genus_genes))
# 
# 
# family<- read_excel("./DNA/summary_tables/study_wide_significant/gwscan/family_results_DNA_SW.xlsx", sheet="all")
# family_genes<- as.character(unlist(sapply(family$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# family_genes <- unique(family_genes)
# 
# all_genes<- unique(append(all_genes, family_genes))
# 
# order<- read_excel("./DNA/summary_tables/study_wide_significant/gwscan/order_results_DNA_SW.xlsx", sheet=1)
# order_genes<- as.character(unlist(sapply(order$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# order_genes <- unique(order_genes)
# 
# all_genes<- unique(append(all_genes, order_genes))
# 
# class<- read_excel("./DNA/summary_tables/study_wide_significant/gwscan/class_results_DNA_SW.xlsx", sheet="all")
# class_genes<- as.character(unlist(sapply(class$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# class_genes <- unique(class_genes)
# 
# all_genes<- unique(append(all_genes, class_genes))
# 
# phylum<- read_excel("./DNA/summary_tables/study_wide_significant/gwscan/phylum_results_DNA_SW.xlsx", sheet=2)
# phylum_genes<- as.character(unlist(sapply(phylum$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# phylum_genes <- unique(phylum_genes)
# 
# all_genes<- unique(append(all_genes, phylum_genes))
# 
# write(ASV_genes,file="ASV_genes_DNA.txt")
# write(genus_genes,file="genus_genes_DNA.txt")
# write(family_genes,file="family_genes_DNA.txt")
# write(order_genes,file="order_genes_DNA.txt")
# write(class_genes,file="class_genes_DNA.txt")
# write(phylum_genes,file="phylum_genes_DNA.txt")
# write(all_genes,file="all_genes_DNA.txt")
# all_genes_DNA<-all_genes


####### RNA ######
# ASV<- read_excel("./RNA/summary_tables/study_wide_significant/gwscan/SW_results_SV.xlsx",skip=1, sheet="all")
# ASV_genes<- as.character(unlist(sapply(ASV$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# ASV_genes<- unique(ASV_genes)
# 
# all_genes<- ASV_genes
# 
# genus<- read_excel("./RNA/summary_tables/study_wide_significant/gwscan/SW_results_genus.xlsx",skip=1, sheet="all")
# genus_genes<- as.character(unlist(sapply(genus$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# genus_genes <- unique(genus_genes)
# 
# all_genes<- unique(append(all_genes, genus_genes))
# 
# 
# family<- read_excel("./RNA/summary_tables/study_wide_significant/gwscan/SW_results_family.xlsx", skip=1,sheet="all")
# family_genes<- as.character(unlist(sapply(family$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# family_genes <- unique(family_genes)
# 
# all_genes<- unique(append(all_genes, family_genes))
# 
# order<- read_excel("./RNA/summary_tables/study_wide_significant/gwscan/SW_results_order.xlsx",skip=1, sheet="all")
# order_genes<- as.character(unlist(sapply(order$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# order_genes <- unique(order_genes)
# 
# all_genes<- unique(append(all_genes, order_genes))
# 
# class<- read_excel("./RNA/summary_tables/study_wide_significant/gwscan/SW_results_class.xlsx")
# class_genes<- as.character(unlist(sapply(class$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# class_genes <- unique(class_genes)
# 
# all_genes<- unique(append(all_genes, class_genes))
# 
# phylum<- read_excel("./RNA/summary_tables/study_wide_significant/gwscan/SW_results_phylum.xlsx")
# phylum_genes<- as.character(unlist(sapply(phylum$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
# phylum_genes <- unique(phylum_genes)
# 
# all_genes<- unique(append(all_genes, phylum_genes))
# 
# write(ASV_genes,file="ASV_genes_RNA.txt")
# write(genus_genes,file="genus_genes_RNA.txt")
# write(family_genes,file="family_genes_RNA.txt")
# write(order_genes,file="order_genes_RNA.txt")
# write(class_genes,file="class_genes_RNA.txt")
# write(phylum_genes,file="phylum_genes_RNA.txt")
# write(all_genes,file="all_genes_RNA.txt")
# 
# DNA_RNA_genes <- unique(append(all_genes, all_genes_DNA))
# write(DNA_RNA_genes,file="all_genes_DNA_RNA.txt")
####################

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")
library(tidyverse)
library(stringr)
library(readxl)

infile <- read.csv("../../Shared/all_markers_RNA_DNA-with-genes_SW.csv", sep=",")
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

