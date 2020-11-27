setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/")
library(readxl)
library(tidyverse)
ASV<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/SV_results_DNA.xlsx", skip = 3)
ASV_genes<- as.character(unlist(sapply(ASV$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
ASV_genes<- unique(ASV_genes)

all_genes<- ASV_genes

genus<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/genus_results_DNA.xlsx", sheet="all")
genus_genes<- as.character(unlist(sapply(genus$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
genus_genes <- unique(genus_genes)

all_genes<- unique(append(all_genes, genus_genes))


family<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/family_results_DNA.xlsx", sheet="all")
family_genes<- as.character(unlist(sapply(family$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
family_genes <- unique(family_genes)

all_genes<- unique(append(all_genes, family_genes))

order<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/order_results_DNA.xlsx", sheet="all")
order_genes<- as.character(unlist(sapply(order$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
order_genes <- unique(order_genes)

all_genes<- unique(append(all_genes, order_genes))

class<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/class_results_DNA.xlsx", sheet="all")
class_genes<- as.character(unlist(sapply(class$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
class_genes <- unique(class_genes)

all_genes<- unique(append(all_genes, class_genes))

phylum<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/phylum_results_DNA.xlsx", sheet="all")
phylum_genes<- as.character(unlist(sapply(phylum$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
phylum_genes <- unique(phylum_genes)

all_genes<- unique(append(all_genes, phylum_genes))

write(ASV_genes,file="ASV_genes_DNA.txt")
write(genus_genes,file="genus_genes_DNA.txt")
write(family_genes,file="family_genes_DNA.txt")
write(order_genes,file="order_genes_DNA.txt")
write(class_genes,file="class_genes_DNA.txt")
write(phylum_genes,file="phylum_genes_DNA.txt")
write(all_genes,file="all_genes_DNA.txt")
all_genes_DNA<-all_genes


####### RNA ######
ASV<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/SV_results_RNA.xlsx",skip=2, sheet="all")
ASV_genes<- as.character(unlist(sapply(ASV$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
ASV_genes<- unique(ASV_genes)

all_genes<- ASV_genes

genus<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/genus_results_RNA.xlsx", sheet="all")
genus_genes<- as.character(unlist(sapply(genus$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
genus_genes <- unique(genus_genes)

all_genes<- unique(append(all_genes, genus_genes))


family<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/family_results_RNA.xlsx", sheet="all")
family_genes<- as.character(unlist(sapply(family$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
family_genes <- unique(family_genes)

all_genes<- unique(append(all_genes, family_genes))

order<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/order_results_RNA.xlsx", sheet="all")
order_genes<- as.character(unlist(sapply(order$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
order_genes <- unique(order_genes)

all_genes<- unique(append(all_genes, order_genes))

class<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/class_results_RNA.xlsx")
class_genes<- as.character(unlist(sapply(class$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
class_genes <- unique(class_genes)

all_genes<- unique(append(all_genes, class_genes))

phylum<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/phylum_results_RNA.xlsx")
phylum_genes<- as.character(unlist(sapply(phylum$list_total_genes, FUN=function(x) strsplit(x, ",", fixed = T))))
phylum_genes <- unique(phylum_genes)

all_genes<- unique(append(all_genes, phylum_genes))

write(ASV_genes,file="ASV_genes_RNA.txt")
write(genus_genes,file="genus_genes_RNA.txt")
write(family_genes,file="family_genes_RNA.txt")
write(order_genes,file="order_genes_RNA.txt")
write(class_genes,file="class_genes_RNA.txt")
write(phylum_genes,file="phylum_genes_RNA.txt")
write(all_genes,file="all_genes_RNA.txt")

DNA_RNA_genes <- unique(append(all_genes, all_genes_DNA))
write(DNA_RNA_genes,file="all_genes_DNA_RNA.txt")


####################3

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")
library(tidyverse)
library(stringr)
library(readxl)

infile <- read.csv("../../Shared/all_markers_RNA_DNA-with-genes_SW.csv")
snp_genes<- as.character(unlist(sapply(infile$all_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
snp_genes<- unique(snp_genes)
preced_genes <-as.character(unlist(sapply(infile$all_preceding_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
preced_genes<- unique(preced_genes)
follow_genes <- as.character(unlist(sapply(infile$all_following_genes, FUN=function(x) strsplit(x, " | ", fixed = T))))
follow_genes<- unique(follow_genes)

all_genes <- c(snp_genes, preced_genes, follow_genes)
all_genes<- unique(all_genes)

write(all_genes, "all_genes_snps.txt")
write(snp_genes, "all_snp_genes_snps.txt")

