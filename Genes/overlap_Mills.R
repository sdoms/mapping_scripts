# This script will find the genes from SW threshold and the overlap with diff expressed protein mills 2020
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

# ASV_DNA<- read_excel("../../../../DNA/summary_tables/study_wide_significant/gwscan/SV_results_DNA_SW.xlsx", skip=2)
# genus_DNA<- read_excel("../../../../DNA/summary_tables/study_wide_significant/gwscan/results_genus_DNA_SW.xlsx", sheet=1)
# family_DNA<- read_excel("../../../../DNA/summary_tables/study_wide_significant/gwscan/family_results_DNA_SW.xlsx", sheet="all")
# order_DNA<- read_excel("../../../../DNA/summary_tables/study_wide_significant/gwscan/order_results_DNA_SW.xlsx", sheet=1)
# class_DNA<- read_excel("../../../../DNA/summary_tables/study_wide_significant/gwscan/class_results_DNA_SW.xlsx", sheet="all")
# phylum_DNA<- read_excel("../../../../DNA/summary_tables/study_wide_significant/gwscan/phylum_results_DNA_SW.xlsx", sheet=2)
# 
# genes_DNA<-rbind(genus_DNA,family_DNA)
# genes_DNA<- rbind(family_DNA,genes_DNA)
# genes_DNA<-rbind(order_DNA,genes_DNA)
# genes_DNA<-rbind(class_DNA,genes_DNA)
# genes_DNA<-rbind(phylum_DNA,genes_DNA)
# genes_DNA$Genus<-"NA"
# genes_DNA<- rbind(genes_DNA,ASV_DNA)
# genes_DNA$DNAorRNA <- "DNA"
# 
# ASV_RNA<- read_excel("../../../../RNA/summary_tables/study_wide_significant/gwscan/SW_results_SV.xlsx",skip=1, sheet="all")
# genus_RNA<- read_excel("../../../../RNA/summary_tables/study_wide_significant/gwscan/SW_results_genus.xlsx",skip=1, sheet="all")
# family_RNA<- read_excel("../../../../RNA/summary_tables/study_wide_significant/gwscan/SW_results_family.xlsx", skip=1,sheet="all")
# order_RNA<- read_excel("../../../../RNA/summary_tables/study_wide_significant/gwscan/SW_results_order.xlsx",skip=1, sheet="all")
# class_RNA<- read_excel("../../../../RNA/summary_tables/study_wide_significant/gwscan/SW_results_class.xlsx")
# phylum_RNA<- read_excel("../../../../RNA/summary_tables/study_wide_significant/gwscan/SW_results_phylum.xlsx")
# 
# genes_RNA<-rbind(order_RNA,family_RNA)
# genes_RNA<- rbind(class_RNA,genes_RNA)
# genes_RNA<-rbind(phylum_RNA,genes_RNA)
# genes_RNA$Genus<-"NA"
# genes_RNA<- rbind(genes_RNA,ASV_RNA)
# genes_RNA<- rbind(genes_RNA,genus_RNA)
# 
# genes_RNA$DNAorRNA <- "RNA"
# 
# genes_DNA_RNA<-rbind(genes_DNA,genes_RNA) %>% 
#   drop_na(name)
save(input_SW, file="../../../all_genes_intervals/SW/genes_DNA_RNA_intervals_SW.Rdata")

con_si<- unique(read_xlsx("../../../overexpressed_protein_GF_con.xlsx",sheet="CON_SI", col_names = c("gene", "score"))[,1])
gf_si <- unique(read_xlsx("../../../overexpressed_protein_GF_con.xlsx",sheet="GF_SI", col_names = c("gene", "score"))[,1])
gf_col <- unique(read_xlsx("../../../overexpressed_protein_GF_con.xlsx",sheet="GF", col_names = c("gene")))
con_col <- unique(na.omit(read_xlsx("../../../overexpressed_protein_GF_con.xlsx",sheet="con", col_names = c("gene"))))

genes_con <- unique(rbind(con_si, con_col))
genes_gf <- unique(rbind(gf_si, gf_col))
mills_genes <- unique((rbind(genes_con, genes_gf)))


# conventional overexpressed
overlap_con_intervals<-data.frame()
for (gene in genes_con$gene){
  xx<-str_which(all_genes, paste0("\\b",gene, "\\b$"))
  if (length(xx)>0){
    result<-input_SW[str_which(input_SW$list_total_genes,regex(paste0("\\b",gene, "\\b[, ]"),ignore_case = T)),]
    result$gene <- gene
    overlap_con_intervals<- rbind(overlap_con_intervals,result)
  }
}
overlap_con_intervals <- overlap_con_intervals[,c(33, 1:32)]
write_excel_csv(overlap_con_intervals, "intervals/overlap_con_intervals_SW_curated.csv")

# 5 mb
overlap_con_intervals_filt5mb <- overlap_con_intervals %>% 
  filter(as.numeric(as.character(length))<5)

overlap_con_intervals_filt_genes<- overlap_con_intervals_filt5mb %>% 
  distinct(gene)
overlap_con_intervals_genes<- overlap_con_intervals %>% 
  distinct(gene)

# germ-free overexpressed

overlap_gf_intervals<-data.frame()
for (gene in genes_gf$gene){
  xx<-str_which(all_genes, paste0("\\b",gene, "\\b$"))
  if (length(xx)>0){
    result<-input_SW[str_which(input_SW$list_total_genes,regex(paste0("\\b",gene, "\\b[, ]"),ignore_case = T)),]
    result$gene <- gene
    overlap_gf_intervals<- rbind(overlap_gf_intervals,result)
  }
}
overlap_gf_intervals <- overlap_gf_intervals[,c(33, 1:32)]
write_csv2(overlap_gf_intervals, "intervals/overlap_gf_intervals_SW_curated.csv")

overlap_all_intervals <- rbind(overlap_con_intervals, overlap_gf_intervals)
write_csv2(overlap_all_intervals, "intervals/overlap_all_intervals_SW_curated.csv")
# read in the curated file, removed wrong searches
# overlap_curated <- read.csv2("overlap_all_intervals_curated.csv")

# count the number of taxa
overlap_all_intervals$minP<- apply(overlap_all_intervals[,c("peak.P.type", "add.P_peak_snp", "dom.P_peak_snp")],1, FUN = function(x) min(as.numeric(as.character(x)), na.rm=T))

pvals_all_overlap_intervals <- overlap_all_intervals%>% 
  dplyr::group_by(gene, trait) %>% 
  dplyr::mutate(P_vals=paste(P.type, collapse=(" | "))) %>% 
  distinct(gene, trait, P_vals, minP, length, DNAorRNA) 

taxa_all_overlap_intervals<- pvals_all_overlap_intervals%>% 
  dplyr::group_by(gene, DNAorRNA) %>% 
  add_tally() %>% 
  dplyr::mutate(taxa=paste(trait,"(", P_vals, ")", collapse = " | ")) %>% 
  dplyr::mutate(minP=min(minP), min_int=min(length)) %>% 
distinct(gene, taxa, minP, n, min_int, DNAorRNA) %>% ungroup() %>% ungroup

write_csv(taxa_all_overlap_intervals, "intervals/taxa_all_overlap_intervals_curated_SW.csv")
write_csv(pvals_all_overlap_intervals, "intervals/pvals_all_overlap_intervals_curated_SW.csv")

genes_overlap_intervals_curated <- taxa_all_overlap_intervals %>% 
  distinct(gene)
write(genes_overlap_intervals_curated$gene, "intervals/genes_overlap_intervals_all_curated_SW.txt")
taxa_dna_rna<- taxa_all_overlap_intervals %>% 
  group_by(gene) %>% 
add_tally(name="dna_rna_count") %>% 
  dplyr::mutate(DNA_both_RNA=ifelse(dna_rna_count==2, "both", DNAorRNA), n_tot=sum(n), 
                taxa_tot=paste(DNAorRNA, ": ", taxa, collapse="; "), minimum_interval=min(min_int), min_P=min(minP)) %>% 
  distinct(gene, DNA_both_RNA, n_tot, taxa_tot, min_P, minimum_interval)
write_csv(taxa_dna_rna, "intervals/taxa_dna_rna_overlap_intervals_curated_SW.csv")
# 5 mb
overlap_gf_intervals_filt5mb <- overlap_gf_intervals %>% 
  filter(as.numeric(as.character(length))<5)

overlap_gf_intervals_filt_genes<- overlap_gf_intervals_filt5mb %>% 
  distinct(gene)
overlap_gf_intervals_genes<- overlap_gf_intervals %>% 
  distinct(gene)

tot_overlap_intervals <- unique(rbind(overlap_con_intervals_genes, overlap_gf_intervals_genes))
tot_overlap_intervals_5mb <- unique(rbind(overlap_con_intervals_filt_genes, overlap_gf_intervals_filt_genes))

write(tot_overlap_intervals$gene, file="intervals/tot_overlap_intervals_SW.txt")
write(tot_overlap_intervals_5mb$gene, file="intervals/tot_overlap_intervals_filt_5mb_SW.txt")
write(overlap_con_intervals_filt_genes$gene, file="intervals/overlap_intervals_con_filt_5mb_genes_SW.txt")
write(overlap_gf_intervals_filt_genes$gene, file="intervals/overlap_intervals_gf_filt_5mb_genes_SW.txt")
write(overlap_gf_intervals_genes$gene, file="intervals/overlap_intervals_gf_genes_SW.txt")
write(overlap_con_intervals_genes$gene, file="intervals/overlap_intervals_con_genes_SW.txt")


####### snps ######
# file from all_distinct_genes_DNARNA.R
all_genes_snps <- read.delim("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/all_genes_snps_SW.txt", header = F, col.names = "gene")
all_genes_snps$gene <- splus2R::upperCase(all_genes_snps$gene)
overlap_gf_all_genes_snps<-data.frame()
for (gene in genes_gf$gene){
  xx<-str_which(all_genes_snps$gene, gene)
  if (length(xx)>0){
    result<-input_SW[str_which(input_SW$list_total_genes,regex(gene,ignore_case = T)),]
    result$gene <- gene
    overlap_gf_all_genes_snps<- rbind(overlap_gf_all_genes_snps,result)
  }
}
write_excel_csv(overlap_gf_all_genes_snps, "snps/overlap_gf_all_genes_snps.csv")

overlap_gf_all_genes_snps_genes<- overlap_gf_all_genes_snps %>% 
  distinct(gene)

overlap_con_all_genes_snps<-data.frame()
for (gene in genes_con$gene){
  xx<-str_which(all_genes_snps$gene, gene)
  if (length(xx)>0){
    result<-input_SW[str_which(input_SW$list_total_genes,regex(gene,ignore_case = T)),]
    result$gene <- gene
    overlap_con_all_genes_snps<- rbind(overlap_con_all_genes_snps,result)
  }
}
write_csv2(overlap_con_all_genes_snps, "snps/overlap_con_all_genes_snps.csv")

overlap_con_all_genes_snps_genes<- overlap_con_all_genes_snps %>% 
  distinct(gene)

overlap_all_genes_snps_total <- unique(rbind(overlap_con_all_genes_snps_genes, overlap_gf_all_genes_snps_genes))

write(overlap_all_genes_snps_total$gene, file="snps/overlap_all_genes_snps_total.txt")
write(overlap_con_all_genes_snps$gene, file="snps/overlap_con_all_genes_snps.txt")
write(overlap_gf_all_genes_snps$gene, file="snps/overlap_gf_all_genes_snps.txt")

