setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/")

library(biomaRt)
library(readxl)
library(tidyverse)
mills_genes <- read_excel("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/overexpressed_protein_GF_con.xlsx", sheet="ALL_DF", col_names = F)
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = 'www') # can take long



attributes.1 = c("chromosome_name","start_position", "end_position","strand", "mgi_symbol")
results.1 = getBM(attributes = attributes.1, filters = c("mgi_symbol"), values = mills_genes$...1, mart = gene.ensembl)

bed_results <- results.1 %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  dplyr::select(chr, start_position, end_position) 
write.table(bed_results, "mills_genes_bed.txt")
  

load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/mills/genes_DNA_RNA_intervals_SW.Rdata")
input_SW$length<- (input_SW$stop.LD.pos-input_SW$start.LD.pos)/1e6
input_SW<- transform(input_SW, length_num=as.numeric(length))
input_SW_unique_int <- input_SW %>% 
  group_by(chr) %>% 
  distinct(start.LD.pos, stop.LD.pos) %>% 
  mutate(chr=paste0("chr",chr))

write.table(input_SW_unique_int, file="all_intervals_our_study.txt", row.names = F, col.names = F)

system("~/poverlap/poverlap2.py poverlap --a mills_genes_bed.txt --b all_intervals_our_study.txt -g ../mouse.mm10.genome --n 9999")
system("bedtools intersect -wo -a mills_genes_bed.txt -b all_intervals_our_study.txt > bed_intersect_out_mills.txt")

mills_intersect <- read_tsv("bed_intersect_out_mills.txt", col_names = F)

genes_overlap <- results.1 %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  inner_join(mills_intersect, by=c("chr"="X1","start_position"="X2" ))

genes_mills <- unique(genes_overlap$mgi_symbol)
## 190 genes -> fine 



