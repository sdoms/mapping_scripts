setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/")

library(biomaRt)
library(readxl)
library(tidyverse)
ibd_genes <- read_excel("../../Inflammtory Bowel Disease.xlsx", sheet=1, col_names = T)
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = 'www') # can take long

snps_genes<- read_delim("../../all_marker_genes.txt", col_names =F, delim=" ")


# Basic function to convert mouse to human gene names


convertToMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
humanx<- convertToMouseGeneList(ibd_genes$`HGNC gene symbol`)

convertToHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


attributes.1 = c("chromosome_name","start_position", "end_position","strand", "mgi_symbol")
results.1 = getBM(attributes = attributes.1, filters = c("mgi_symbol"), values = humanx, mart = gene.ensembl)

bed_results <- results.1 %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  dplyr::select(chr, start_position, end_position) 
write.table(bed_results, "ibd_genes_bed.txt")


load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_intervals/SW/genes_DNA_RNA_intervals_SW.Rdata")
input_SW$length<- (input_SW$stop.LD.pos-input_SW$start.LD.pos)/1e6
input_SW<- transform(input_SW, length_num=as.numeric(length))
input_SW_unique_int <- input_SW %>% 
  group_by(chr) %>% 
  distinct(start.LD.pos, stop.LD.pos) %>% 
  mutate(chr=paste0("chr",chr))

write.table(input_SW_unique_int, file="all_intervals_our_study.txt", row.names = F, col.names = F)

system("~/poverlap/poverlap2.py poverlap --a ibd_genes_bed.txt --b all_intervals_our_study.txt -g ../mouse.mm10.genome --n 9999")
# {"'wc -l'": {"observed": 328.0, "shuffle_cmd": "bedtools intersect -wo -a ibd_genes_bed.txt -b <(bedtools shuffle   -i all_intervals_our_study.txt -g ../mouse.mm10.genome )", 
# "metric": "'wc -l'", "simulated mean metric": 288.98399839984, "simulated_p": 0.1513}}

system("bedtools intersect -wo -a ibd_genes_bed.txt -b all_intervals_our_study.txt > bed_intersect_out_ibd.txt")

ibd_intersect <- read_tsv("bed_intersect_out_ibd.txt", col_names = F)

genes_overlap <- results.1 %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  inner_join(ibd_intersect, by=c("chr"="X1","start_position"="X2" ))

genes_ibd <- unique(genes_overlap$mgi_symbol)
## 169 genes in our intervals
# add info from ibd_genes
human_ibd <- convertToHumanGeneList(genes_ibd)
genes_ibd_info <- ibd_genes[which(ibd_genes$`HGNC gene symbol`%in%human_ibd),]
write_delim(genes_ibd_info, file="human_IBD_genes_in_interval_genes.csv", delim=";")

#### snps genes ####

results.snps = getBM(attributes = attributes.1, filters = c("mgi_symbol"), values = snps_genes, mart = gene.ensembl)

bed_results.snps <- results.snps %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  dplyr::select(chr, start_position, end_position) 
write.table(bed_results, "snps_genes_bed.txt")

system("~/poverlap/poverlap2.py poverlap --a ibd_genes_bed.txt --b snps_genes_bed.txt -g ../mouse.mm10.genome --n 9999")
#{"'wc -l'": {"observed": 327.0, "shuffle_cmd": "bedtools intersect -wo -a ibd_genes_bed.txt -b <(bedtools shuffle   -i snps_genes_bed.txt -g ../mouse.mm10.genome )", 
# "metric": "'wc -l'", "simulated mean metric": 2.705770577057706, "simulated_p": 0.0001}}

system("bedtools intersect -wo -a ibd_genes_bed.txt -b snps_genes_bed.txt > bed_intersect_out_ibd_with_snp_genes.txt")

ibd_intersect_snps <- read_tsv("bed_intersect_out_ibd_with_snp_genes.txt", col_names = F)

genes_overlap_snps_with_ibd <- results.snps %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  inner_join(ibd_intersect_snps, by=c("chr"="X1","start_position"="X2" ))

genes_ibd_with_snps <- unique(genes_overlap_snps_with_ibd$mgi_symbol)
# add info from ibd_genes
human_ibd <- convertToHumanGeneList(genes_ibd_with_snps)
genes_ibd_with_snps_info <- ibd_genes[which(ibd_genes$`HGNC gene symbol`%in%human_ibd),]
write_delim(genes_ibd_with_snps_info, file="human_IBD_genes_in_closest_genes.csv", delim=";")


sum(toupper(humanx) %in% snps_genes$X1)

fisher.test(matrix(data=c(14,912,305,59037), nrow=2))

# get info on taxa

all_markers_genes <- read.delim("../../all_genes_snps/all_markers_RNA_DNA-with-genes_SW.csv", sep=";")
info_marker_genes_ibd <- all_markers_genes[which(all_markers_genes$all_genes%in% genes_ibd_with_snps |all_markers_genes$all_preceding_genes%in% genes_ibd_with_snps | all_markers_genes$all_following_genes%in% genes_ibd_with_snps ),]
write_delim(info_marker_genes_ibd, file="ibd_genes_snps_with_taxa.csv", delim=";")
