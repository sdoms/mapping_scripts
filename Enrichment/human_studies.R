setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/")

library(biomaRt)
library(readxl)
library(tidyverse)
human_studies <- read_excel("../other_studies/human_studies.xlsx", sheet="Genes", col_names = T)

all_human_genes <- c()
for (i in human_studies$Genes){
  all_genes <- unlist(strsplit(i, ","))
  all_human_genes <- c(all_human_genes, trimws(all_genes))
}
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = 'www') # can take long

convertToMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = trimws(x) , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
humanx<- convertToMouseGeneList(trimws(all_human_genes))

attributes.1 = c("chromosome_name","start_position", "end_position","strand", "mgi_symbol")
results.1 = getBM(attributes = attributes.1, filters = c("mgi_symbol"), values = humanx, mart = gene.ensembl)

bed_results <- results.1 %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  dplyr::select(chr, start_position, end_position) 
write.table(bed_results, "human_studies_genes_bed.txt")


load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/mills/genes_DNA_RNA_intervals_SW.Rdata")
input_SW$length<- (input_SW$stop.LD.pos-input_SW$start.LD.pos)/1e6
input_SW<- transform(input_SW, length_num=as.numeric(length))
input_SW_unique_int <- input_SW %>% 
  group_by(chr) %>% 
  distinct(start.LD.pos, stop.LD.pos) %>% 
  mutate(chr=paste0("chr",chr))

write.table(input_SW_unique_int, file="all_intervals_our_study.txt", row.names = F, col.names = F)

system("~/poverlap/poverlap2.py poverlap --a human_studies_genes_bed.txt --b all_intervals_our_study.txt -g ../mouse.mm10.genome --n 9999")
#{"'wc -l'": {"observed": 403.0, "shuffle_cmd": "bedtools intersect -wo -a human_studies_genes_bed.txt -b <(bedtools shuffle   -i all_intervals_our_study.txt -g ../mouse.mm10.genome )", 
  #"metric": "'wc -l'", "simulated mean metric": 367.93099309930994, "simulated_p": 0.156}}

system("bedtools intersect -wo -a human_studies_genes_bed.txt -b all_intervals_our_study.txt > bed_intersect_out_human_genes.txt")

human_intersect <- read_tsv("bed_intersect_out_human_genes.txt", col_names = F)

genes_overlap <- results.1 %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  inner_join(human_intersect, by=c("chr"="X1","start_position"="X2" ))

genes_human <- unique(genes_overlap$mgi_symbol)
write_delim(as.data.frame(genes_human), file="human_genes_in_interval_genes.csv", delim=";")



### snps genes


#### snps genes ####

results.snps = getBM(attributes = attributes.1, filters = c("mgi_symbol"), values = snps_genes, mart = gene.ensembl)

bed_results.snps <- results.snps %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  dplyr::select(chr, start_position, end_position) 
write.table(bed_results.snps, "snps_genes_bed.txt")

system("~/poverlap/poverlap2.py poverlap --a human_studies_genes_bed.txt --b snps_genes_bed.txt -g ../mouse.mm10.genome --n 9999")
#{"'wc -l'": {"observed": 37.0, "shuffle_cmd": "bedtools intersect -wo -a human_studies_genes_bed.txt -b <(bedtools shuffle   -i snps_genes_bed.txt -g ../mouse.mm10.genome )",
#"metric": "'wc -l'", "simulated mean metric": 39.42414241424142, "simulated_p": 0.6181}}

system("bedtools intersect -wo -a human_studies_genes_bed.txt -b snps_genes_bed.txt > bed_intersect_out_human_genes_with_snp_genes.txt")

human_intersect_snps <- read_tsv("bed_intersect_out_human_genes_with_snp_genes.txt", col_names = F)

genes_overlap_snps_with_human <- results.1 %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  inner_join(human_intersect_snps, by=c("chr"="X1","start_position"="X2" ))

genes_human_with_snps <- unique(genes_overlap_snps_with_human$mgi_symbol)

write_delim(as.data.frame(genes_human_with_snps), file="human_genes_in_closest_genes.csv", delim=";")


inner_join(as.data.frame(toupper(humanx)), snps_genes, by=c("toupper(humanx)"="X1"))

fisher.test(matrix(data=c(32,894,261,59099), nrow=2))


# malte's genes 
maltesgenes <- read_excel("../other_studies/human_studies.xlsx", sheet="ruehlemann_2021", col_names = T)

all_maltes_genes <- c()
for (i in maltesgenes$`Genesinlocus(±100kb)`){
  all_genes <- unlist(strsplit(i, ","))
  all_maltes_genes <- c(all_maltes_genes, trimws(all_genes))
}

humanx_malte<- convertToMouseGeneList(trimws(all_maltes_genes))

attributes.1 = c("chromosome_name","start_position", "end_position","strand", "mgi_symbol")
results.malte = getBM(attributes = attributes.1, filters = c("mgi_symbol"), values = humanx_malte, mart = gene.ensembl)

bed_results <- results.malte %>% 
  mutate(chr=paste0("chr", chromosome_name)) %>% 
  dplyr::select(chr, start_position, end_position) 
write.table(bed_results, "human_studies_malte_genes_bed.txt")
#inner_join(as.data.frame(toupper(humanx_malte)), snps_genes, by=c("toupper(humanx_malte)"="X1"))

system("~/poverlap/poverlap2.py poverlap --a human_studies_malte_genes_bed.txt --b snps_genes_bed.txt -g ../mouse.mm10.genome --n 9999")
system("~/poverlap/poverlap2.py poverlap --a human_studies_malte_bed.txt --b snps_genes_bed.txt -g ../mouse.mm10.genome --n 9999")
