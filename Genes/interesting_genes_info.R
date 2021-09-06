setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/")
library(readxl)
library(stringr)
library(readr)
library(tidyverse)
library(biomaRt)
#load("all_genes_snps/enrichment/All_info_markers_genes_intervals_enrichment.Rdata")
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version=102) 
genes <- read_excel("Interesting_genes.xlsx", sheet = "all_genes", col_names = "Genes")

genes.col <- stringr::str_to_title(genes$Genes)
dir.create("interesting_genes")
load("freq_overlap.Rdata")
freq_overlaps_other_studies <- output
hub_genes_neighbours_snps <- c("Kcnj3", "Fzd4", "Plcb1", "Adcy7", "Wnt11", "Adcy2", "Pth", "Calca", "Adcyap1", "Grik3", "Itpr1", "Pik3r1", "Galr1", "Grm8", "Cxcl12", "Oxgr1","Htr1f", "Grm7",
                               "C5ar2", "Serpina3b", "Dync1h1", "Nhlrc3", "Cfb", "Igfbp5", "C4b", "Fbn1", "Mfge8", "Spp2", "Itih2")
hub_genes_snps <- c("Gng12", "Nmur2", "Mchr1", "Xcr1", "Chrm3", "Prok2", "Tacr3", "Ptgfr", "Oxtr", "C3")

hub_genes_overlap <- c("Gstm5", "Gsta1", "Gsta4", "Cyp2c55", "Cyp2c65", "Cyp2b10", "Hspa1b", "Hspa1a1", "Dnajb11")
hub_genes_neighbours_overlap <- c("Alox15", "Cyp2d26", "Ephx2", "Gclc", "Alpi", "Timm50", "Tomm22", "Phb", "Hspb1", "Hspa1l", "Ostc", "Sparcl1", "Lamb1", "Mia3", "Ero1l", "Actr10", "Lyz1", "Hexb")

genes_list_SW_in_overlaps <-read.table("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/genes_intervals_in_other_studies-SW.txt")
all_genes<- genes_list_SW_in_overlaps$V1
marker_genes<- read.csv("../Shared/all_markers_RNA_DNA-with-genes_SW.csv", sep=",")

marker_genes$total_genes <- paste(marker_genes$all_genes, marker_genes$all_following_genes, 
                                  marker_genes$all_preceding_genes, sep=", ")

all_marker_genes <- read.table( "all_marker_genes.txt")

all_mills_genes <- read_excel("overexpressed_protein_GF_con.xlsx", "ALL_DF", col_names = F)
genes_mills<- all_mills_genes$...1

genes.col <- unique(str_to_title(unlist(c(genes.col, hub_genes_neighbours_overlap, hub_genes_neighbours_snps, hub_genes_overlap, hub_genes_snps ))))

out<- data.frame(matrix(ncol=16))
colnames(out)<- c("gene","snps", "mills", "other studies","hub genes network snps", "first neighbours hubs snps network", "hub genes network mills", "first neighbours hub genes network mills", "frequency overlaps with other studies", "lowest p type", "lowest p val","min_interval", "number of taxa intervals","taxa interval", "number of taxa marker","taxa marker")
all_dna_rna_SW$length<- all_dna_rna_SW$stop.LD.pos-all_dna_rna_SW$start.LD.pos
for (gene in genes.col){
  xx<-str_which(str_to_title(all_dna_rna_SW$list_total_genes), regex(paste0("\\b",gene,"\\b(?!-)"), ignore_case = T))
  yy<- str_which(str_to_title(all_marker_genes), regex(paste0("\\b",gene,"\\b(?!-)"), ignore_case = T))
  
  if (length(xx)>0){
    result_i<-all_dna_rna_SW[str_which(all_dna_rna_SW$list_total_genes,regex(paste0("\\b",gene,"\\b(?!-)"),ignore_case = T)),]
    write_delim(result_i,file=paste0("interesting_genes/",gene,"_intervals_results.csv"), delim=";")
    min_int <- min(result_i$length, na.rm = T)
    lowest_p <- min(result_i$peak.P.type, na.rm=T)
    lowest_p_type <- result_i$P.type[which.min(result_i$peak.P.type)]
    taxa_i<- paste(result_i$dna.rna,":", result_i$trait, collapse = " | ")
    count_taxa_I <- length(unique(result_i$trait))
  } else {
    min_int <- NA
    taxa_i <- NA
    count_taxa_I <- 0
  }
  if (length(yy)>0){
    result_m<-marker_genes[str_which(marker_genes$total_genes,regex(paste0("\\b",gene, "\\b(?!-)"),ignore_case = T)),]
    write_delim(result_m,file=paste0("interesting_genes/",gene,"_marker_results.csv"), delim=";")
    taxa_m<- paste(result_m$all_taxa, collapse = " | ")
    ll <- trimws(unlist(str_split(result_m$all_taxa, "[:|]")))
    count_taxa_m <- length(unique(ll[which(! (ll %in% c("DNA", "RNA")))]))
    
    result_i2_all<- data.frame()
    for (marker in result_m$marker){
      result_i2<-all_dna_rna_SW[str_which(all_dna_rna_SW$markers,regex(marker,ignore_case = T)),]
      result_i2_all <- rbind(result_i2_all, result_i2)
    } 
    
    result_i2_all <- result_i2_all %>% 
      distinct()
    
    lowest_p2 <- min(result_i2_all$peak.P.type, na.rm=T)
    lowest_p_type2 <- result_i2_all$P.type[which.min(result_i2_all$peak.P.type)]
    if (is_empty(xx)){
      
      result_i<-result_i2_all
      min_int<- min(result_i$length, na.rm = T)
      
      taxa_i<- paste(result_i$dna.rna, result_i$trait, collapse = " | ")
      count_taxa_I <- length(unique(result_i$trait))
      write_delim(result_i,file=paste0("interesting_genes/",gene,"_intervalsM_results.csv"), delim=";")
      
    }
    
    lowest_p<- min(lowest_p, lowest_p2, na.rm = T)
    lowest_p_type <- c(lowest_p_type, lowest_p_type2)[which.min(c(lowest_p,lowest_p2))]
  } else {
    taxa_m<-NA
    count_taxa_m <-0
  }
  if (str_to_upper(gene) %in% genes_mills){
    mills<-"yes"
  } else {
    mills <- "no"
  }
  if (gene %in% genes_list_SW_in_overlaps$V1){
    other_studies <- "yes"
    freq_OV_OS <- tryCatch(max(freq_overlaps_other_studies[str_which(freq_overlaps_other_studies$all.genes,regex(paste0("\\b",gene, "\\b(?!-)"),ignore_case = T)),"Freq"]), error=function(x) NA)
  } else {
    other_studies <- "no"
    freq_OV_OS <- 0
  }
  if (str_to_upper(gene) %in% all_marker_genes$V1){
    snps <- "yes"
  } else {
    snps <- "no"
  }
if (gene %in% hub_genes_snps){
  in_hub_snps <- "yes"
} else {
  in_hub_snps <- "no"
}
  if (gene %in% hub_genes_overlap){
    in_hub_overlap <- "yes"
  } else {
    in_hub_overlap <- "no"
  }
  
  if (gene %in% hub_genes_neighbours_snps){
    in_hub_neighbours_snps <- "yes"
  } else {
    in_hub_neighbours_snps <- "no"
  }
  
  if (gene %in% hub_genes_neighbours_overlap){
    in_hub_neighbours_overlap <- "yes"
  } else {
    in_hub_neighbours_overlap <- "no"
  }
  
  output<- c(gene,snps, mills,other_studies, in_hub_snps, in_hub_neighbours_snps, in_hub_overlap, in_hub_neighbours_overlap,freq_OV_OS,lowest_p_type, lowest_p,min_int,count_taxa_I,taxa_i, count_taxa_m, taxa_m )
  out<- rbind(out, output)
}

out <- out[-1,]



out.bm.genes.region <- getBM(
  attributes = c('chromosome_name','start_position','end_position','ensembl_gene_id','external_gene_name',  "description"), 
  filters = c('external_gene_name'), 
  values = genes.col, 
  mart = gene.ensembl)

out_info <- out %>% 
  inner_join(out.bm.genes.region, by=c("gene"="external_gene_name"))
write_delim(out_info, "interesting_genes_infoR.csv", delim=";")
DT::datatable(out_info)




### mills genes output file for network ####3
out <- data.frame()
for (gene in genes_mills) {
  xx<-str_which(all_dna_rna_SW$list_total_genes, regex(paste0("\\b",gene,"\\b(?!-)"), ignore_case = T))
  if (length(xx)>0){
    if (str_to_title(gene) %in% genes_list_SW_in_overlaps$V1){
      other_studies <- "yes"
      freq_OV_OS <- tryCatch(max(freq_overlaps_other_studies[str_which(freq_overlaps_other_studies$all.genes,regex(paste0("\\b",gene, "\\b(?!-)"),ignore_case = T)),"Freq"]), error=function(x) NA)
    } else {
      other_studies <- "no"
      freq_OV_OS <- 0
    }
    if (gene %in% all_marker_genes$V1){
      snps <- "yes"
    } else {
      snps <- "no"
    }
    if (snps=="yes"&other_studies=="yes"){
      both<- "both"
    } else if (snps=="yes"&other_studies=="no"){
      both <- "snps"
    } else if (snps=="no"&other_studies=="yes"){
      both <- "overlap"
    } else{
      both <- "mills"
    }
    result_i<-all_dna_rna_SW[str_which(all_dna_rna_SW$list_total_genes,regex(paste0("\\b",gene,"\\b(?!-)"),ignore_case = T)),]
    if (nrow(result_i)<1){
      both<- "no"
      count_taxa <-NA
    } else{
      count_taxa<- length(unique(result_i$trait))
    }
    
    DNA_orRNA<- unique(result_i$dna.rna)
    if (length(DNA_orRNA)>1){
      DNA_orRNA<- "both"
    }
    output<- c(gene, both, DNA_orRNA, count_taxa)
    out<- rbind(out, output)
  }
 
}
colnames(out)<- c("gene", "both", "DNA_or_RNA", "count_taxa")

reactable::reactable(out)
write.csv(out, file="genes_overlap/SW/info_mills_genes_for_network.csv")

for (i in 1:nrow(marker_genes)){
  ll<-trimws(unlist(str_split(marker_genes$all_taxa[i], "[:|]")))
  marker_genes$number_of_taxa[i]<-length(unique(ll[which(! (ll %in% c("DNA", "RNA")))]))
}
top20<- slice_max(marker_genes, number_of_taxa, n=20)
genes_top20 <- unique(unlist(str_split(top20$all_genes, " | ")))
 genes_top20 <- genes_top20[-c(2:4)]                     

 top10_overlaps<- slice_max(freq_overlaps_other_studies, Freq, n=10) 
 genes_top10overlaps <- unique(unlist(str_split(top10_overlaps$protein.coding.genes, ",")))
 