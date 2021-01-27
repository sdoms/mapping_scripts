setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/genes_overlap/")
library(tidyverse)
library(stringr)
library(readxl)
all_genes <- read.delim("../Results/Bacterial traits/Genes/all_genes_intervals/all_genes_DNA_RNA.txt", header=F)
all_genes_snps <- drop_na(read.delim("../Results/Bacterial traits/Genes/all_genes_snps/all_genes_snps.txt", header=F))
all_genes$x <- splus2R::upperCase(all_genes$V1)
all_genes_snps$x <- splus2R::upperCase(all_genes_snps$V1)

load("../Results/Bacterial traits/genes_DNA_RNA.Rdata")
# load gene list of other studies
turpin<- read_xlsx("Other_studies.xlsx", sheet="Turpin_2016", skip=1)
turpin$...1<- NULL
turpin_genes <- unlist(str_split(turpin$`Gene(s) in the region`, "/"))
turpin_genes <- unlist(str_split(turpin_genes, ","))


dir.create("other_studies")
dir.create("other_studies/turpin")
overlap<-data.frame()
for (gene in turpin_genes){
  xx<-str_which(all_genes$x, gene)
  if (length(xx)>0){
    
    result<-genes_DNA_RNA[str_which(genes_DNA_RNA$list_total_genes,regex(gene,ignore_case = T)),]
    write_delim(result,file=paste0("other_studies/turpin/Gene_",gene,"_results.txt"))
    result$gene <- gene
    overlap<- rbind(overlap,result)
  }
}

overlap_filt <- overlap %>% 
  filter(as.numeric(length)<5)

gggg<- overlap_filt %>% 
  distinct(gene)
write.table(gggg, "other_studies/turpin/overlap_genes.txt")

# check with genes closest to significant snps
dir.create("other_studies/turpin/snps")
overlap_snps<-data.frame()
for (gene in turpin_genes){
  xx<-str_which(all_genes_snps$x, gene)
  if (length(xx)>0){
    
    result<-genes_DNA_RNA[str_which(genes_DNA_RNA$list_total_genes,regex(gene,ignore_case = T)),]
    write_delim(result,file=paste0("other_studies/turpin/snps/Gene_",gene,"_results.txt"))
    result$gene <- gene
    overlap_snps<- rbind(overlap_snps,result)
  }
}

overlap_snps_filt <- overlap_snps %>% 
  filter(as.numeric(length)<5)

gggg_snps<- overlap_snps_filt %>% 
  distinct(gene)
write.table(gggg_snps, "other_studies/turpin/snps/overlap_genes_snps.txt")


#### Wang_2015 ####
# load gene list of other studies
wang<- read_xlsx("Other_studies.xlsx", sheet="Wang_2015")
load("../../Genotype_data/GM_snps.Rdata")
gm

dir.create("other_studies/wang")

##### Leamy 2014 ####
leamy<- read_xlsx("Other_studies.xlsx", sheet="Leamy_2014")

dir.create("other_studies/leamy")
overlap<-data.frame()
line_nrs<- c()

for (line in 1:nrow(leamy)){
  interval <- str_split(leamy[line,'CI (Mb)'], "-")
  start_pos<- as.numeric(as.character(interval[[1]][1]))
  end_pos<- as.numeric(as.character(interval[[1]][2]))
  chr <- as.character(leamy[line, "Ch"])
  ss<-genes_DNA_RNA[genes_DNA_RNA$chr==chr,]%>% 
    drop_na(name)
  
  result <- ss %>% 
    filter(position_peak_snp<end_pos & position_peak_snp>start_pos)

    overlap<- rbind(overlap,result)
    if (dim(result)>0) {
      line_nrs<- append(line_nrs, line)
    }
}

overlap_hits<- leamy[line_nrs,]
overlap_filt <- overlap %>% 
  filter(as.numeric(length)<5)

gggg<- overlap_filt %>% 
  distinct(peak_snp)
write.table(overlap, "other_studies/leamy/overlap_regions.txt")
write.table(overlap_filt, "other_studies/leamy/overlap_regions_filtered.txt")
write.table(overlap_hits, "other_studies/leamy/overlap_hits.txt")


#### Org 2015 ####
org<- read_xlsx("Other_studies.xlsx", sheet="Org_2015")
org<- org[3:9, -c(5:9)]
dir.create("other_studies/org")
overlap<-data.frame()
line_nrs<- c()

for (line in 1:nrow(org)){
  interval <- str_split(org[line,'LD block (Mb)'], "-")
  start_pos<- as.numeric(as.character(interval[[1]][1]))
  end_pos<- as.numeric(as.character(interval[[1]][2]))
  chr <- as.character(org[line, "Chromosome"])
  ss<-genes_DNA_RNA[genes_DNA_RNA$chr==chr,]%>% 
    drop_na(name)
  
  result <- ss %>% 
    filter(position_peak_snp<end_pos & position_peak_snp>start_pos)
  
  overlap<- rbind(overlap,result)
  if (dim(result)>0) {
    line_nrs<- append(line_nrs, line)
  }
}

overlap_hits<- org[line_nrs,]
overlap_filt <- overlap %>% 
  filter(as.numeric(length)<5)

gggg<- overlap_filt %>% 
  distinct(peak_snp)
write.table(overlap, "other_studies/org/overlap_regions.txt")
write.table(overlap_filt, "other_studies/org/overlap_regions_filtered.txt")
write.table(overlap_hits, "other_studies/org/overlap_hits.txt")

###################################################

bedfile<- drop_na(genes_DNA_RNA[, c("chr", "start", "stop", "length")])
bedfile<- apply(bedfile, 2, function(x) as.numeric(as.character(x)))
for (i in 1:nrow(bedfile)){
  if (bedfile[i, "length"]==0){
    bedfile[i, "start"]<-bedfile[i, "start"]-1.45e6
    bedfile[i, "stop"]<-bedfile[i, "stop"] +1.45e6
  }
}
bedfile<- as.data.frame(bedfile)
bedfile<- bedfile[order(bedfile$chr,bedfile$start),]
bedfile$chr<- paste0("chr", bedfile$chr)
bedfile<- unique(bedfile[,-4])
write_excel_csv(bedfile, file="my_intervals_DNA_RNA.csv")

# centromeres mouse
gaps<- read.table("gap.txt")
centromeres<- gaps[which(gaps$V8=="centromere"),-1]
centromeres<- centromeres[gtools::mixedorder(centromeres$V2),]
write.table(centromeres[1:3], file="centromeres.bed", row.names = F, col.names = F, sep = "\t", quote = F)

# studies ##
wang_bed<- wang[, c("chromsome", "start", "end")]
wang_bed<- wang_bed[order(wang_bed$chromsome, wang_bed$start),]
wang_bed$start <- wang_bed$start*1e6
wang_bed$end <- wang_bed$end*1e6
wang_bed$chromsome <- paste0("chr", wang_bed$chromsome)
write.table(wang_bed, file="wang.bed", row.names = F, col.names=F, sep = "\t", quote = F)


