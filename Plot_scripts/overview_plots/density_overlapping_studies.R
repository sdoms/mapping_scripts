setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/")

library(tidyverse)
library(biomaRt)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
mm10 <- read.table("../mouse.mm10.genome")
colnames(mm10)<- c("chr", "chr_len")


overlaps <- read_delim("bed_intersect_out.txt",delim = "\t", col_names = F)
colnames(overlaps)<- c("chrA", "startA", "endA","chrB", "startB", "endB", "width")
overlap_intervals <- overlaps %>% 
  rowwise() %>% 
  mutate(start=max(startA, startB), end=min(endA, endB)) %>% 
  mutate(length=end-start)

overlap_intervals_uniq <- overlap_intervals %>% 
  distinct(chrA, start, end)


goodChrOrder <- c(1:19)

output <- data.frame()
for (chr in goodChrOrder){
  df.tmp <- overlap_intervals[overlap_intervals$chrA==paste0("chr",chr),]
  dd<- c()
  for (x in 1:nrow(df.tmp)){
    dd<- c(dd,seq(df.tmp$start[x],df.tmp$end[x], 1e6))
  }
  res<-as.data.frame(table(cut_interval(dd,length=1e6)))
  
  #res <- as.data.frame(table(cut_width(df.tmp$start/1e6,1, boundary = 0)))

  res$start <- as.numeric(as.character(substr(res$Var1,2,regexpr(",",res$Var1)-1)))
  res$end <- res$start +1e6
  res$chr <- chr
  output <- rbind(output,res)
}

# output <- data.frame()
# for (chr in goodChrOrder){
#   df.tmp <- overlaps[overlaps$chr==paste0("chr",chr),]
#   res <- as.data.frame(table(cut_width(df.tmp$start/1e6,1, boundary = 0)))
#   res2<-as.data.frame(table(cut_width(df.tmp$end/1e6,1, boundary = 0)))
#   res$start <- as.numeric(as.character(substr(res$Var1,2,regexpr(",",res$Var1)-1)))
#   res$end <- res$start +1
#   res$chr <- chr
#   output <- rbind(output,res)
# }
mean_freq <- mean(output$Freq, na.rm=T)

output <- output %>% filter(Freq != 0)

chrSizes <- mm10[1:19,]
chrSizes$chr<- 1:19

output$start<- as.numeric(as.character(output$start))/1e6
output$end<- as.numeric(as.character(output$end))/1e6
library(viridis)
snp_plot <- ggplot()+
  geom_segment(data = chrSizes,
               aes(x = chr, xend = chr, y = 0, yend = chr_len/1e6),
               lineend = "round", color = "#DEE1E3", size = 5) + 
  
  geom_segment(data=output, 
               aes(x=chr, xend=chr, y=start, yend=end, color=Freq), size=5)+
  coord_flip() + scale_color_viridis_c(option="magma", direction = -1, limits=c(0,10)) + 
  labs(y="Position (Mb)", x= "Chromosome", color="Times found in other studies") +theme_test()
snp_plot + scale_x_discrete(limits=factor(goodChrOrder))
ggsave("overlapping_studies_intervals_SW.pdf")
write_delim(output, "freq_overlaps.csv", delim=";")

### get genes in intervals #### 
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www") # can take long
top20<- output %>% slice_max(Freq,n=20)
for (i in 1:nrow(top20)){
  chr <- top20$chr[i]
  start <- top20$start[i]*1e6
  stop <- top20$end[i] *1e6
  out.bm.genes.region <- getBM(
    attributes = c('external_gene_name'), 
    filters = c('chromosome_name','start','end', 'biotype'), 
    values = list(chr,start, stop, "protein_coding"), 
    mart = gene.ensembl)
  out.bm.genes.region.all <- getBM(
    attributes = c('external_gene_name', 'gene_biotype'), 
    filters = c('chromosome_name','start','end'), 
    values = list(as.character(chr),start, stop), 
    mart = gene.ensembl)
  if (nrow(out.bm.genes.region>0)){
    top20$protein.coding.genes[i]<- paste(out.bm.genes.region$external_gene_name, collapse=',')
  }
  if (nrow(out.bm.genes.region.all)>0){
    top20$all.genes[i] <- paste(out.bm.genes.region.all$external_gene_name, collapse = ',')
    top20$all.genes.biotype[i]<- paste(out.bm.genes.region.all$gene_biotype, collapse = ',')
  }
  
}

#### all ####
for (i in 1:nrow(output)){
  chr <- output$chr[i]
  start <- output$start[i]*1e6
  stop <- output$end[i] *1e6
  out.bm.genes.region <- getBM(
    attributes = c('external_gene_name'), 
    filters = c('chromosome_name','start','end', 'biotype'), 
    values = list(chr,start, stop, "protein_coding"), 
    mart = gene.ensembl)
  out.bm.genes.region.all <- getBM(
    attributes = c('external_gene_name', 'gene_biotype'), 
    filters = c('chromosome_name','start','end'), 
    values = list(as.character(chr),start, stop), 
    mart = gene.ensembl)
  if (nrow(out.bm.genes.region>0)){
    output$protein.coding.genes[i]<- paste(out.bm.genes.region$external_gene_name, collapse=',')
  }
  if (nrow(out.bm.genes.region.all)>0){
    output$all.genes[i] <- paste(out.bm.genes.region.all$external_gene_name, collapse = ',')
    output$all.genes.biotype[i]<- paste(out.bm.genes.region.all$gene_biotype, collapse = ',')
  }
  
}
save(output, file="freq_overlap.Rdata")
