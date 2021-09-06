setwd("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_intervals/SW/")
load("../../../../../Cleaning_snps/clean_snps.Rdata")

#load("genes_DNA_RNA_intervals_SW.Rdata")
load(file="/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_02032021_GW.Rdata")
load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_03032021_SW.Rdata")

# system("~/poverlap/poverlap2.py poverlap --a all_intervals_mouse_studies_10mb.txt --b my_intervals_DNA_RNA_10mb.txt -g mouse.mm10.genome --n 9999")
# system("bedtools intersect -wo -a all_intervals_mouse_studies_10mb.txt -b my_intervals_DNA_RNA_10mb.txt >bed_intersect_out.txt")

library(tidyverse)
library(biomaRt)
mm10 <- read.table("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/mouse.mm10.genome")
colnames(mm10)<- c("chr", "chr_len")

# put the chromosomes in the good order: chr1, chr2, chr22, chrX
goodChrOrder <- c(1:20)


# results file 

# input_file <- read_csv2("all_sig_RNA.csv")
# input_file <- input_file %>%  drop_na(trait)
# input_file_RNA$DNAorRNA<-"RNA"
# input_file_DNA$DNAorRNA<-"DNA"
# 
# input_file_DNA<- read_csv2("../../DNA/sig.summaries/all_sig_DNA.csv")
# input_file_RNA <- input_file
# input_file<- rbind(input_file_DNA,input_file_RNA)

input_file <- all_dna_rna_SW
input_DNA<- all_dna_rna_SW %>% 
  filter(dna.rna=="DNA")
input_RNA<- all_dna_rna_SW %>% 
  filter(dna.rna=="RNA")
med_DNA<- median(input_DNA$length[input_DNA$length>0], na.rm=T)
med_RNA<- median(input_RNA$length[input_RNA$length>0], na.rm=T)
med_tot <-median(input_file$length[input_file$length>0], na.rm=T)
#start plot
output <- data.frame()
for (chr in goodChrOrder){
  df.tmp <- input_file[input_file$chr.num==as.character(chr),] %>% drop_na(start.LD.pos)
  if (nrow(df.tmp>0)){
  dd<- c()
  for (x in 1:nrow(df.tmp)){
    dd<- c(dd,seq(df.tmp$start.LD.pos[x],df.tmp$stop.LD.pos[x], 1e6))
  }
  res<-as.data.frame(table(cut_interval(dd,length=1e6)))
  
  #res <- as.data.frame(table(cut_width(df.tmp$start/1e6,1, boundary = 0)))
  
  res$start <- as.numeric(as.character(substr(res$Var1,2,regexpr(",",res$Var1)-1)))
  res$end <- res$start +1e6
  res$chr <- chr
  output <- rbind(output,res)}
}



mean_freq <- mean(output$Freq, na.rm=T)

output <- output %>% filter(Freq != 0)
write_delim(output, "frequencies_intervals_SW.csv", delim=";")



chrSizes <- mm10[1:20,]
chrSizes$chr<- 1:20

output$start<- as.numeric(as.character(output$start))/1e6
output$end<- as.numeric(as.character(output$end))/1e6
library(viridis)
snp_plot <- ggplot()+
  geom_segment(data = chrSizes,
               aes(x = chr, xend = chr, y = 0, yend = chr_len/1e6),
               lineend = "round", color = "#DEE1E3", size = 5) + 
  
  geom_segment(data=output, 
               aes(x=chr, xend=chr, y=start, yend=end, color=Freq), size=5)+
  coord_flip() + scale_color_viridis_c(option="magma", direction = -1, limits=c(0,30)) + 
  labs(y="Position (Mb)", x= "Chromosome", color="Overlap count") +theme_test()
interval_plot<-snp_plot + scale_x_discrete(limits=factor(goodChrOrder))
interval_plot
ggsave("intervals_DNA_RNA_plot_SW.pdf")
# study-wide
SW_threshold =0.05/32625/118.2177

median_interval_size <- median(input_file$length, na.rm=T)

unique_int <- all_dna_rna %>% 
  filter(peak.P.type<SW_threshold) %>% 
  group_by(chr) %>% 
  distinct(start.LD.pos, stop.LD.pos)
median_uniq <- median(unique_int$stop.LD.pos-unique_int$start.LD.pos, nr.rm=T)/1e6
un_marker<- unique(all_dna_rna$peak.snp)

### combine with plot from density_results_snps.R
load("../../all_genes_snps/snp_plot.Rdata")
library(patchwork)
interval_plot+snp_plot + plot_annotation(tag_levels = "A")
ggsave("../density_results_intervalsANDsnps.pdf")



# get genes in intervals 

gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www") # can take long
for (i in 1:nrow(output)){
  chr <- output$chr[i]
  start <- output$start[i]
  stop <- output$end[i] 
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
top50<- output %>% slice_max(Freq,n=50)

for (i in 1:nrow(top50)){
  chr <- top50$chr[i]
  start <- top50$start[i]
  stop <- top50$end[i] 
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
    top50$protein.coding.genes[i]<- paste(out.bm.genes.region$external_gene_name, collapse=',')
  }
  if (nrow(out.bm.genes.region.all)>0){
    top50$all.genes[i] <- paste(out.bm.genes.region.all$external_gene_name, collapse = ',')
    top50$all.genes.biotype[i]<- paste(out.bm.genes.region.all$gene_biotype, collapse = ',')
  }
  
}
write_delim(top50, file="top50_overlaps.csv", delim=";")
