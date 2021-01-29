setwd("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_intervals/SW/")
load("../../../../../Cleaning_snps/clean_snps.Rdata")

load("genes_DNA_RNA_intervals_SW.Rdata")

# system("~/poverlap/poverlap2.py poverlap --a all_intervals_mouse_studies_10mb.txt --b my_intervals_DNA_RNA_10mb.txt -g mouse.mm10.genome --n 9999")
# system("bedtools intersect -wo -a all_intervals_mouse_studies_10mb.txt -b my_intervals_DNA_RNA_10mb.txt >bed_intersect_out.txt")

library(tidyverse)
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

input_file <- input_SW
input_DNA<- input_SW %>% 
  filter(DNAorRNA=="DNA")
input_RNA<- input_SW %>% 
  filter(DNAorRNA=="RNA")
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

unique_int <- input_SW %>% 
  group_by(chr) %>% 
  distinct(start.LD.pos, stop.LD.pos)
median_uniq <- median(unique_int$stop.LD.pos-unique_int$start.LD.pos, nr.rm=T)/1e6
un_marker<- unique(input_SW$peak.snp)

### combine with plot from density_results_snps.R
load("../../all_genes_snps/snp_plot.Rdata")
interval_plot+snp_plot + plot_annotation(tag_levels = "A")
ggsave("../density_results_intervalsANDsnps.pdf")
