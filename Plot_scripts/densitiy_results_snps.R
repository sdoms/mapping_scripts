setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")

load("../../../../Cleaning_snps/clean_snps.Rdata")

library(tidyverse)
# format df
chrSizes <- snps %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos))

# put the chromosomes in the good order: chr1, chr2, chr22, chrX
goodChrOrder <- c(1:19,"X")


# results file 


inputfile <- read.csv("all_markers_RNA_DNA-with-genes_SW.csv", sep=";")
output <- data.frame()
for (chr in goodChrOrder){
  df.tmp <- inputfile[inputfile$chr==chr,]
  if(nrow(df.tmp)>0){
    res <- as.data.frame(table(cut_width(df.tmp$pos,1, boundary = 0)))
    res$start <- grep("^[/(/[][[:digit:]]", res$Var1)
    res$end <- res$start +1
    res$chr <- chr
    output <- rbind(output,res)
  }

}
mean_freq <- mean(output$Freq, na.rm=T)

output <- output %>% filter(Freq != 0)

chrSizes <- na.omit(chrSizes)
library(viridis)
snp_plot <- ggplot()+
  geom_segment(data = chrSizes,
               aes(x = chr, xend = chr, y = 0, yend = chr_len),
               lineend = "round", color = "#DEE1E3", size = 5) + 
  
  geom_segment(data=output, 
               aes(x=chr, xend=chr, y=start, yend=end, color=Freq), size=4)+
  coord_flip() + scale_color_viridis(option="magma", direction = -1) + 
  labs(y="Position (Mb)", x= "Chromosome", color="Times SNP significant") +theme_test()
snp_plot<-snp_plot + scale_x_discrete(limits=goodChrOrder)
snp_plot
ggsave("SNP_density_results_SW.pdf")
save(snp_plot, file="snp_plot.Rdata")
library(patchwork)


#### genome-wide

inputfile <- read.csv("all_markers_RNA_DNA-with-genes_GW.csv", sep = ";")
output <- data.frame()
for (chr in goodChrOrder){
  df.tmp <- inputfile[inputfile$chr==chr,]
  res <- as.data.frame(table(cut_width(df.tmp$pos,1, boundary = 0)))
  res$start <- grep("^[/(/[][[:digit:]]", res$Var1)
  res$end <- res$start +1
  res$chr <- chr
  output <- rbind(output,res)
}
mean_freq <- mean(output$Freq, na.rm=T)

output <- output %>% filter(Freq != 0)

chrSizes <- na.omit(chrSizes)
library(viridis)
snp_plot <- ggplot()+
  geom_segment(data = chrSizes,
               aes(x = chr, xend = chr, y = 0, yend = chr_len),
               lineend = "round", color = "#DEE1E3", size = 5) + 
  
  geom_segment(data=output, 
               aes(x=chr, xend=chr, y=start, yend=end, color=Freq), size=4)+
  coord_flip() + scale_color_viridis(option="magma", direction = -1) + 
  labs(y="Position (Mb)", x= "Chromosome", color="QTL Count") +theme_test()
snp_plot + scale_x_discrete(limits=goodChrOrder)
ggsave("Results/Plots/SNP_density_results_GW.pdf")
