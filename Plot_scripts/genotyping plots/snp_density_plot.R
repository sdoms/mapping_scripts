setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/")

load("Cleaning_snps/clean_snps.Rdata")

library(tidyverse)
# format df
chrSizes <- snps %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos))

# put the chromosomes in the good order: chr1, chr2, chr22, chrX
goodChrOrder <- c(1:19,"X","Y", "M")

output <- data.frame()
for (chr in goodChrOrder){
  df.tmp <- snps[snps$chr==chr,]
  res <- as.data.frame(table(cut_width(df.tmp$pos,1, boundary = 0)))
  res$start <- grep("^[/(/[][[:digit:]]", res$Var1)
  res$end <- res$start +1
  res$chr <- chr
  output <- rbind(output,res)
}

output <- output %>% filter(Freq != 0)
chrSizes <- na.omit(chrSizes)
library(viridis)
snp_plot <- ggplot()+
geom_segment(data = chrSizes,
             aes(x = chr, xend = chr, y = 0, yend = chr_len),
             lineend = "round", color = "#f9f9f9", size = 5) + 
  geom_segment(data=output, 
             aes(x=chr, xend=chr, y=start, yend=end, color=Freq), size=4)+
  coord_flip() + scale_color_viridis(option="magma", direction = -1) + labs(y="Position (Mb)", x= "Chromosome", color="SNP Count") +theme_classic()
snp_plot
ggsave("Results/Plots/SNP_density_plot_without_grid.pdf")
