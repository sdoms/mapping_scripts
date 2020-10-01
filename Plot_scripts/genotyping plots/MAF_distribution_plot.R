setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/")

maf_freq_unfil <- read.table("unfiltered_maf.out.frq", header =TRUE, as.is=T)
maf_freq_fil <- read.table("filtered_maf.out.frq", header =TRUE, as.is=T) 

maf_freq_fil$filtered <- "filtered"
maf_freq_unfil$filtered <- "unfiltered"
maf_freq <- rbind(maf_freq_unfil, maf_freq_fil)

library(ggplot2)
library(ggsci)
ggplot(maf_freq, aes(MAF, color=filtered)) + geom_freqpoly() + theme_bw() + labs(y="Number of SNPs", color="") +scale_color_tron()
ggsave("MAF_distribution_unfiltered_filtered.pdf")
