setwd("/Users/turner/Dropbox/BainesTurnerShared2.0/Results/Bacula/")
require(biomaRt)
require(stringr)
require(forcats)
require(RColorBrewer)
require(viridis)
require(plyr)
require(patchwork)
require(tidyverse)



trait_list=read.csv("trait.list.csv")

i=1
j=1



results.comb <- data.frame()
#chrom_range <- c(1:19)
for (chr in 1:19){
  current.chr=readRDS(paste0(trait_list$Trait[i],"/",trait_list$Trait[i],"_chr_",chr,"with_add_dom.rds"))
  results.comb <- rbind(results.comb, current.chr)
}
results.comb$Q=p.adjust(results.comb$P,method = "BH")

  X_chrom <- readRDS(paste0(trait_list$Trait[i],"/",trait_list$Trait[i], "_chrX.rds"))
  #X_chrom$chr <- "20"
  results.comb <- dplyr::bind_rows(results.comb, X_chrom)
  chrom_range <- c(1:19, "X")
  
}

results.comb[results.comb$chr=="X" & results.comb$index<6,]