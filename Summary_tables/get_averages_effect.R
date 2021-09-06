#allP_all_minmaj<- read_csv2("allP_all_minmaj.csv")
load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_geno.Rdata")
clean_geno <- as.data.frame(clean_geno)
colnames(clean_geno) <- gsub(x = colnames(clean_geno), pattern = "\\/", replacement = ".")  

for (row in 335:nrow(allP_all_minmaj)){
  marker <- as.character(allP_all_minmaj[row, "marker"])
  taxon <- allP_all_minmaj[row, "taxon"]
  DNAorRNA<- allP_all_minmaj[row, "DNAorRNA"]
  taxa <- as.character(allP_all_minmaj[row, "tax"])
  aa_min <- as.character(allP_all_minmaj[row,"AA_min"])
  bb_maj <- as.character(allP_all_minmaj[row,"BB_maj"])
  

  all_mark <-clean_geno[marker,]
  genotypes <- as.data.frame(t(all_mark))
  genotypes <- genotypes %>% 
    dplyr::filter(.[[1]] !="N")
  
  otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/",taxon, "_table_f2_core.csv"))
  rownames(otu_table)<- otu_table$X
  otu_table <- otu_table[which(grepl(paste0("_",DNAorRNA),rownames(otu_table))),] # only do the DNA
  rownames(otu_table) <- substr(rownames(otu_table), 1,15)
  SV <- as.data.frame(otu_table[,taxa])
  rownames(SV) <- rownames(otu_table)
  tot_with_SV <- merge(genotypes, SV, by="row.names")
  colnames(tot_with_SV)<- c("individuals", "marker", "abundance")
  tot_sum <- tot_with_SV %>% 
    group_by(marker) %>% 
    dplyr::summarise(mean_abund=mean(abundance))
  tot_mean<-data.frame( "tot_mean",mean(tot_with_SV$abundance))
  colnames(tot_mean)=c("marker", "mean_abund")
  tot_sum<- rbind(tot_sum, tot_mean)
  tot_sum <- tot_sum %>% 
    mutate(min_maj=ifelse(marker==aa_min, "minor_allele", ifelse(marker==bb_maj, "major_allele",ifelse(marker=="H", "H", NA))))
  # rownames(tot_sum)<- tot_sum$marker
  # tot_sum[aa_min,"min_maj"]<- "minor_allele"
  # rownames(tot_sum)<- tot_sum$marker
  # tot_sum[bb_maj,"min_maj"]<- "major_allele"
  # rownames(tot_sum)<- tot_sum$marker
  # tot_sum["H","min_maj"]<- "H"
  # rownames(tot_sum)<- tot_sum$marker
  if (aa_min %in% tot_sum$marker){
    allP_all_minmaj[row, "min_mean_abund"]<-tot_sum %>% filter(marker==aa_min) %>% pull(mean_abund)
    
  }
  if (bb_maj %in% tot_sum$marker){
    allP_all_minmaj[row, "maj_mean_abund"]<-tot_sum %>% filter(marker==bb_maj) %>% pull(mean_abund)
    
  }
  if ("H" %in% tot_sum$marker){
    allP_all_minmaj[row, "H_mean_abund"]<-tot_sum %>% filter(marker=="H") %>% pull(mean_abund)
    
  }
  allP_all_minmaj[row, "tot_mean_abund"]<-tot_sum %>% filter(marker=="tot_mean") %>% pull(mean_abund)
  allP_all_minmaj[row, "high"]<-tot_sum[which.max(tot_sum$mean_abund[1:3]),"min_maj"]
  allP_all_minmaj[row, "middle"]<-tot_sum[which(!c(1:3) %in% c(which.max(tot_sum$mean_abund[1:3]), which.min(tot_sum$mean_abund[1:3]))),"min_maj"]
  
  allP_all_minmaj[row, "low"]<-tot_sum[which.min(tot_sum$mean_abund[1:3]),"min_maj"]
  
}
write.csv(allP_all_minmaj, file="allP_all_minmaj.csv")
