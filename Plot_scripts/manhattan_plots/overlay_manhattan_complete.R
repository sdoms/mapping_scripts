setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/")

run_manhattan <- function(trait, Xtrue=T){
  require(readr)
  require(ggrepel)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(RColorBrewer)
  # source('~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Plot_scripts/pretty_manhattan_overlay.R')
  source('~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Plot_scripts/annotated_manhattan.R')
  if (trait =="RNA"){
    tax_table <- read_delim("../../Scripts/Summary_tables/trait.list.rna.csv", delim=";")
    
  } else if (trait =="DNA"){
    tax_table <- read_delim("../../Scripts/Summary_tables/trait.list.dna.csv", delim=";")
    
  } else if (trait == "all"){
    tax_table_rna <- read_delim("../../Scripts/Summary_tables/trait.list.rna.csv", delim=";")
    tax_table_dna <- read_delim("../../Scripts/Summary_tables/trait.list.dna.csv", delim=";")
    tax_table <- rbind(tax_table_dna, tax_table_rna)
  } else {
    cat("incorrect trait given")
  }

  
  #gws <- "C_Bacilli"
  #Xtrue <- TRUE
  chrom_range <- c(1:19)
  colnames(tax_table)<- tolower(colnames(tax_table))
  
  # Initialize 
  tax_levels <- c("otu", "genus", "family", "order", "class", "phylum")
  
    gws <-as.character(tax_table[1,tax_levels[1]])
    gwscan <- data.frame()
    for (chr in 1:19){
      tx <- readRDS(paste0(trait,"/", tax_levels[1], "/",gws, "_chr_", chr, "with_add_dom.rds"))
      gwscan <- rbind(gwscan, tx)
    }
    if (file.exists(paste0(trait,"/", tax_levels[1],"/",gws, "_chrX.rds"))){
      X_chrom <- readRDS(paste0(trait,"/", tax_levels[1],"/",gws, "_chrX.rds"))
      #X_chrom$chr <- "20"
      gwscan <- dplyr::bind_rows(gwscan, X_chrom)
      chrom_range <- c(1:19, "X")
      
    }
    dat_gws <- gwscan %>% 
      dplyr::select(marker, chr, pos, P,add.P, dom.P, index )
   
    for (tax_level in tax_levels){ 
      tax_table2 <- tax_table %>% drop_na(tax_level)
    for (gws in pull(tax_table2,tax_level)){
      taxa <- data.frame()
      for (chr in 1:19){
        tx <- readRDS(paste0(trait,"/", tax_level, "/",gws, "_chr_", chr, "with_add_dom.rds"))
        
        tx <- tx %>% 
          dplyr::select(marker, P,add.P, dom.P)
        taxa <- rbind(taxa, tx)
        #dat_gws <- merge(dat_gws, tx, all.x=T, by="marker")
      }
      if (file.exists(paste0(trait,"/", tax_level,"/",gws, "_chrX.rds"))){
        X_chrom <- readRDS(paste0(trait,"/", tax_level,"/",gws, "_chrX.rds"))
        #X_chrom$chr <- "20"
        X_chrom <- X_chrom %>% dplyr::select(marker, P)
        taxa_X <- dplyr::bind_rows(taxa, X_chrom)
        
        chrom_range <- c(1:19, "X")
        
      } else{
        X_chrom<- dat_gws %>% 
          filter(chr=="X") %>% 
          dplyr::select(marker...1) %>% 
          mutate(P=NA)
        taxa_X <- bind_rows(taxa, X_chrom)
      }
      taxa_X$P<- as.numeric(as.character(taxa_X$P))
      taxa_X$add.P<- as.numeric(as.character(taxa_X$add.P))
      taxa_X$dom.P<- as.numeric(as.character(taxa_X$dom.P))
      dat_gws <- dplyr::bind_cols(dat_gws, taxa_X )
    }
    }
    P_cols <- grep("^P", colnames(dat_gws))
    addP_cols <- grep("^add.P", colnames(dat_gws))
    domP_cols <- grep("^dom.P", colnames(dat_gws))
    
    min_tax <- dat_gws %>% 
      rowwise() %>% 
      mutate(minP = min(c_across(all_of(P_cols)),na.rm = T)) %>% 
      mutate(minaddP = min(c_across(all_of(addP_cols)), na.rm = T)) %>% 
      mutate(mindomP = min(c_across(all_of(domP_cols)), na.rm = T)) %>%  
      dplyr::select(marker=marker...1, chr, pos,minP, minaddP, mindomP )
    
 
  
  
  dir.create(path="./overlay_manhattan_plots", showWarnings = F)
  outputdir <- "./overlay_manhattan_plots/"
  
  
  mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
  require(patchwork)
  
  p1 <-gg.manhattan(gwscan = min_tax, SNP="marker", CHR = "chr", BP="pos", 
                    P="minP",col = mypalette, title = paste0(trait, " total SNP effect"))
  
  p2<- gg.manhattan(min_tax,SNP="marker", CHR = "chr", BP="pos", P="minaddP", col = mypalette, title = paste0(trait, " additive SNPs"))
  
  
  p3<-gg.manhattan(min_tax,SNP="marker", CHR = "chr", BP="pos", P="mindomP", col = mypalette, title =paste0(trait, " dominant SNPs"))
  p <- p1 /(p2+p3)
  p
  ggsave(paste0(outputdir, trait, ".pdf"), p, width = 12, height = 8)
  ggsave(paste0(outputdir, trait, ".svg"), p, width = 12, height = 8)
  
  ggsave(paste0(outputdir, trait, "_total.pdf"), p1,  width = 7, height = 8)
  ggsave(paste0(outputdir, trait, "_add.pdf"), p2,  width = 7, height = 8)
  ggsave(paste0(outputdir, trait, "_dom.pdf"), p3,  width = 7, height = 8)
  ggsave(paste0(outputdir, trait, ".png"), p, width = 12, height = 8, dpi=300)
  ggsave(paste0(outputdir, trait, "_total.png"), p1,  width = 7, height = 8, dpi=300)
  ggsave(paste0(outputdir, trait, "_add.png"), p2,  width = 7, height = 8, dpi=300)
  ggsave(paste0(outputdir, trait, "_dom.png"), p3,  width = 7, height = 8, dpi=300)
}
run_manhattan("DNA")
run_manhattan("RNA")

