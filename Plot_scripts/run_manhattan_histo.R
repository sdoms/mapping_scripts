run_manhattan <- function(trait, Xtrue=TRUE, covar ="none", method="bon"){
  require(tidyverse)
  require(dplyr)
  source('~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Plot_scripts/pretty_manhattan.r')
  dir.create(path="./manhattan_plots", showWarnings = F)
  outputdir <- "./manhattan_plots/"
  taxa <- data.frame()
  if (covar!="none"){
    path <- paste0(trait , "/",trait,"_with_covar_", covar, "_chr_")
    pathX <- paste0(trait , "/",trait,"_with_covar_", covar, "_chrX.rds")
    outPath <- paste0(outputdir, trait,"_with_covar_", covar)
  } else {
    path <- paste0(trait , "/",trait, "_chr_")
    pathX <- paste0(trait ,"/",trait, "_chrX.rds")
    outPath <- paste0(outputdir, trait)
  }
  for (chr in 1:19){
    tx <- readRDS(paste0(path, chr, "with_add_dom.rds"))
    taxa <- rbind(taxa, tx)
  }
  if (Xtrue){
    X_chrom <- readRDS(pathX)
    X_chrom$chr <- "20"
    taxa_X <- dplyr::bind_rows(taxa, X_chrom)
  } else{
    taxa_X <- taxa
  }
  

  mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
  require(patchwork)
  
  p1 <-gg.manhattan(gwscan = taxa_X, SNP="marker", CHR = "chr", BP="pos", P="P",col = mypalette, title = paste0(trait, " total SNP effect"), method=method)
  
  p2<- gg.manhattan(taxa,SNP="marker", CHR = "chr", BP="pos", P="add.P", col = mypalette, title = paste0(trait, " additive SNPs"), method=method)
  
  
  p3<-gg.manhattan(taxa,SNP="marker", CHR = "chr", BP="pos", P="dom.P", col = mypalette, title =paste0(trait, " dominant SNPs"), method=method)
  p <- p1 /(p2+p3)
  ggsave(paste0(outPath,"_",method, ".pdf"), p, width = 12, height = 8)
  ggsave(paste0(outPath,"_",method, "_total.pdf"), p1,  width = 7, height = 8)
  ggsave(paste0(outPath,"_",method, "_add.pdf"), p2,  width = 7, height = 8)
  ggsave(paste0(outPath, "_", method,"_dom.pdf"), p3,  width = 7, height = 8)
}

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Sterility/")

run_manhattan("apo_per_tubule")
run_manhattan("arb_per_tubule")
run_manhattan("body_weight", covar = "body_length")
run_manhattan("relative_testis_weight")
run_manhattan("combined_testis_weight", covar="body_weight")
run_manhattan("deg_per_tubule")
run_manhattan("exf_per_tubule")
run_manhattan("mgc_per_tubule")
run_manhattan("sperm_count")
run_manhattan("sperm_count", covar="combined_testis_weight")
run_manhattan("spermatids_per_perimeter")
run_manhattan("spermatocytes_per_perimeter")
run_manhattan("tubule_area")
run_manhattan("vac_per_tubule")

#fdr
run_manhattan("apo_per_tubule", method="fdr")
run_manhattan("arb_per_tubule", method="fdr")
run_manhattan("body_weight", covar = "body_length", method="fdr")
run_manhattan("relative_testis_weight", method="fdr")
run_manhattan("combined_testis_weight", covar="body_weight", method="fdr")
run_manhattan("deg_per_tubule", method="fdr")
run_manhattan("exf_per_tubule", method="fdr")
run_manhattan("mgc_per_tubule", method="fdr")
run_manhattan("sperm_count", method="fdr")
run_manhattan("sperm_count", covar="combined_testis_weight", method="fdr")
run_manhattan("spermatids_per_perimeter", method="fdr")
run_manhattan("spermatocytes_per_perimeter", method="fdr")
run_manhattan("tubule_area", method="fdr")
run_manhattan("vac_per_tubule", method="fdr")


setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacula/")
run_manhattan("bac_PC01")
run_manhattan("bac_PC02")
run_manhattan("bac_PC03")
run_manhattan("bac_PC04")
run_manhattan("bac_PC05")
run_manhattan("bac_PC06")
run_manhattan("bac_PC07")
run_manhattan("bac_PC08")
run_manhattan("bac_PC09")
run_manhattan("bac_PC10")
run_manhattan("bac_PC11")
run_manhattan("bac_PC12")

run_manhattan("bac2d_area", covar="body_length")
run_manhattan("bac_surf_area", covar="body_length")
run_manhattan("bac2d_area", covar="body_length")
run_manhattan("bac_vol", covar="body_length")
run_manhattan("bac2d_base_width", covar="body_length")
run_manhattan("bac2d_length", covar="body_length")
run_manhattan("bac2d_shaft_width", covar="body_length")




#fdr
run_manhattan("bac_PC01", method="fdr")
run_manhattan("bac_PC02", method="fdr")
run_manhattan("bac_PC03", method="fdr")
run_manhattan("bac_PC04", method="fdr")
run_manhattan("bac_PC05", method="fdr")
run_manhattan("bac_PC06", method="fdr")
run_manhattan("bac_PC07", method="fdr")
run_manhattan("bac_PC08", method="fdr")
run_manhattan("bac_PC09", method="fdr")
run_manhattan("bac_PC10", method="fdr")
run_manhattan("bac_PC11", method="fdr")
run_manhattan("bac_PC12", method="fdr")

run_manhattan("bac2d_area", covar="body_length", method="fdr")
run_manhattan("bac_surf_area", covar="body_length", method="fdr")
run_manhattan("bac2d_area", covar="body_length", method="fdr")
run_manhattan("bac_vol", covar="body_length", method="fdr")
run_manhattan("bac2d_base_width", covar="body_length", method="fdr")
run_manhattan("bac2d_length", covar="body_length", method="fdr")
run_manhattan("bac2d_shaft_width", covar="body_length", method="fdr")

