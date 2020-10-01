run_manhattan <- function(trait, Xtrue=TRUE, covar ="none"){
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
  
  p1 <-gg.manhattan(gwscan = taxa_X, SNP="marker", CHR = "chr", BP="pos", P="P",col = mypalette, title = paste0(trait, " total SNP effect"))
  
  p2<- gg.manhattan(taxa,SNP="marker", CHR = "chr", BP="pos", P="add.P", col = mypalette, title = paste0(trait, " additive SNPs"))
  
  
  p3<-gg.manhattan(taxa,SNP="marker", CHR = "chr", BP="pos", P="dom.P", col = mypalette, title =paste0(trait, " dominant SNPs"))
  p <- p1 /(p2+p3)
  ggsave(paste0(outPath, ".pdf"), p, width = 12, height = 8)
  ggsave(paste0(outPath, "_total.pdf"), p1,  width = 7, height = 8)
  ggsave(paste0(outPath, "_add.pdf"), p2,  width = 7, height = 8)
  ggsave(paste0(outPath, "_dom.pdf"), p3,  width = 7, height = 8)
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

