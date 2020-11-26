
#################################################################################
##                  Plot region script and write results table                 ##
##   This script will make a locuszoom-like plot for the significant regions   ##
##            using Bonferroni as the default significance threshold           ##
##      It plots the significant SNPs within their confidance interval,        ##
##  which is calculted through finding snps in high LD with peak SNP (LD>0.9)  ##
#################################################################################


setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Sterility/")

run_pretty_region <- function(directory,trait, Xtrue=T, method="bon"){
  source("~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Plot_scripts/New_pretty_plot_script.r")
  setwd(directory)
  
  out_P <- plot.pretty.region(trait, Xtrue, title="Total SNP effect", method=method)
  out_add.P <- plot.pretty.region(trait, F,P="add.P", title="Additive effect", method=method)
  out_dom.P <- plot.pretty.region(trait, F, P="dom.P", title="Dominant effect", method=method)
  total_out1 <- rbind(out_P, out_add.P)
  total_out2 <- rbind(total_out1, out_dom.P)
  write.table(total_out2, file=paste0("../results/",trait,"_results.txt"))
  setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Sterility/")
  #setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacula/")
}

run_pretty_region("./apo_per_tubule/", "apo_per_tubule")
run_pretty_region("./arb_per_tubule/", "arb_per_tubule")
run_pretty_region("./body_weight/", "body_weight_with_covar_body_length") # to do, but too much
run_pretty_region("./combined_testis_weight/", "combined_testis_weight_with_covar_body_weight")
run_pretty_region("./deg_per_tubule/", "deg_per_tubule")
run_pretty_region("./exf_per_tubule/", "exf_per_tubule")
run_pretty_region("./mgc_per_tubule/", "mgc_per_tubule")
run_pretty_region("./relative_testis_weight/", "relative_testis_weight")
run_pretty_region("./sperm_count/", "sperm_count")
run_pretty_region("./sperm_count/", "sperm_count_with_covar_combined_testis_weight")
run_pretty_region("./spermatids_per_perimeter/", "spermatids_per_perimeter")
run_pretty_region("./spermatids_per_perimeter/", "spermatids_per_perimeter_with_covar_spermatocytes_per_perimeter")
run_pretty_region("./spermatocytes_per_perimeter/", "spermatocytes_per_perimeter")
run_pretty_region("./tubule_area/", "tubule_area")
run_pretty_region("./vac_per_tubule/", "vac_per_tubule")
run_pretty_region("./body_length/", "body_length") 
run_pretty_region("./tail_length_no_outliers/", "tail_length_no_outliers")
#fdr

run_pretty_region("./apo_per_tubule/", "apo_per_tubule", method="fdr")
run_pretty_region("./arb_per_tubule/", "arb_per_tubule", method="fdr")
run_pretty_region("./body_weight/", "body_weight_with_covar_body_length", method="fdr") # to do, but too much
run_pretty_region("./combined_testis_weight/", "combined_testis_weight_with_covar_body_weight", method="fdr")
run_pretty_region("./deg_per_tubule/", "deg_per_tubule", method="fdr")
run_pretty_region("./exf_per_tubule/", "exf_per_tubule", method="fdr")
run_pretty_region("./mgc_per_tubule/", "mgc_per_tubule", method="fdr")
run_pretty_region("./relative_testis_weight/", "relative_testis_weight", method="fdr")
run_pretty_region("./sperm_count/", "sperm_count", method="fdr")
run_pretty_region("./sperm_count/", "sperm_count_with_covar_combined_testis_weight", method="fdr")
# run_pretty_region("./spermatids_per_perimeter/", "spermatids_per_perimeter", method="fdr") # too many
run_pretty_region("./spermatids_per_perimeter/", "spermatids_per_perimeter_with_covar_spermatocytes_per_perimeter", method="fdr")
# run_pretty_region("./spermatocytes_per_perimeter/", "spermatocytes_per_perimeter", method="fdr") # too many
run_pretty_region("./tubule_area/", "tubule_area", method="fdr")
run_pretty_region("./vac_per_tubule/", "vac_per_tubule", method="fdr")
run_pretty_region("./bmi/", "bmi", method="fdr")
# run_pretty_region("./body_length/", "body_length", method="fdr") # too many
# run_pretty_region("./tail_length_no_outliers/", "tail_length_no_outliers", method="fdr") # too many

run_pretty_region("./sterPCA1_PC1/","sterPCA1_PC1", method="fdr")
run_pretty_region("sterPCA1_PC2","sterPCA1_PC2", method="fdr")
run_pretty_region("sterPCA1_PC3","sterPCA1_PC3", method="fdr")
run_pretty_region("sterPCA1_PC4","sterPCA1_PC4", method="fdr")
run_pretty_region("sterPCA1_PC5","sterPCA1_PC5", method="fdr")
run_pretty_region("sterPCA1_PC6","sterPCA1_PC6", method="fdr")
run_pretty_region("sterPCA2_PC1", "sterPCA2_PC1", method="fdr")
run_pretty_region("sterPCA2_PC2","sterPCA2_PC2", method="fdr")
run_pretty_region("sterPCA2_PC3","sterPCA2_PC3", method="fdr")
run_pretty_region("sterPCA2_PC4","sterPCA2_PC4", method="fdr")
run_pretty_regiono("sterPCA2_PC5","sterPCA2_PC5", method="fdr")
run_pretty_region("sterPCA2_PC6","sterPCA2_PC6", method="fdr")
run_pretty_region("sterPCA2_PC7","sterPCA2_PC7", method="fdr")
run_pretty_region("sterPCA2_PC8","sterPCA2_PC8", method="fdr")

run_pretty_region("sterPCA3_PC1","sterPCA3_PC1", method="fdr")
run_pretty_region("sterPCA3_PC2", "sterPCA3_PC2", method="fdr")
run_pretty_region("sterPCA3_PC3","sterPCA3_PC3", method="fdr")
run_pretty_region("sterPCA3_PC4", "sterPCA3_PC4", method="fdr")
run_pretty_region("sterPCA3_PC5","sterPCA3_PC5", method="fdr")
run_pretty_region("sterPCA3_PC6","sterPCA3_PC6", method="fdr")
run_pretty_region("sterPCA3_PC7", "sterPCA3_PC7", method="fdr")
run_pretty_region("sterPCA3_PC8","sterPCA3_PC8", method="fdr")
run_pretty_region("sterPCA3_PC9","sterPCA3_PC9", method="fdr")
run_pretty_region("sterPCA3_PC10","sterPCA3_PC10", method="fdr")

run_pretty_region("sterPCA4_PC1","sterPCA4_PC1", method="fdr")
run_pretty_region("sterPCA4_PC2", "sterPCA4_PC2", method="fdr")
run_pretty_region("sterPCA4_PC3", "sterPCA4_PC3", method="fdr")
run_pretty_region("sterPCA4_PC4","sterPCA4_PC4", method="fdr")
run_pretty_region("sterPCA4_PC5", "sterPCA4_PC5", method="fdr")
run_pretty_region("sterPCA4_PC6","sterPCA4_PC6", method="fdr")

run_pretty_region("sterPCA5_PC1","sterPCA5_PC1", method="fdr")
run_pretty_region("sterPCA5_PC2", "sterPCA5_PC2", method="fdr")
run_pretty_region("sterPCA5_PC3", "sterPCA5_PC3", method="fdr")
run_pretty_region("sterPCA5_PC4", "sterPCA5_PC4", method="fdr")
run_pretty_region("sterPCA5_PC5", "sterPCA5_PC5", method="fdr")
run_pretty_region("sterPCA5_PC6", "sterPCA5_PC6", method="fdr")
run_pretty_region("sterPCA5_PC7","sterPCA5_PC7", method="fdr")
run_pretty_region("sterPCA5_PC8","sterPCA5_PC8", method="fdr")

run_pretty_region("fecal_carb_perc","fecal_carb_perc", method="fdr")
run_pretty_region("fecal_protein_perc", "fecal_protein_perc", method="fdr")


setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacula/")
run_pretty_region("bac_PC01/", "bac_PC01")
run_pretty_region("bac_PC02/", "bac_PC02")
run_pretty_region("bac_PC03/", "bac_PC03")
run_pretty_region("bac_PC04/", "bac_PC04") # too many
run_pretty_region("bac_PC05/", "bac_PC05") # too many
run_pretty_region("bac_PC06/", "bac_PC06")
run_pretty_region("bac_PC07/", "bac_PC07")
run_pretty_region("bac_PC08/", "bac_PC08")
run_pretty_region("bac_PC09/", "bac_PC09")
run_pretty_region("bac_PC10/", "bac_PC10")
run_pretty_region("bac_PC11/", "bac_PC11")
run_pretty_region("bac_PC12/", "bac_PC12")

run_pretty_region("bac2d_area/", "bac2d_area_with_covar_body_length")
run_pretty_region("bac2d_base_width/", "bac2d_base_width_with_covar_body_length")
run_pretty_region("bac2d_length/", "bac2d_length_with_covar_body_length")
run_pretty_region("bac2d_shaft_width/", "bac2d_shaft_width_with_covar_body_length")
run_pretty_region("bac_surf_area/", "bac_surf_area_with_covar_body_length")
run_pretty_region("bac_vol/", "bac_vol_with_covar_body_length")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/")
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/")
run_pretty_region("alpha/", "chao1", method="fdr")
run_pretty_region("alpha/", "shannon", method="fdr")
