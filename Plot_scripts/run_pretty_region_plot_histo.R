setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Sterility/")

run_pretty_region <- function(directory,trait, Xtrue=T){
  source("~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Plot_scripts/New_pretty_plot_script.r")
  setwd(directory)
  
  out_P <- plot.pretty.region(trait, Xtrue, title="Total SNP effect")
  out_add.P <- plot.pretty.region(trait, F,P="add.P", title="Additive effect")
  out_dom.P <- plot.pretty.region(trait, F, P="dom.P", title="Dominant effect")
  total_out1 <- rbind(out_P, out_add.P)
  total_out2 <- rbind(total_out1, out_dom.P)
  write.table(total_out2, file=paste0("../results/",trait,"_results.txt"))
  setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Sterility/")
  
}

run_pretty_region("./apo_per_tubule/", "apo_per_tubule")
run_pretty_region("./arb_per_tubule/", "arb_per_tubule")
run_pretty_region("./body_weight/", "body_weight_with_covar_body_length")
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


