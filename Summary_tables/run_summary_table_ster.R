setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Sterility/")
source("~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Summary_tables/summary_function.R")
summary_tables_histo("apo_per_tubule")
summary_tables_histo("arb_per_tubule")
summary_tables_histo("deg_per_tubule")
summary_tables_histo("exf_per_tubule")
summary_tables_histo("mgc_per_tubule")
summary_tables_histo("relative_testis_weight")
summary_tables_histo("sperm_count")
summary_tables_histo("spermatids_per_perimeter")
summary_tables_histo("spermatocytes_per_perimeter")
summary_tables_histo("spermatids_per_tubule")
summary_tables_histo("spermatocytes_per_tubule")
summary_tables_histo("tubule_area")
summary_tables_histo("vac_per_tubule")
summary_tables_histo("body_length",method= "bon")

summary_tables_histo("tail_length_no_outliers", method="bon")
# bmi problem

summary_tables_histo("combined_testis_weight", covar="body_weight")
summary_tables_histo("sperm_count", covar="combined_testis_weight")
summary_tables_histo("spermatids_per_perimeter", covar="spermatocytes_per_perimeter")
summary_tables_histo("body_weight", covar="body_length", method="bon")



setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacula")
summary_tables_histo("bac_PC01")
summary_tables_histo("bac_PC02")
summary_tables_histo("bac_PC03")
summary_tables_histo("bac_PC04")
summary_tables_histo("bac_PC05")
summary_tables_histo("bac_PC06")
summary_tables_histo("bac_PC07")
summary_tables_histo("bac_PC08")
summary_tables_histo("bac_PC09")
summary_tables_histo("bac_PC10")
summary_tables_histo("bac_PC11")
summary_tables_histo("bac_PC12")
summary_tables_histo("bac2d_area", covar="body_length")
summary_tables_histo("bac2d_base_width", covar="body_length")
summary_tables_histo("bac2d_length", covar="body_length")
summary_tables_histo("bac2d_shaft_width", covar="body_length")
summary_tables_histo("bac_surf_area", covar="body_length")
summary_tables_histo("bac_vol", covar="body_length")



summary_tables_histo("bac_PC01", method="bon")
summary_tables_histo("bac_PC02", method="bon")
summary_tables_histo("bac_PC03", method="bon")
summary_tables_histo("bac_PC04", method="bon")
summary_tables_histo("bac_PC05", method="bon")
summary_tables_histo("bac_PC06", method="bon")
summary_tables_histo("bac_PC07", method="bon")
summary_tables_histo("bac_PC08", method="bon")
summary_tables_histo("bac_PC09", method= "bon")
summary_tables_histo("bac_PC10", method= "bon")
summary_tables_histo("bac_PC11", method="bon")
summary_tables_histo("bac_PC12", method="bon")







