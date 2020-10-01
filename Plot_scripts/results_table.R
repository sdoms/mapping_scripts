
make_total_results_file <- function(directory, taxa){
  require(filesstrings)
  require(dplyr)
  #directory <- "../results"
  filenames <- list.files(directory, pattern="*_results.txt", full.names=T)
  df_list <- lapply(filenames, function(x) {
    if (!file.size(x) == 0) {
      read.table(x,header=T)
    }
  })
  ##
  out <- do.call("rbind", df_list)
  out$SV <- sub("_.*", "", out$name)
  if (taxa=="SV"){
  tax_table<- read_csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/otu_tax_table_f2_core.csv")
  tax_table <- tax_table %>%
    select(Genus, X1)
  colnames(tax_table )<- c("Genus", "SV")
  joined_data <- left_join(out, tax_table, by = "SV")
  out <- joined_data
  }
  
  openxlsx::write.xlsx(out, file=paste0(directory,"/",taxa, "_results.xlsx"))
  
  
}

directory <- "~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/results/Phylum"
make_total_results_file(directory, "SV")
make_total_results_file(directory, "genus")
make_total_results_file(directory, "family")
make_total_results_file(directory, "order")
make_total_results_file(directory, "class")
make_total_results_file(directory, "phylum")
