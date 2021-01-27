source('~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Plot_scripts/plot_region.R')
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = 'www') # can take long

for (i in 1:nrow(input_SW)){
  plot.region(trait = as.character(input_SW[i,"trait"]), tax_level = as.character(input_SW[i,"tax_level"]), 
              chr = as.character(input_SW[i,"chr.num"]), start = as.numeric(input_SW[i,"start.LD.pos"]), 
              stop = as.numeric(input_SW[i,"stop.LD.pos"]), 
              peak_snp = as.character(input_SW[i,"peak.snp"]), P=as.character(input_SW[i,"P.type"]), dnaorrna = "DNA", 
              outputdir = "../new_figures/", gene.ensembl)
}
