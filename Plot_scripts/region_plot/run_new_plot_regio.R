source('~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Plot_scripts/plot_region.R')
require(biomaRt)
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version=102) # can take long
# load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_intervals/SW/genes_DNA_RNA_intervals_SW.Rdata")
load(file="/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_18032021_SW.Rdata")


for (i in 1:nrow(all_dna_rna_SW)){
  plot.region(trait = as.character(all_dna_rna_SW[i,"trait"]), tax_level = as.character(all_dna_rna_SW[i,"tax_level"]), 
              chr = as.character(all_dna_rna_SW[i,"chr.num"]), start = as.numeric(all_dna_rna_SW[i,"start.LD.pos"]), 
              stop = as.numeric(all_dna_rna_SW[i,"stop.LD.pos"]), 
              peak_snp = as.character(all_dna_rna_SW[i,"peak.snp"]), P=as.character(all_dna_rna_SW[i,"P.type"]), 
              dnaorrna = as.character(all_dna_rna_SW[i,"dna.rna"]), 
              outputdir = "~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_new_fig_DNA_RNA/", gene.ensembl)
}

# plot.region(trait = "SV184", tax_level = "otu",chr = 4, start = 67067748,stop = 72515350,peak_snp = "UNC7414459", P="add.P", dnaorrna = "RNA",outputdir = "~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_new_fig_DNA_RNA/", gene.ensembl)
# plot.region(trait = "SV184", tax_level = "otu",chr = 15, start = 94362391,stop = 94362391,peak_snp = "UNC26145702", P="add.P", dnaorrna = "RNA",outputdir = "~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_new_fig_DNA_RNA/", gene.ensembl)
