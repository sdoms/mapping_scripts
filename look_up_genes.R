setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/")
library(tidyverse)
library(stringr)
library(readxl)
all_genes <- read.delim("all_genes_DNA_RNA.txt", header=F)
all_genes$x <- splus2R::upperCase(all_genes$V1)

ASV_DNA<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/SV_results_DNA.xlsx", skip = 3)

genus_DNA<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/genus_results_DNA.xlsx", sheet="all")
family_DNA<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/family_results_DNA.xlsx", sheet="all")
order_DNA<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/order_results_DNA.xlsx", sheet="all")
class_DNA<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/class_results_DNA.xlsx", sheet="all")
phylum_DNA<- read_excel("./DNA/summary_tables/genome_wide_significant/gwscan/phylum_results_DNA.xlsx", sheet="all")


genes_DNA<-rbind(genus_DNA,family_DNA)
genes_DNA<- rbind(family_DNA,genes_DNA)
genes_DNA<-rbind(order_DNA,genes_DNA)
genes_DNA<-rbind(class_DNA,genes_DNA)
genes_DNA<-rbind(phylum_DNA,genes_DNA)
genes_DNA$Genus<-"NA"
genes_DNA<- rbind(genes_DNA,ASV_DNA)
genes_DNA$DNAorRNA <- "DNA"

ASV_RNA<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/SV_results_RNA.xlsx",skip=2, sheet="all")
genus_RNA<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/genus_results_RNA.xlsx", sheet="all")
family_RNA<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/family_results_RNA.xlsx", sheet="all")
order_RNA<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/order_results_RNA.xlsx", sheet="all")
class_RNA<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/class_results_RNA.xlsx")
phylum_RNA<- read_excel("./RNA/summary_tables/genome_wide_significant/gwscan/phylum_results_RNA.xlsx")

genes_RNA<-rbind(genus_RNA,family_RNA)
genes_RNA<- rbind(family_RNA,genes_RNA)
genes_RNA<-rbind(order_RNA,genes_RNA)
genes_RNA<-rbind(class_RNA,genes_RNA)
genes_RNA<-rbind(phylum_RNA,genes_RNA)
genes_RNA$Genus<-"NA"
genes_RNA<- rbind(genes_RNA,ASV_RNA)
genes_RNA$DNAorRNA <- "RNA"

genes_DNA_RNA<-rbind(genes_DNA,genes_RNA)
save(genes_DNA_RNA, file="genes_DNA_RNA.Rdata")


oscilli <- read.delim("Arunabh_genes/all_genes_oscillibacter.txt",sep = "")


for (gene in oscilli$x){
  xx<-str_which(all_genes$x, gene)
  if (length(xx)>0){
    result<-genes_DNA_RNA[str_which(genes_DNA_RNA$list_total_genes,regex(gene,ignore_case = T)),]
    write_delim(result,file=paste0("Arunabh_genes/Results/Oscillibacter/Gene_",gene,"_results.txt"))
  }
}

taxa<- c("faecalibacterium","oscillibacter","prevotella","roseburia","shannon")
for (taxon in taxa){

input<- read.delim(paste0("Arunabh_genes/all_genes_",taxon,".txt"), sep=" ")
dir.create(paste0("Arunabh_genes/Results/",taxon),showWarnings = F)

for (gene in input$x){
  xx<-str_which(all_genes$x, gene)
  if (length(xx)>0){
    result<-genes_DNA_RNA[str_which(genes_DNA_RNA$list_total_genes,regex(gene,ignore_case = T)),]
    #write_delim(result,file=paste0("Arunabh_genes/Results/",taxon,"/Gene_",gene,"_results.txt"))
    write_csv2(result,file=paste0("Arunabh_genes/Results/",taxon,"/Gene_",gene,"_results.csv"))
    }
}
}


setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/protein_overexpression/")
# protein overexpression
input<- read_xlsx("overexpressed_protein_GF_con.xlsx",sheet=2, col_names = F)


for (gene in input$...1){
  xx<-str_which(all_genes$x, gene)
  if (length(xx)>0){
    result<-genes_DNA_RNA[str_which(genes_DNA_RNA$list_total_genes,regex(gene,ignore_case = T)),]
    write_csv2(result,file=paste0("overexpressed_GF_",gene,".csv"))
  }
}



input<- read_xlsx("overexpressed_protein_GF_con.xlsx",sheet=3, col_names = F)

found_genes<-character()
for (gene in input$...1){
  xx<-str_which(all_genes$x, gene)
 
  if (length(xx)>0){
    found_genes<- append(found_genes,gene)
    result<-genes_DNA_RNA[str_which(genes_DNA_RNA$list_total_genes,regex(gene,ignore_case = T)),]
    write_csv2(result,file=paste0("overexpressed_Conventional_",gene,".csv"))
  }
}


# only < 5mb intervals 
input<- read_xlsx("overexpressed_protein_GF_con.xlsx",sheet=2, col_names = F)

overexpressed_con<-data.frame()
for (gene in input$...1){
  xx<-str_which(all_genes$x, gene)
  if (length(xx)>0){
    result<-genes_DNA_RNA[str_which(genes_DNA_RNA$list_total_genes,regex(gene,ignore_case = T)),]
    result$gene <- gene
    overexpressed_con<- rbind(overexpressed_con,result)
  }
}

overexpressed_con_filt <- overexpressed_con %>% 
  filter(as.numeric(length)<5)

gggg<- overexpressed_con_filt %>% 
  distinct(gene)

input<- read_xlsx("overexpressed_protein_GF_con.xlsx",sheet=3, col_names = F)

overexpressed_con_o<-data.frame()
for (gene in input$...1){
  xx<-str_which(all_genes$x, gene)
  if (length(xx)>0){
    result<-genes_DNA_RNA[str_which(genes_DNA_RNA$list_total_genes,regex(gene,ignore_case = T)),]
    result$gene <- gene
    overexpressed_con_o<- rbind(overexpressed_con_o,result)
  }
}

overexpressed_con_o_filt <- overexpressed_con_o %>% 
  filter(as.numeric(length)<5)


gggg_o<- overexpressed_con_o_filt %>% 
  distinct(gene)
dd<-overexpressed_con_o %>% 
  distinct(gene)
ee<- overexpressed_con %>% 
  distinct(gene)
