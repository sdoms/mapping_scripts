setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Protein_Carb_content/")

trait <- "genus"
DNAorRNA<-"DNA"
outputdir<- "./output"

library(ggplot2)
library(ggpubr)
corTestAndPlot <- function(trait, DNAorRNA, outputdir="./output"){
  dir.create(outputdir, showWarnings = F)
pheno <- read.csv("NewPhenotypesCombined20201120.csv", sep=";")  
pheno$SampleID<- gsub(x = pheno$SampleID, pattern = "\\/", replacement = ".")  
otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/",trait, "_table_f2_core.csv"))
rownames(otu_table)<- otu_table$X
otu_table <- otu_table[which(grepl(paste0("_",DNAorRNA),rownames(otu_table))),] # only do the DNA
rownames(otu_table) <- substr(rownames(otu_table), 1,15)
otu_table$X<- NULL
otu_table$SampleID<- rownames(otu_table)
intable <- merge(pheno[,2:4], otu_table, by="SampleID")

otu_table$SampleID <- NULL
#intable <- intable[-113,] # remove the outlier with 2 and 7 (was the only one)

output_table <- data.frame()
for (tax in colnames(otu_table)){
  cp<-cor.test(intable$fecal_protein_perc,intable[,tax], use = "complete.obs", method="spearman", exact = F)
  output_table[tax,"fecal_protein_perc_p-val"]<-cp$p.value
  output_table[tax,"fecal_protein_perc_corr"]<-cp$estimate
  if (cp$p.value<0.05){
    ggscatter(intable, x = "fecal_protein_perc", y = tax, 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman", xlab="Fecal protein percentage", ylab="Relative abundance", title = tax)
    ggsave(paste0(outputdir,"/", DNAorRNA, "_",trait,"_", tax, "_protein_perc.pdf"))
  }
  cc<-cor.test(intable$fecal_carb_perc,intable[,tax], use = "complete.obs", method = "spearman", exact = F)
  if (cc$p.value<0.05){
    ggscatter(intable, x = "fecal_carb_perc", y = tax, 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "spearman", xlab="Fecal carbohydrate percentage", ylab="Relative abundance", title = tax)
    ggsave(paste0(outputdir,"/", DNAorRNA, "_",trait,"_", tax, "_carb_perc.pdf"))
  }
  output_table[tax,"fecal_carb_perc_p-val"]<-cc$p.value
  output_table[tax,"fecal_carb_perc_corr"]<-cc$estimate

}
output_table<- cbind(taxa=rownames(output_table),output_table)
write_excel_csv2(output_table, file = paste0(DNAorRNA,"_", trait, ".csv"))
}


corTestAndPlot("phylum", "DNA")
corTestAndPlot("phylum", "RNA")
corTestAndPlot("class", "DNA")
corTestAndPlot("class", "RNA")
corTestAndPlot("order", "DNA")
corTestAndPlot("order", "RNA")
corTestAndPlot("family", "DNA")
corTestAndPlot("family", "RNA")
corTestAndPlot("genus", "DNA")
corTestAndPlot("genus", "RNA")
corTestAndPlot("otu", "DNA")
corTestAndPlot("otu", "RNA")
# 

# ggscatter(intable, x = "fecal_protein_perc", y = "P_Bacteroidetes", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "spearman")
# 
# ggscatter(intable, x = "fecal_carb_perc", y = "P_Proteobacteria", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "spearman")

lmod <- lm(G_Butyricicoccus~fecal_protein_perc+fecal_carb_perc, data=intable)
summary(lmod)





