setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/")

pheno <- read.csv("../Protein_Carb_content/NewPhenotypesCombined20201120.csv", sep=";")  
pheno$SampleID<- gsub(x = pheno$SampleID, pattern = "\\/", replacement = ".")  
pheno2<- read.csv("../Phenotypes_histology_new.csv", sep=";")

tax_table <- read.csv(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/phylum_table_f2_core.csv"))
rownames(tax_table)<- tax_table$X
tax_table$F_B_ratio <- tax_table$P_Firmicutes/tax_table$P_Bacteroidetes
tax_table_DNA <- tax_table[which(grepl("_DNA",rownames(tax_table))),] # only do the DNA
rownames(tax_table_DNA) <- substr(rownames(tax_table_DNA), 1,15)
tax_table$X<- NULL
tax_table_DNA$X<- NULL
tax_table_DNA$SampleID<- rownames(tax_table_DNA)

tax_table_RNA <- tax_table[which(grepl("_RNA",rownames(tax_table))),] # only do the RNA
rownames(tax_table_RNA) <- substr(rownames(tax_table_RNA), 1,15)
tax_table_RNA$X<- NULL
tax_table_RNA$SampleID<- rownames(tax_table_RNA)

intable_DNA <- merge(pheno[,2:4], tax_table_DNA, by="SampleID")
intable_RNA <- merge(pheno[,2:4], tax_table_RNA, by="SampleID")
intable_DNA <- merge(pheno2, intable_DNA,by.x="Mouse_Name", by.y="SampleID")
intable_RNA <- merge(pheno2, intable_RNA,by.x="Mouse_Name", by="SampleID")

rownames(intable_DNA)<- intable_DNA$Mouse_Name
intable_DNA<- intable_DNA[,c(5:9, 44:50)]

rownames(intable_RNA)<- intable_RNA$Mouse_Name
intable_RNA<- intable_RNA[,c(5:9, 44:50)]
intable_DNA<- apply(intable_DNA, 2, function(x) as.numeric(as.character(x)))
cp <- Hmisc::rcorr(intable_DNA,type = "spearman")
p_cor<- round(cp$P, 4)

summary(intable_DNA)
