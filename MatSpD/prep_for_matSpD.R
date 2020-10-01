# make correlation matrix 
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/")
otu <- read.csv(file = "otu_table_f2_core.csv", row.names = "X")
genus <- read.csv("genus_table_f2_core.csv")
OG <- cbind(otu, genus)
OG$X <- NULL
family <- read.csv("family_table_f2_core.csv")
OGF <- cbind(OG, family)
OGF$X <- NULL
order <- read.csv("order_table_f2_core.csv")
OGFO <- cbind(OGF, order)
OGFO$X <- NULL
class <- read.csv(("class_table_f2_core.csv"))
OGFOC <- cbind(OGFO, class)
OGFOC$X <- NULL
phylum <- read.csv("phylum_table_f2_core.csv")
OGFOCP <- cbind(OGFOC, phylum)
OGFOCP$X <- NULL




dna <- OGFOCP[which(grepl("DNA",rownames(OGFOCP))),]
rna <- OGFOCP[which(grepl("RNA",rownames(OGFOCP))),]
rownames(dna)<- substr(rownames(dna), 1,nchar(rownames(dna))-4)
rownames(rna)<- substr(rownames(rna), 1,nchar(rownames(rna))-4)
colnames(rna)<-paste0(colnames(rna), "_RNA")
#remove SV9 in rna 
rna[,"SV9_RNA"] <- NULL
taxa <- merge(dna, rna, all.y=T, by="row.names")
rownames(taxa)<- taxa$Row.names
taxa$Row.names <- NULL
correlation.matrix <- cor(taxa, method="spearman", use="complete.obs")

write.table(correlation.matrix, file="correlation.matrix.txt", row.names = F, col.names = F)
