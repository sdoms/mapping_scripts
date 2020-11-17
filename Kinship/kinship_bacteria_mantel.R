gemma <-"/usr/local/bin/gemma"
cat("Reading map file.\n")
load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata")
snps$cM <- NULL
snps$tier <- NULL
snps$pos <- snps$pos * 10e5

analysis <- "otu"
DNAorRNA <- "DNA"

taxa <- read.pheno(paste0("../../Phenotyping_27.02.2020/tables_core/",analysis,"_table_f2_core.csv"))
# 
taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),] # only do the DNA
rownames(taxa) <- substr(rownames(taxa), 1,15) # remove the _DNA

# pheno table 
write.table(data.frame(rownames(taxa),1),"pheno.txt", quote=F, row.names=F, col.names=F)

gts <- as.data.frame(t(taxa), check.names=F)

added_cols <- data.frame("marker"=rownames(gts), "A1"="A", "A2"="C")
gts <- cbind(added_cols,gts)

write.table(gts,"kinship_taxa.txt",sep = " ",quote = FALSE,row.names = FALSE,
            col.names = FALSE)


# Compute the kinship matrices LOCO
system(paste0(gemma, " -g kinship_taxa.txt -gk 1 -p pheno.txt -debug -outdir kinship_microbes -notsnp"))

kinship_DNA_micro <- read.table("kinship_microbes/result.cXX.txt")
rownames(kinship_DNA_micro)<- colnames(gts)[-c(1,2,3)]
colnames(kinship_DNA_micro)<- colnames(gts)[-c(1,2,3)]

# genotype kinship
F2 <- read.plink("../../Cleaning_snps/clean_f2")
F2 <- recode(F2, "relative")

colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")  

#Select the individuals from the genotype file 
F2 <- as.matrix(as.data.frame(F2))
ind <- rownames(taxa)

individuals <- ncol(F2)
F2_part <- apply(F2[,7:individuals], 2, as.numeric)# make sure that the genotype information(0,1,2) is numeric

F2_part <- F2_part[,colnames(F2_part) %in% ind]

F2 <- cbind(map[,c(1,5,6)],F2_part)

write.table(F2,"geno.txt",sep = " ",quote = FALSE,row.names = FALSE,
            col.names = FALSE)

#write.gemma.map("map.txt", map)

# # compute kinship matrix

system(paste(gemma, "-g geno.txt -gk 1 -p pheno.txt -debug -o kinship_DNA" ))

# load kinship matrix geno
kinship_DNA <- read.table("kinship_DNA.cXX.txt")
rownames(kinship_DNA)<- colnames(F2)[-c(1:3)]
colnames(kinship_DNA)<- colnames(F2)[-c(1:3)]

rownames(kinship_DNA)==rownames(kinship_DNA_micro)
# make sure in right order
kinship_DNA2 <- arrange(kinship_DNA,rownames(kinship_DNA))
kinship_DNA2<- kinship_DNA2[,rownames(kinship_DNA2)]

kinship_DNA_micro2 <- arrange(kinship_DNA_micro,rownames(kinship_DNA_micro))
kinship_DNA_micro2<- kinship_DNA_micro2[,rownames(kinship_DNA_micro2)]
rownames(kinship_DNA2)==rownames(kinship_DNA_micro2)
library(vegan)
m_DNA<-mantel(kinship_DNA2,kinship_DNA_micro2,method="spearman", permutations = 9999)

library(ggplot2)
library(ggpubr)
library(patchwork)
meltKK_DNA <- reshape2::melt(cbind(rownames(kinship_DNA2),kinship_DNA2))
meltMM_DNA <- reshape2::melt(cbind(rownames(kinship_DNA_micro2),kinship_DNA_micro2))
meltMM_DNA$id <- paste0(meltMM_DNA$`rownames(kinship_DNA_micro2)`,"_",meltMM_DNA$variable)
meltKK_DNA$id <- paste0(meltKK_DNA$`rownames(kinship_DNA2)`,"_",meltKK_DNA$variable)
meltTT_DNA <- merge(meltKK_DNA, meltMM_DNA, by="id")

colnames(meltTT_DNA)<- c("id", "idgeno1", "idgeno2","geno_kinship", "idmicro1", "idmicro2", "geno_micro")

g_DNA <- ggplot(data=meltTT_DNA, aes(y=geno_micro, x=geno_kinship))+geom_point(shape=1)  +  geom_smooth(method=lm,   # Add linear regression line
                                                                              se=FALSE) +   # Don't add shaded confidence region
  theme_test()+labs(y="Microbiome-based kinship",x= "Host genotype kinship")+annotate("text",label= paste("Mantel test: p=",m_DNA$signif),x=0.2,y=0.005)
g_DNA
ggsave("DNA_micro_geno_kinship_mantel.pdf")
beta <- as.matrix(vegdist(taxa))
meltTT$fid <- substr(meltTT$id1, 1,7)
rownames(beta)==rownames(kinship_DNA2)
mantel(kinship_DNA2,beta,method="spearman", permutations = 9999)

DNAorRNA <- "RNA"

taxa <- read.pheno(paste0("../../Phenotyping_27.02.2020/tables_core/",analysis,"_table_f2_core.csv"))
# 
taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),] # only do the DNA
rownames(taxa) <- substr(rownames(taxa), 1,15) # remove the _DNA

# pheno table 
write.table(data.frame(rownames(taxa),1),"pheno.txt", quote=F, row.names=F, col.names=F)

gts <- as.data.frame(t(taxa), check.names=F)
added_cols <- data.frame("marker"=rownames(gts), "A1"="A", "A2"="C")
gts <- cbind(added_cols,gts)

write.table(gts,"kinship_taxa.txt",sep = " ",quote = FALSE,row.names = FALSE,
            col.names = FALSE)


# Compute the kinship matrices 
system(paste0(gemma, " -g kinship_taxa.txt -gk 1 -p pheno.txt -debug -outdir kinship_microbes -notsnp"))

# genotype kinship
F2 <- read.plink("../../Cleaning_snps/clean_f2")
F2 <- recode(F2, "relative")

colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")  

#Select the individuals from the genotype file 
F2 <- as.matrix(as.data.frame(F2))
individuals <- ncol(F2)
F2_part <- apply(F2[,7:individuals], 2, as.numeric)# make sure that the genotype information(0,1,2) is numeric

F2_part <- F2_part[,colnames(F2_part) %in% ind]

F2 <- cbind(map[,1:6],F2_part)

write.gemma.geno("geno.txt",F2)
#write.gemma.map("map.txt", map)

# # compute kinship matrix

system(paste(gemma, "-g geno.txt -gk 1 -p pheno.txt -outdir . -o kinship_RNA" ))

kinship_RNA_micro <- read.table("kinship_microbes/result.cXX.txt")
rownames(kinship_RNA_micro)<- colnames(gts)[-c(1,2,3)]
colnames(kinship_RNA_micro)<- colnames(gts)[-c(1,2,3)]
# load kinship matrix geno
kinship_RNA <- read.table("kinship_RNA.cXX.txt")
rownames(kinship_RNA)<- colnames(F2)[-c(1:6)]
colnames(kinship_RNA)<- colnames(F2)[-c(1:6)]

rownames(kinship_RNA)==rownames(kinship_RNA_micro)
# make sure in right order
kinship_RNA2 <- arrange(kinship_RNA,rownames(kinship_RNA))
kinship_RNA2<- kinship_RNA2[,rownames(kinship_RNA2)]

kinship_RNA_micro2 <- arrange(kinship_RNA_micro,rownames(kinship_RNA_micro))
kinship_RNA_micro2<- kinship_RNA_micro2[,rownames(kinship_RNA_micro2)]
rownames(kinship_RNA2)==rownames(kinship_RNA_micro2)
library(vegan)
m<-mantel(kinship_RNA2,kinship_RNA_micro2,method="spearman", permutations = 999)

library(ggplot2)
library(ggpubr)

meltKK <- reshape2::melt(cbind(rownames(kinship_RNA2),kinship_RNA2))
meltMM <- reshape2::melt(cbind(rownames(kinship_RNA_micro2),kinship_RNA_micro2))
meltMM$id <- paste0(meltMM$`rownames(kinship_RNA_micro2)`,"_",meltMM$variable)
meltKK$id <- paste0(meltKK$`rownames(kinship_RNA2)`,"_",meltKK$variable)
meltTT <- merge(meltKK, meltMM, by="id")
colnames(meltTT)<- c("id", "idgeno1", "idgeno2","geno_kinship", "idmicro1", "idmicro2", "geno_micro")
g_RNA <- ggplot(data=meltTT, aes(y=geno_micro, x=geno_kinship))+geom_point(shape=1)  +  geom_smooth(method=lm,   # Add linear regression line
                                                                                       se=FALSE) +   # Don't add shaded confidence region
  theme_test()+labs(y="Microbiome-based kinship",x= "Host genotype kinship")+annotate("text",label= paste("Mantel test: p=",m$signif),x=0.2,y=0.0075)

g_RNA

g_DNA +g_RNA +plot_annotation(tag_levels = "A")
ggsave("kinship_geno_micro_mantel.pdf")

melt_tot <- merge(meltTT, meltTT_DNA, by="id", all.x=T)
meltTT2 <- merge(meltKK, meltMM_DNA, by="id")
colnames(meltTT2)<- c("id", "idgeno1", "idgeno2","geno_kinship", "idmicro1", "idmicro2", "geno_micro")
ggplot(data=meltTT2, aes(y=geno_micro, x=geno_kinship))+geom_point(shape=1)  +  geom_smooth(method=lm,   # Add linear regression line
                                                                                           se=FALSE) +   # Don't add shaded confidence region
  theme_test()+labs(y="Microbiome-based kinship",x= "Host genotype kinship")+annotate("text",label= paste("Mantel test: p=",m$signif),x=0.2,y=0.0075)

