setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/")

library(lme4qtl)
library(BEDMatrix)
library(MASS)
library(parallel)
library(gtools)
library(lme4)
library(plyr)
library(lmerTest)
library(car)
library(MuMIn)
#library(qqman)

# Parameters ####
DNAorRNA <- "RNA"
# read in combined sig.summaries file, the script will calculate for every marker-trait association
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/all_RNA.Rdata")
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/all_DNA.Rdata")

# the loaded file is called all_files (just an rbind of sig.summaries results)
cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("./Phenotypes_new.csv", sep=";", header=T)
rownames(pheno)<-pheno$Mouse_Name
pheno$id <- pheno$Mouse_Name

# Genotypes ------------
cat("Reading in genotypes. \n")
geno_in<- BEDMatrix("Cleaning_snps/clean_f2")
rownames(geno_in)<- substr(rownames(geno_in),5 ,19)
rownames(geno_in) <- gsub(x = rownames(geno_in), pattern = "\\/", replacement = ".")  
colnames(geno_in)<- substr(colnames(geno_in),1 ,nchar(colnames(geno_in))-2)

load("Cleaning_snps/clean_snps.Rdata")

# recode from major allele count coding to -1(homo min), 0(hetero), 1(homo ref)
add.matrix <-ifelse(geno_in[,] == 0, 1, ifelse(geno_in[,] == 1, 0,  ifelse(geno_in[,] == 2, -1, NA)))

# Dominance
dom.matrix <- ifelse(abs(geno_in[,]) == 2, 0, ifelse(geno_in[,] == 1, 1,  ifelse(geno_in[,] == 0, 0, NA)))
rownames(dom.matrix)<- rownames(add.matrix)
colnames(dom.matrix)<- colnames(dom.matrix)

out<-data.frame(all_files[2:32],r2_glmm=NA,r2_lr=NA,r2_lm=NA)
# loops over all SNP-trait combo's
for (i in 1:nrow(all_files)){
  # parameters
  trait <- all_files$tax_level[i]
  snp<- all_files$peak.snp[i]
  tx <- all_files$trait[i]
  chr<- all_files$chr[i]
  # taxa abundances
  taxa <-read.csv(paste0("./Phenotyping_27.02.2020/tables_core/",trait,"_table_f2_core.csv"), header=T, row.names = 1)
  taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),]
  names_taxa<- substr(rownames(taxa), 1,nchar(rownames(taxa))-4)
  # Logit Transform 
  taxa <- as.data.frame(sapply(taxa, function(x) inv.logit(x), USE.NAMES = T))
  rownames(taxa)<- names_taxa
  
  da <- merge(taxa[1], pheno, by="row.names")
  individuals <-da$Row.names
  
  geno <- geno_in[individuals,]
  indi_all <- rownames(geno)
  add.mat<- add.matrix[individuals,]
  dom.mat <- dom.matrix[individuals,]
  
  taxa2 <- taxa[tx]
  
  tax <- names(taxa2)
  cat("Running model for ", tax, ".\n")
  names(taxa2)<-"tax"
  
  lmm.data <- merge(taxa2, pheno, by="row.names")
  rownames(lmm.data)<- lmm.data$Row.names
  lmm.data$Row.names<- NULL

  
  individuals <-rownames(lmm.data)
  
  # kinship matrix 
  if (DNAorRNA=="DNA"){
    kinship <- read.table(paste0("Kinship_DNA/kinship_chr",chr,".cXX.txt"))
  }else if (DNAorRNA=="RNA"){
    kinship <- read.table(paste0("./Kinship_RNA/kinship_chr",chr,".cXX.txt"))
  }
  kinship<- as.matrix(kinship)
  rownames(kinship)<- indi_all
  colnames(kinship)<-indi_all
  
  df<-data.frame(lmm.data,ad=add.mat[,snp],dom=dom.mat[,snp],gt=geno[,snp])
  df<-df[is.na(df$ad)==F & is.na(df$tax)==F &is.na(df$gt)==F & is.na(df$dom)==F,]
  null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  model <- relmatLmer(tax ~ ad+dom+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  
  out[i,"r2_glmm"]<-r.squaredGLMM(model,null_model)[1] ### this is the good one
  out[i,"r2_lr"]<-r.squaredLR(model,null = null_model)[1]
  out[i,"r2_lm"]<-summary(lm(tax~ad+dom, df))$adj.r.squared
  
  
  }


save(out, file="pve_snps_rna.Rdata")
  save(out, file="pve_snps_dna.Rdata")
