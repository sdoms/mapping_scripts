setwd("/home/doms/glm")
.libPaths("/home/doms/R/x86_64-pc-linux-gnu-library/3.5")

library(lme4qtl)
library(BEDMatrix)
library(MASS)
library(parallel)
library(gtools)
library(lme4)
library(plyr)
library(lmerTest)
library(car)

#library(qqman)

args <- commandArgs(TRUE)

#tx <- as.integer(args[1])
# Parameters ####
#DNAorRNA <- args[2]
link_family <- args[1]
trait <- args[2]
covar <- args[3]

#source("kinship/make_gemma_kinship.R")

# Reading in data ####
# Phenotypes ----------
cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("./input/Phenotypes_histology_new.csv", sep=";", header=T, stringsAsFactors = F)
pheno <- pheno[1:327,]
rownames(pheno)<-pheno$mouse_name
pheno$id <- pheno$mouse_name
#individuals <- as.character(rownames(pheno))

# Genotypes ------------
cat("Reading in genotypes. \n")
geno<- BEDMatrix("input/clean_f2")
rownames(geno)<- substr(rownames(geno),4 ,18)
rownames(geno) <- gsub(x = rownames(geno), pattern = "\\/", replacement = ".")  
colnames(geno)<- substr(colnames(geno),1 ,nchar(colnames(geno))-2)
indi_all <- rownames(geno)
load("input/clean_snps.Rdata")
pheno <- pheno[indi_all, ]
individuals <- as.character(rownames(pheno))
geno <- geno[individuals,]
indi_all <- rownames(geno)
# make kinship matrices
#writeKinship("input/clean_f2", pheno,"kinship/", individuals)



# recode from major allele count coding to -1(homo min), 0(hetero), 1(homo ref)
add.mat <-ifelse(geno[,] == 0, 1, ifelse(geno[,] == 1, 0,  ifelse(geno[,] == 2, -1, NA)))

# Dominance
dom.mat <- ifelse(abs(geno[,]) == 2, 0, ifelse(geno[,] == 1, 1,  ifelse(geno[,] == 0, 0, NA)))
rownames(dom.mat)<- rownames(add.mat)
colnames(dom.mat)<- colnames(dom.mat)


taxa2 <- pheno[,c(trait, covar)]


cat("Running model for ", trait, " with ", covar, " as covariate.\n")
names(taxa2)<-c("tax", "covar")
tax <- trait 
lmm.data <- merge(taxa2, pheno, by="row.names")
rownames(lmm.data)<- lmm.data$Row.names
lmm.data$Row.names<- NULL


individuals <-rownames(lmm.data)


# model
gwasResults <- data.frame()
for (chr in 1:19){
  cat("Running model on chromosome",chr,".\n")
  kinship <- read.table(paste0("./kinship_RNA/kinship_chr",chr,".cXX.txt"))
  
  
  kinship<- as.matrix(kinship)
  rownames(kinship)<- indi_all
  colnames(kinship)<-indi_all
  #kinship <- kinship[individuals,individuals]
  marker_chr <- snps[which(snps$chr==chr),1]
  gts<-add.mat[individuals,marker_chr]
  sub<-colnames(gts)[colMeans(gts,na.rm=T)/2>0.025 & colMeans(gts,na.rm=T)/2<0.975 & is.na(colMeans(gts,na.rm=T))==F]
  gts<-gts[,sub]
  out<-data.frame(snps[sub,1:6],tax=NA,n=NA,AA=NA,AB=NA,BB=NA,add.Beta=NA,add.StdErr=NA,add.T=NA, dom.Beta=NA,dom.StdErr=NA, dom.T=NA,P=NA, add.P=NA, dom.P=NA)
  out$tax<-tax
  
  
  
  #snp <- sub[1]
  df<-data.frame(lmm.data)
  #null_model <- relmatLmer(tax ~  covar + (1|id), df,relmat=list(id=kinship))
  null_model <- relmatGlmer(tax ~ covar + (1|id), df,relmat=list(id=kinship),family=link_family)
  size <- nrow(df)
  
  
  system.time(f<-mclapply(as.list(sub), function(snp){
    df<-data.frame(lmm.data,ad=gts[,snp],dom=dom.mat[,snp],gt=geno[,snp])
    df<-df[is.na(df$ad)==F & is.na(df$tax)==F & is.na(df$covar)==F,]
    if (nrow(df)!=size){
      null_model <- relmatGlmer(tax ~ covar + (1|id), df,relmat=list(id=kinship),family=link_family)
    }

    model <- relmatGlmer(tax ~ ad+dom +covar+ (1|id), df,relmat=list(id=kinship),family=link_family)
    res<-c(nrow(df),table(factor(df$ad,levels=c(-1,0,1))),tryCatch(summary(model)$coefficients[2,], error=function(x) return(rep(NA,3))),tryCatch(summary(model)$coefficients[3,], error=function(x) return(rep(NA,3))),tryCatch(anova(null_model,model)[2,8], error=function(x) return(NA)), tryCatch(Anova(model)[1:2,3],error=function(x) return(rep(NA, 2)) ))
    names(res)<-c("n","AA","AB","BB","add.Beta","add.StdErr","add.T","dom.Beta", "dom.StdErr", "dom.T", "P", "add.P", "dom.P")
    
    return(res)},mc.cores=getOption("mc.cores", 20)))
  out[,c("n","AA","AB","BB","add.Beta","add.StdErr","add.T","dom.Beta", "dom.StdErr","dom.T","P","add.P", "dom.P")]<-data.frame(do.call(rbind, f))
  
  # Add a column with the marker index.
  n      <- nrow(out)
  out <- cbind(out,index = 1:n)
  #out$log10p <- -1 * log10(out$P)
  gwasResults <- rbind(gwasResults, out)
  dir.create(path=paste0("./out/",tax, "/"), showWarnings = F )
  
  saveRDS(out,paste0("./out/",tax, "/",tax,"_with_covar_", covar,  "_chr_",chr,"with_add_dom.rds"))
  head(out)
  # }
  
  
  
}
saveRDS(out,paste0("./out/",tax, "/",tax,"_with_covar_", covar, "_with_add_dom.rds"))






