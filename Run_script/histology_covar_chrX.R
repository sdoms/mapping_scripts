
###############################################################################
##  This script will run the model with one covariate for the X chromosome.  ##
##                  Here we do not have a dominance effect.                  ##
###############################################################################

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
#library(qqman)

args <- commandArgs(TRUE)


trait <- args[1]
covar <- args[2]
if (length(args)>2){
  data_trans <- args[3]
} else {
  data_trans <- FALSE
}


##----------------------------------------------------------------
##                        1. Data import                        --
##----------------------------------------------------------------


# Phenotypes ----------
cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("input/Phenotypes_histology_new.csv", sep=";", header=T)
pheno <- pheno[1:327,]
rownames(pheno)<-pheno$Mouse_Name
pheno$id <- pheno$Mouse_Name

#individuals <- as.character(pheno$id)
colnames(pheno)[colnames(pheno)==trait] <- "trait"
colnames(pheno)[colnames(pheno)==covar] <- "covar"
# transform counts
if (data_trans){
  pheno$trait <-  inv.logit(pheno$trait)

}

# Genotypes ------------
cat("Reading in genotypes. \n")
geno<- BEDMatrix("input/clean_f2")
rownames(geno)<- substr(rownames(geno),4 ,18)
rownames(geno) <- gsub(x = rownames(geno), pattern = "\\/", replacement = ".")  
colnames(geno)<- substr(colnames(geno),1 ,nchar(colnames(geno))-2)

load("input/clean_snps.Rdata")
marker_chr <- snps[which(snps$chr=="X"),1]
indi_all <- rownames(geno)

pheno <- pheno[indi_all, ]
individuals <- as.character(rownames(pheno))
geno <- geno[individuals,marker_chr]

indi_all <- rownames(geno)

# recode from major allele count coding to -1(homo min), 0(hetero), 1(homo ref)
add.mat <-ifelse(geno[,] == 0, 1, ifelse(geno[,] == 2, -1, NA))



##----------------------------------------------------------------
##                         2. Run model                         --
##----------------------------------------------------------------



cat("Running model for ", trait, ".\n")
  # model
  cat("Running model on chromosome X.\n")
  kinship <- read.table("./kinship_RNA/kinship_chrX.cXX.txt")
  kinship<- as.matrix(kinship)
  rownames(kinship)<- indi_all
  colnames(kinship)<-indi_all
  #kinship <- kinship[individuals,individuals]

  gts<-add.mat[individuals,]
  sub<-colnames(gts)[colMeans(gts,na.rm=T)/2>0.025 & colMeans(gts,na.rm=T)/2<0.975 & is.na(colMeans(gts,na.rm=T))==F]
  gts<-gts[,sub]
  out<-data.frame(snps[sub,1:6],trait=NA,n=NA,AA=NA,AB=NA,BB=NA,add.Beta=NA,add.StdErr=NA,add.T=NA,P=NA)
  out$trait<-trait
  
  
  
  #snp <- sub[1]
  df<-data.frame(pheno)
  size <- nrow(df)
  null_model <- relmatLmer(trait ~ covar+ (1|id), df,relmat=list(id=kinship))
  
  system.time(f<-mclapply(as.list(sub), function(snp){
    df<-data.frame(pheno,ad=gts[,snp],gt=geno[,snp])
    df<-df[is.na(df$ad)==F & is.na(df$trait)==F,]
    #solves the NA problem
    if (nrow(df)!=size){
      null_model <- relmatLmer(trait ~ covar + (1|id), df,relmat=list(id=kinship))
    }
    model <- relmatLmer(trait ~ ad + covar + (1|id), df,relmat=list(id=kinship))
    res<-c(nrow(df),table(factor(df$ad,levels=c(-1,0,1))),tryCatch(summary(model)$coefficients[2,], error=function(x) return(rep(NA,3))),tryCatch(anova(null_model,model)[2,8], error=function(x) return(NA)))
    names(res)<-c("n","AA","AB","BB","add.Beta","add.StdErr","add.T", "P")
    
    return(res)},mc.cores=getOption("mc.cores", 20)))
  out[,c("n","AA","AB","BB","add.Beta","add.StdErr","add.T","P")]<-data.frame(do.call(rbind, f))
  
  # Add a column with the marker index.
  n      <- nrow(out)
  out <- cbind(out,index = 1:n)
  
  saveRDS(out,paste0("./out/",trait,"/", trait, "_with_covar_", covar, "_chrX.rds"))
  


#}
