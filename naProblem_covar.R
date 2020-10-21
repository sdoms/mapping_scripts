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



# Parameters ####
args <- commandArgs(TRUE)

# Parameters ####
trait <- args[1]
covar <- args[2]



cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("./input/Phenotypes_histology_new.csv", sep=";", header=T)
rownames(pheno)<-pheno$Mouse_Name
pheno$id <- pheno$Mouse_Name



taxa2 <- pheno[,c(trait, covar)]


cat("Running model for ", trait, " with ", covar, " as covariate.\n")
names(taxa2)<-c("tax", "covar")
tax <- trait 
lmm.data <- merge(taxa2, pheno, by="row.names")
rownames(lmm.data)<- lmm.data$Row.names
lmm.data$Row.names<- NULL
individuals <- rownames(lmm.data)
if (!file.exists(paste0("./out/",tax, "/",tax,"_with_covar_", covar,  "_chr_1with_add_dom.rds"))){cat(tax, "does not exist!\n") 
  next}
  

  

  
  # Genotypes
  geno<- BEDMatrix("input/clean_f2")
  rownames(geno)<- substr(rownames(geno),4 ,18)
  rownames(geno) <- gsub(x = rownames(geno), pattern = "\\/", replacement = ".")  
  colnames(geno)<- substr(colnames(geno),1 ,nchar(colnames(geno))-2)
  
  load("input/clean_snps.Rdata")
  

  indi_all <- rownames(geno)
  lmm.data <- lmm.data[indi_all,]
  individuals <-rownames(lmm.data)
  geno <- geno[individuals,]
  # recode from major allele count coding to -1(homo min), 0(hetero), 1(homo ref)
  add.mat <-ifelse(geno[,] == 0, 1, ifelse(geno[,] == 1, 0,  ifelse(geno[,] == 2, -1, NA)))
  
  # Dominance
  dom.mat <- ifelse(abs(geno[,]) == 2, 0, ifelse(geno[,] == 1, 1,  ifelse(geno[,] == 0, 0, NA)))

  
  gwasResults <- data.frame()
  for (chr in 1:19){
    out <- readRDS(paste0("./out/",tax,"/",tax,"_with_covar_", covar,  "_chr_",chr,"with_add_dom.rds"))
    
      kinship <- read.table(paste0("./kinship_RNA/kinship_chr",chr,".cXX.txt"))
    
  
    kinship<- as.matrix(kinship)
    rownames(kinship)<- indi_all
    colnames(kinship)<-indi_all
    #kinship <- kinship[individuals,individuals]
    marker_chr <- snps[which(snps$chr==chr),1]
    gts<-add.mat[individuals,marker_chr]
    sub<-out[is.na(out$P),"marker"]
    gts<-gts[,sub]
    
    system.time(f<-mclapply(sub, function(snp){
      df<-data.frame(lmm.data,ad=gts[,snp],dom=dom.mat[,snp],gt=geno[,snp])
      df<-df[is.na(df$ad)==F & is.na(df$tax)==F &is.na(df$gt)==F & is.na(df$dom)==F,]
      null_model <- relmatLmer(tax ~ covar+ (1|id), df,relmat=list(id=kinship))
      model <- relmatLmer(tax ~ ad+dom +covar + (1|id), df,relmat=list(id=kinship))
      res<-c(snp, nrow(df),table(factor(df$ad,levels=c(-1,0,1))),tryCatch(summary(model)$coefficients[2,], error=function(x) return(rep(NA,3))),tryCatch(summary(model)$coefficients[3,], error=function(x) return(rep(NA,3))),tryCatch(anova(null_model,model)[2,8], error=function(x) return(NA)), tryCatch(Anova(model)[1:2,3],error=function(x) return(rep(NA, 2)) ))
      
      names(res)<-c("marker","n","AA","AB","BB","add.Beta","add.StdErr","add.T","dom.Beta", "dom.StdErr", "dom.T", "P", "add.P", "dom.P")
      #out[snp,8:18]<-res[2:12]
      return(res)},mc.cores=getOption("mc.cores", 20)))
    
    f2 <-data.frame(do.call(rbind, f),stringsAsFactors = F)
    out[out$marker %in% f2$marker,8:20] <- data.frame(lapply(f2[,2:14],as.numeric))
    gwasResults <- rbind(gwasResults, out)
    saveRDS(out,paste0("out/",tax,"/",tax,"_with_",covar,"_chr_",chr,"_new.rds"))
    head(out)
  }
  saveRDS(gwasResults,paste0("out/",tax,"/",tax,"_with_",covar,"_new.rds"))

