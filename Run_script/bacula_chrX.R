
##################################################################
##     This script will run the model for the X chromosome.     ##
##            Here we do not have a dominance effect.           ##
##################################################################

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


args <- commandArgs(TRUE)

tx <- as.integer(args[1])
# Parameters ####
#DNAorRNA <- args[2]
#trait <- args[3]


##----------------------------------------------------------------
##                        1. Data import                        --
##----------------------------------------------------------------


cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("./input/BaculaPhenosG2sMap20200921.csv", sep=";", header=T)
rownames(pheno) <- pheno$mouse_name
pheno$id <- pheno$mouse_name

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
marker_chr <- snps[which(snps$chr=="X"),1]

geno <- geno[individuals,marker_chr]
indi_all <- rownames(geno)

# recode from major allele count coding to -1(homo min), 0(hetero), 1(homo ref)
add.mat <-ifelse(geno[,] == 0, 1, ifelse(geno[,] == 2, -1, NA))


taxa2 <- pheno[tx]

tax <- names(taxa2)
cat("Running model for", tax, "on chromosome X.\n")
names(taxa2)<-"tax"

lmm.data <- merge(taxa2, pheno, by="row.names")
rownames(lmm.data)<- lmm.data$Row.names
lmm.data$Row.names<- NULL
individuals <-rownames(lmm.data)


##----------------------------------------------------------------
##                    2. Load kinship matrix                    --
##----------------------------------------------------------------


kinship <- read.table("./kinship_RNA/kinship_chrX.cXX.txt")


kinship<- as.matrix(kinship)
rownames(kinship)<- indi_all
colnames(kinship)<-indi_all

gts<-add.mat
sub<-colnames(gts)[colMeans(gts,na.rm=T)/2>0.025 & colMeans(gts,na.rm=T)/2<0.975 & is.na(colMeans(gts,na.rm=T))==F]
gts<-gts[,sub]
out<-data.frame(snps[sub,1:6],tax=NA,n=NA,AA=NA,AB=NA,BB=NA,add.Beta=NA,add.StdErr=NA,add.T=NA,P=NA)
out$tax<-tax



##----------------------------------------------------------------
##                         3. Run model                         --
##----------------------------------------------------------------


#snp <- sub[1]
df<-data.frame(lmm.data)
size <- nrow(df)
null_model <- relmatLmer(tax ~  (1|id), df,relmat=list(id=kinship))

system.time(f<-mclapply(as.list(sub), function(snp){
  df<-data.frame(lmm.data,ad=gts[,snp],gt=geno[,snp])
  df<-df<-df[is.na(df$ad)==F & is.na(df$tax)==F &is.na(df$gt)==F,]
  #solves the NA problem
  if (nrow(df)!=size){
    null_model <- relmatLmer(tax ~  (1|id), df,relmat=list(id=kinship))
  }
  model <- relmatLmer(tax ~ ad + (1|id), df,relmat=list(id=kinship))
  res<-c(nrow(df),table(factor(df$ad,levels=c(-1,0,1))),tryCatch(summary(model)$coefficients[2,], error=function(x) return(rep(NA,3))),tryCatch(anova(null_model,model)[2,8], error=function(x) return(NA)))
  names(res)<-c("n","AA","AB","BB","add.Beta","add.StdErr","add.T", "P")
  
  return(res)},mc.cores=getOption("mc.cores", 20)))
out[,c("n","AA","AB","BB","add.Beta","add.StdErr","add.T","P")]<-data.frame(do.call(rbind, f))

# Add a column with the marker index.
n      <- nrow(out)
out <- cbind(out,index = 1:n)
dir.create(path=paste0("./out/bacula/",tax, "/"), showWarnings = F )
saveRDS(out,paste0("./out/bacula/",tax,"/",tax, "_chrX.rds"))



