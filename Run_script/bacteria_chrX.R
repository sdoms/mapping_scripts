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
DNAorRNA <- args[2]
trait <- args[3]

# Phenotypes ----------
cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("./input/Phenotypes_new.csv", sep=";", header=T)
rownames(pheno)<-pheno$Mouse_Name
pheno$id <- pheno$Mouse_Name

taxa <-read.csv(paste0("./input/",trait,"_table_f2_core.csv"), header=T, row.names = 1)
taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),]
names_taxa<- substr(rownames(taxa), 1,nchar(rownames(taxa))-4)
# Logit Transform 
taxa <- as.data.frame(sapply(taxa, function(x) inv.logit(x), USE.NAMES = T))
rownames(taxa)<- names_taxa

da <- merge(taxa[1], pheno, by="row.names")
individuals <-da$Row.names

# Genotypes ------------
cat("Reading in genotypes. \n")
geno<- BEDMatrix("input/clean_f2")
rownames(geno)<- substr(rownames(geno),4 ,18)
rownames(geno) <- gsub(x = rownames(geno), pattern = "\\/", replacement = ".")  
colnames(geno)<- substr(colnames(geno),1 ,nchar(colnames(geno))-2)

load("input/clean_snps.Rdata")

marker_chr <- snps[which(snps$chr=="X"),1]
geno <- geno[individuals,marker_chr]
indi_all <- rownames(geno)
# recode from major allele count coding to -1(homo min), 0(hetero), 1(homo ref)
add.mat <-ifelse(geno[,] == 0, 1, ifelse(geno[,] == 2, -1, NA))

cols<- c(1:ncol(taxa))

#for (tx in cols){
taxa2 <- taxa[tx]

tax <- names(taxa2)
cat("Running model for", tax, "on chromosome X.\n")
names(taxa2)<-"tax"

lmm.data <- merge(taxa2, pheno, by="row.names")
rownames(lmm.data)<- lmm.data$Row.names
lmm.data$Row.names<- NULL
individuals <-rownames(lmm.data)
if (DNAorRNA=="DNA"){
  kinship <- read.table("./kinship_DNA/kinship_chrX.cXX.txt")
}else if (DNAorRNA=="RNA"){
  kinship <- read.table("./kinship_RNA/kinship_chrX.cXX.txt")
}

kinship<- as.matrix(kinship)
rownames(kinship)<- indi_all
colnames(kinship)<-indi_all

gts<-add.mat[individuals,]
sub<-colnames(gts)[colMeans(gts,na.rm=T)/2>0.025 & colMeans(gts,na.rm=T)/2<0.975 & is.na(colMeans(gts,na.rm=T))==F]
gts<-gts[,sub]
out<-data.frame(snps[sub,1:6],tax=NA,n=NA,AA=NA,AB=NA,BB=NA,add.Beta=NA,add.StdErr=NA,add.T=NA,P=NA)
out$tax<-tax



#snp <- sub[1]
df<-data.frame(lmm.data)
size <- nrow(df)
null_model <- relmatLmer(tax ~  (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))

system.time(f<-mclapply(as.list(sub), function(snp){
  df<-data.frame(lmm.data,ad=gts[,snp],gt=geno[,snp])
  df<-df[is.na(df$ad)==F & is.na(df$tax)==F,]
  #solves the NA problem
  if (nrow(df)!=size){
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  }
  model <- relmatLmer(tax ~ ad +(1|mating.pair)+  (1|id), df,relmat=list(id=kinship))
  res<-c(nrow(df),table(factor(df$ad,levels=c(-1,0,1))),tryCatch(summary(model)$coefficients[2,], error=function(x) return(rep(NA,3))),tryCatch(anova(null_model,model)[2,8], error=function(x) return(NA)))
  names(res)<-c("n","AA","AB","BB","add.Beta","add.StdErr","add.T", "P")
  
  return(res)},mc.cores=getOption("mc.cores", 20)))
out[,c("n","AA","AB","BB","add.Beta","add.StdErr","add.T","P")]<-data.frame(do.call(rbind, f))

# Add a column with the marker index.
n      <- nrow(out)
out <- cbind(out,index = 1:n)

saveRDS(out,paste0("./out/",DNAorRNA,"/",tax, "_chrX.rds"))



