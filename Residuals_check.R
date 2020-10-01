# This script will 
# - run the model locally
# - test the assumptions of normality
# - calculate the heritability

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



trait <- "body_weight"
covar <- "body_length"

# Phenotypes ----------
cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("Phenotypes_histology_new.csv", sep=";", header=T)
pheno <- pheno[1:327,] # remove the last NA lines
rownames(pheno)<-pheno$mouse_name
pheno$id <- pheno$mouse_name

# Genotypes ------------
cat("Reading in genotypes. \n")
geno<- BEDMatrix("Cleaning_snps/clean_f2")
rownames(geno)<- substr(rownames(geno),4 ,18)
rownames(geno) <- gsub(x = rownames(geno), pattern = "\\/", replacement = ".")  
colnames(geno)<- substr(colnames(geno),1 ,nchar(colnames(geno))-2)
indi_all <- rownames(geno)
load("Cleaning_snps/clean_snps.Rdata")
pheno <- pheno[indi_all, ]
individuals <- as.character(rownames(pheno))
geno <- geno[individuals,]
indi_all <- rownames(geno)

# recode from major allele count coding to -1(homo min), 0(hetero), 1(homo ref)
add.mat <-ifelse(geno[,] == 0, 1, ifelse(geno[,] == 1, 0,  ifelse(geno[,] == 2, -1, NA)))
# add.mat <- geno
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

##### lme4qtl #######
# we focus on chromosome 1
chr <- 1
cat("Running model on chromosome ",chr, ".\n")
kinship <- read.table(paste0("./Kinship_RNA/kinship_chr",chr,".cXX.txt"))


kinship<- as.matrix(kinship)
rownames(kinship)<- indi_all
colnames(kinship)<-indi_all
#kinship <- kinship[individuals,individuals]
marker_chr <- snps[which(snps$chr==chr),1]
gts<-add.mat[individuals,marker_chr]
sub<-colnames(gts)[colMeans(gts,na.rm=T)/2>0.025 & colMeans(gts,na.rm=T)/2<0.975 & is.na(colMeans(gts,na.rm=T))==F]
gts<-gts[,sub]
# sub <- sub[1:100]
out<-data.frame(snps[sub,1:6],tax=NA,n=NA,AA=NA,AB=NA,BB=NA,add.Beta=NA,add.StdErr=NA,add.T=NA, dom.Beta=NA,dom.StdErr=NA, dom.T=NA,P=NA, add.P=NA, dom.P=NA)
out$tax<-tax




df<-data.frame(lmm.data)
null_model <- relmatLmer(tax ~  covar + (1|id), df,relmat=list(id=kinship))
size <- nrow(df)

# heritability
vf <- as.data.frame(VarCorr(null_model))[, c("grp", "vcov")]
vf$prop <- with(vf, vcov / sum(vcov)) # prop gives the heritability estimate

# residuals 
r1 <- residuals(null_model)

# qqplot
qqnorm(r1)
qqline(r1)

# histogram
hist(r1, breaks = 30)

# inference
step(null_model, alpha.random = 1, alpha.fixed = 1) # doesn't work


#confidence interval
prof <- profile(null_model, which = "theta_", prof.scale = "varcov")
# `?lme4qtl::varpropProf`
prof_prop <- varpropProf(prof)
ci <- confint(prof_prop, level = 0.95)
ci

# residuals vs fitted
plot(null_model)

# shapiro
shapiro.test(r1) # when significant, not normally distributed

# Anderson-Darling
library(nortest)
ad.test(r1)

library(moments)
skewness(r1) # should be equal to 0
kurtosis(r1)

# independency on residuals ----
durbinWatsonTest(r1)

# Here starts the loop over all the snps: 

# snp <- sub[1]
# sub <- sub[1:100]
system.time(f<-mclapply(as.list(sub), function(snp){
  df<-data.frame(lmm.data,ad=gts[,snp],dom=dom.mat[,snp],gt=geno[,snp])
  df<-df[is.na(df$ad)==F & is.na(df$tax)==F & is.na(df$covar)==F,]
  if (nrow(df)!=size){
    null_model <- relmatLmer(tax ~  covar +   (1|id), df,relmat=list(id=kinship))
  }
  model <- relmatLmer(tax ~ ad+dom +  covar + (1|id), df,relmat=list(id=kinship))
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


