setwd("~/Documents/research/Experiments/Final_QTL_mapping/")

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
library(tidyverse)

# Parameters ####
DNAorRNA <- "RNA"
# read in combined sig.summaries file, the script will calculate for every marker-trait association
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/all_RNA.Rdata")
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



markers_by_trait<-all_files %>% filter(P.type=="P") %>% group_by(trait) %>%
  summarise(all_markers = paste0(markers, collapse = ", "), tax_level=tax_level) %>% mutate(uniq_markers=NA, num_uniq_markers=NA) %>% distinct()
for (i in 1:nrow(markers_by_trait)){
  mm <- unique(str_split(markers_by_trait$all_markers[i], ", "))
  markers_by_trait$uniq_markers[i] <- paste0(mm[[1]], collapse = ", ")
  markers_by_trait$num_uniq_markers[i] <- length(mm[[1]])
  
}


ov_tab<-table(markers_by_trait$num_uniq_markers)

# with 2 snps 
#snp2 <- markers_by_trait %>% filter(num_uniq_markers==2) %>% ungroup()
out<-data.frame(markers_by_trait,r2_glmm=NA,r2_lr=NA,r2_lm=NA)

# loops over all SNP-trait combo's
# for (i in 1:nrow(markers_by_trait)){
  for (i in which(is.na(out$r2_glmm))){
  # parameters
  trait <- markers_by_trait$tax_level[i]
  snps<- str_split(markers_by_trait$uniq_markers[i], ", ")[[1]]
  tx <- markers_by_trait$trait[i]
  # chr<- snp2$chr[i]
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
    kinship <- read.table(paste0("Results/snp_heritability/kinship_DNA.cXX.txt"))
  }else if (DNAorRNA=="RNA"){
    kinship <- read.table(paste0("Results/snp_heritability/kinship_RNA.cXX.txt"))
  }
  kinship<- as.matrix(kinship)
  rownames(kinship)<- indi_all
  colnames(kinship)<-indi_all
  
  
  
  if (length(snps)==1){
    df<-data.frame(lmm.data,ad=add.mat[,snps[1]],dom=dom.mat[,snps[1]],gt=geno[,snps[1]])
    df<-df[is.na(df$ad)==F & is.na(df$tax)==F &is.na(df$gt)==F & is.na(df$dom)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad+dom+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  } else if (length(snps)==2){
    df<-data.frame(lmm.data,ad1=add.mat[,snps[1]],ad2=add.mat[,snps[2]], dom1=dom.mat[,snps[1]], dom2=dom.mat[,snps[2]],gt1=geno[,snps[1]], gt2=geno[,snps[2]])
    df<-df[is.na(df$ad1)==F & is.na(df$ad2)==F & is.na(df$tax)==F & is.na(df$gt1)==F &is.na(df$gt2)==F& is.na(df$dom1)==F & is.na(df$dom2)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+dom1+dom2+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  } else if (length(snps)==3){
    df<-data.frame(lmm.data,ad1=add.mat[,snps[1]],ad2=add.mat[,snps[2]], ad3=add.mat[,snps[3]],
                   dom1=dom.mat[,snps[1]], dom2=dom.mat[,snps[2]],dom3=dom.mat[,snps[3]],
                   gt1=geno[,snps[1]], gt2=geno[,snps[2]], gt3=geno[,snps[3]])
    df<-df[is.na(df$ad1)==F & is.na(df$ad2)==F & is.na(df$ad3)==F & is.na(df$tax)==F & is.na(df$gt1)==F &
             is.na(df$gt2)==F& is.na(df$dom1)==F & is.na(df$dom2)==F & is.na(df$dom3)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+dom1+dom2+dom3+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  } else if (length(snps)==4){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:4)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:4)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+dom1+dom2+dom3+dom4+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  } else if (length(snps)==5){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:5)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:5)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+dom1+dom2+dom3+dom4+dom5+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  } else if (length(snps)==6){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:6)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:6)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+dom1+dom2+dom3+dom4+dom5+dom6+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  } else if (length(snps)==7){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:7)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:7)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+dom1+dom2+dom3+dom4+dom5+dom6+dom7+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  } else if (length(snps)==8){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:8)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:8)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))    
  } else if (length(snps)==9){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:9)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:9)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))   
  } else if (length(snps)==10){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:10)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:10)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))   
  } else if (length(snps)==11){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:11)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:11)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))   
  } else if (length(snps)==12){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:12)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:12)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))   
  } else if (length(snps)==13){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:13)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:13)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))   
  } else if (length(snps)==14){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:14)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:14)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))   
  } else if (length(snps)==17){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:17)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:17)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))   
  } else if (length(snps)==18){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:18)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:18)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))   
  } else if (length(snps)==19){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:19)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:19)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))  
  } else if (length(snps)==20){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:20)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:20)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))  
  } else if (length(snps)==22){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:22)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:22)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))  
  } else if (length(snps)==23){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:23)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:23)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))  
  } else if (length(snps)==24){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:24)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:24)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))    
  } else if (length(snps)==25){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:25)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:25)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+dom25+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==26){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:26)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:26)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+dom25+dom26+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==28){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:28)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:28)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==29){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:29)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:29)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==30){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:30)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:30)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==31){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:31)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:31)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==34){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:34)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:34)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==37){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:37)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:37)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==45){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:45)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:45)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==55){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:55)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:55)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==56){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:56)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:56)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+ad56+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+dom56+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==58){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:58)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:58)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+ad56+
                          ad57+ad58+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+dom56+dom57+dom58+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==63){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:63)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:63)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+ad56+
                          ad57+ad58+ad59+ad60+ad61+ad62+ad63+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+dom56+dom57+dom58+dom59+dom60+dom61+dom62+dom63+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==68){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:68)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:68)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+ad56+
                          ad57+ad58+ad59+ad60+ad61+ad62+ad63+ad64+ad65+ad66+ad67+ad68+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+dom56+dom57+dom58+dom59+dom60+dom61+dom62+dom63+dom64+dom65+dom66+dom68+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==74){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:74)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:74)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+ad56+
                          ad57+ad58+ad59+ad60+ad61+ad62+ad63+ad64+ad65+ad66+ad67+ad68+ad69+ad70+ad71+ad72+ad73+ad74+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+dom56+dom57+dom58+dom59+dom60+dom61+dom62+dom63+dom64+dom65+dom66+dom68+dom69+dom70+
                          dom71+dom72+dom73+dom74+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==77){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:77)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:77)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+ad56+
                          ad57+ad58+ad59+ad60+ad61+ad62+ad63+ad64+ad65+ad66+ad67+ad68+ad69+ad70+ad71+ad72+ad73+ad74+ad75+ad76+ad77+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+dom56+dom57+dom58+dom59+dom60+dom61+dom62+dom63+dom64+dom65+dom66+dom68+dom69+dom70+
                          dom71+dom72+dom73+dom74+dom75+dom76+dom77+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==86){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:86)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:86)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+ad56+
                          ad57+ad58+ad59+ad60+ad61+ad62+ad63+ad64+ad65+ad66+ad67+ad68+ad69+ad70+ad71+ad72+ad73+ad74+ad75+ad76+ad77+ad79+ad80+ad81+ad82+ad83+ad84+ad85+ad86+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+dom56+dom57+dom58+dom59+dom60+dom61+dom62+dom63+dom64+dom65+dom66+dom68+dom69+dom70+
                          dom71+dom72+dom73+dom74+dom75+dom76+dom77+dom78+dom79+dom80+dom81+dom82+dom83+dom84+dom85+dom86+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else if (length(snps)==111){
    addi <- add.mat[,snps]
    colnames(addi)<- paste0("ad", 1:111)
    domi <- dom.mat[,snps]
    colnames(domi)<- paste0("dom", 1:111)
    df<-data.frame(lmm.data,addi, domi)
    df<-df[is.na(df$tax)==F,]
    null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
    model <- relmatLmer(tax ~ ad1+ad2+ad3+ad4+ad5+ad6+ad7+ad8+ad9+ad10+ad11+ad12+ad13+ad14+ad15+ad16+ad17+ad18+ad19+ad20+ad21+ad22+ad23+ad24+ad25+ad26+ad27+ad28+ad29+
                          ad30+ad31+ad32+ad33+ad34+ad35+ad36+ad37+ad38+ad39+ad40+ad41+ad42+ad43+ad44+ad45+ad46+ad47+ad48+ad49+ad50+ad51+ad52+ad53+ad54+ad55+ad56+
                          ad57+ad58+ad59+ad60+ad61+ad62+ad63+ad64+ad65+ad66+ad67+ad68+ad69+ad70+ad71+ad72+ad73+ad74+ad75+ad76+ad77+ad79+ad80+ad81+ad82+ad83+ad84+ad85+ad86+
                          ad87+ad88+ad89+ad90+ad91+ad92+ad93+ad94+ad95+ad96+ad97+ad98+ad99+ad100+ad101+ad102+ad103+ad104+ad105+ad106+ad107+ad108+ad109+ad110+ad111+
                          dom1+dom2+dom3+dom4+dom5+dom6+dom7+dom8+dom9+dom10+dom11+dom12+dom13+dom14+dom15+dom16+dom17+dom18+dom19+dom20+dom21+dom22+dom23+dom24+
                          dom25+dom26+dom27+dom28+dom29+dom30+dom31+dom32+dom33+dom34+dom35+dom36+dom37+dom38+dom39+dom40+dom41+dom42+dom43+dom44+dom45+dom46+
                          dom47+dom48+dom49+dom50+dom51+dom52+dom53+dom54+dom55+dom56+dom57+dom58+dom59+dom60+dom61+dom62+dom63+dom64+dom65+dom66+dom68+dom69+dom70+
                          dom71+dom72+dom73+dom74+dom75+dom76+dom77+dom78+dom79+dom80+dom81+dom82+dom83+dom84+dom85+dom86+dom87+dom88+dom89+dom90+dom91+dom92+dom93+dom94+
                          dom95+dom95+dom96+dom97+dom98+dom99+dom100+dom101+dom102+dom103+dom104+dom105+dom106+dom107+dom108+dom109+dom110+dom111+
                          (1|mating.pair)+ (1|id), df,relmat=list(id=kinship)) 
  } else {
    next 
  }
    out[i,"r2_glmm"]<-r.squaredGLMM(model,null_model)[1] ### this is the good one
  # out[i,"r2_lr"]<-r.squaredLR(model,null = null_model)[1]
  # out[i,"r2_lm"]<-summary(lm(tax~ad1+ad2+dom1+dom2, df))$adj.r.squared
  
  
}


save(out, file="Results/Bacterial traits/Revisions/pve_all_snps_rna.Rdata")
# save(out, file="pve_snps_dna.Rdata")

out_rna <- out
out_dna <- out

load("pve_snps_rna.Rdata")
load("pve_snps_dna.Rdata")
load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/otu/SV17_gwscan_allchr.Rdata")
# to get the counts 
input<- all_files[,1:23] %>% 
  left_join(SV17_gwscan[,1:12], by=c("peak.snp"="marker"))

input$MAF<- ((2*input$AA/input$n) +input$AB/input$n)/2

# additive effect ####
input$add.num <- 2*(input$add.Beta_peak_snp^2)*input$MAF * (1-input$MAF)
input$add.denom <- 2*(input$add.Beta_peak_snp^2)*input$MAF * (1-input$MAF) + (input$add.StdErr_peak_snp^2)* 2 * input$n *input$MAF * (1-input$MAF)

input$add.pve.lead.snp <- input$add.num/input$add.denom

input$dom.num <- 2*(input$dom.Beta_peak_snp^2)*input$MAF * (1-input$MAF)
input$dom.denom <- 2*(input$dom.Beta_peak_snp^2)*input$MAF * (1-input$MAF) + (input$dom.StdErr_peak_snp^2)* 2 * input$n *input$MAF * (1-input$MAF)

input$dom.pve.lead.snp <- input$dom.num/input$dom.denom
all<- input %>% inner_join(out, by="peak.snp")

library(ggpmisc)
ggplot(all, aes(x=add.pve.lead.snp, y=r2_glmm))+geom_point() +geom_smooth(method = "lm")+
  stat_poly_eq( formula=y ~ x,
                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                parse = TRUE)
sum_pve_rna <- out_rna %>% 
  group_by(trait, P.type) %>%
  dplyr::summarise(tot=sum(r2_glmm), n=n())
sum_pve_dna <- out_dna %>% 
  group_by(trait, P.type) %>%
  dplyr::summarise(tot=sum(r2_glmm), n=n())

##### heritability #####
herit <- read_delim("Results//snp_heritability/lme4qtl/sign_heritability_rna.csv", delim=";")
library(ggrepel)
herit.pve.rna <- herit %>% 
  inner_join(sum_pve_rna, by=c("taxa"="trait"))
herit.pve.rna %>% 
  filter(P.type=="P") %>% 
  ggplot(aes(x=h2_id, y=tot)) + geom_point() +geom_smooth(method="lm" ) +geom_text_repel(aes(label=taxa))+
  stat_poly_eq( formula=y ~ x, aes(label = paste(..eq.label.., ..rr.label.., p.value.label, sep = "~~~")), 
                parse = TRUE)+theme_pubr()
ggsave("Results/snp_heritability/lme4qtl/pve_vs_herit_rna.pdf")


herit <- read_delim("Results//snp_heritability/lme4qtl/sign_heritability_dna.csv", delim=";")
library(ggrepel)
library(ggsci)
herit.pve.dna <- herit %>% 
  inner_join(sum_pve_dna, by=c("taxa"="trait"))
herit.pve.dna %>% 
  filter(P.type=="P") %>% 
  ggplot(aes(x=h2_id, y=tot)) + geom_point() +geom_smooth(method="lm" ) +geom_text_repel(aes(label=taxa))+
  stat_poly_eq( formula=y ~ x, aes(label = paste(..eq.label.., ..rr.label.., p.value.label, sep = "~~~")), 
                parse = TRUE)+theme_pubr()
ggsave("Results/snp_heritability/lme4qtl/pve_vs_herit_dna.pdf")

herit.pve.dna$dna.rna<- "DNA"
herit.pve.rna$dna.rna<- "RNA"
herit.pve.dna$lrt<- NULL
herit.pve.dna$Sample_coef_variation <- NULL
herit.pve <- rbind(herit.pve.dna, herit.pve.rna)

out_all <- rbind(out_dna, out_rna)
uniq_all<- out_all %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
uniq_peak_snp_all <- out_all %>% 
  distinct(peak.snp)


load("pve_snps_dna.Rdata")
out_dna<- out
load("pve_snps_rna.Rdata")
out_rna<-out
out_rna$dna.rna <- "RNA"
out_dna$dna.rna <- "DNA"
out <- rbind(out_dna, out_rna)
ggplot(out)+geom_histogram(aes(x=r2_glmm), fill="white", color="black")

tot_plot<-ggplot(out, aes(r2_glmm, fill = dna.rna)) +
  geom_histogram(bins=50,aes(y = stat(count) / sum(count)), color="black")+
  scale_y_continuous(labels = scales::percent,breaks=seq(0, 0.2, by=0.025)) +
  scale_x_continuous(labels = scales::percent,breaks=seq(0, 0.7, by=0.1)) +
  labs(x="PVE by lead SNP per significant locus", y="Frequency", fill="") +scale_fill_d3()+facet_wrap(~dna.rna)+
  theme_pubr()+theme(legend.position = "none")
tot_plot
ggsave("Results/Bacterial traits/PVE_snps/PVE_per_locus_wrap.pdf")

sum_pve <- out %>% 
  group_by(trait, P.type,dna.rna) %>%
  dplyr::summarise(tot=sum(r2_glmm), n=n())
sum_plot<-ggplot(sum_pve, aes(tot, fill = dna.rna)) +
  geom_histogram(bins=50,aes(y = stat(count) / sum(count)), color="black")+
  scale_y_continuous(labels = scales::percent, breaks=seq(0, 0.13, by=0.02)) +
  scale_x_continuous(labels = scales::percent, breaks=seq(0, 2.6, by=0.25)) +
  labs(x="Total PVE by lead SNPs of significant loci per taxon and P type", y="Frequency", fill="") +scale_fill_d3()+
  facet_wrap(.~dna.rna)+theme_pubr()+theme(legend.position = "none")
sum_plot
ggsave("Results/Bacterial traits/PVE_snps/sum_PVE_wrap.pdf")
library(patchwork)
tot_plot/sum_plot+plot_annotation(tag_levels = "A")
ggsave("Results/Bacterial traits/PVE_snps/sum_and_count_PVE_wrap.pdf")


herit.pve <- herit.pve %>% mutate(larger=ifelse(tot>h2_id, 1, 0))
sum(herit.pve$larger)
herit.pve %>%  filter(P.type=="P") %>%  summarise(sum(larger))
herit.pve %>%  filter(P.type=="add.P") %>%  summarise(sum(larger))
herit.pve %>%  filter(P.type=="dom.P") %>%  summarise(sum(larger))
herit.pve %>% distinct(dna.rna, taxa) %>% count()


# ggplot(sum_pve, aes(tot)) +
#   geom_histogram(bins=50,aes(y = stat(count) / sum(count)), color="black")+
#   scale_y_continuous(labels = scales::percent) +
#   scale_x_continuous(labels = scales::percent) + facet_wrap(~dna.rna)+
#   labs(x="PVE by lead SNP per significant locus", y="Frequency", fill="") +scale_fill_d3()
