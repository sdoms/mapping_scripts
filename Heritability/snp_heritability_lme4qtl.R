setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/snp_heritability/lme4qtl/")

library(lme4qtl)
library(BEDMatrix)
library(MASS)
library(parallel)
library(gtools)
library(lme4)
library(plyr)
library(lmerTest)
library(car)
library(EnvStats)
library(RLRsim)

source("~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/function_for_gemma.r")
gemma <-"/usr/local/bin/gemma"
DNAorRNA="DNA"

for (DNAorRNA in c("DNA", "RNA")){
##### Phenotypes ####
pheno<-read.csv("../../../Phenotypes_histology_new.csv", sep=";", header=T)
rownames(pheno)<-pheno$Mouse_Name
pheno$id <- pheno$Mouse_Name

#### Genotypes ####
geno<- BEDMatrix("../../../Cleaning_snps/clean_f2")
rownames(geno)<- substr(rownames(geno),4 ,18)
rownames(geno) <- gsub(x = rownames(geno), pattern = "\\/", replacement = ".")  
colnames(geno)<- substr(colnames(geno),1 ,nchar(colnames(geno))-2)
indi_all <- rownames(geno)
load("../../../Cleaning_snps/clean_snps.Rdata")
pheno <- pheno[indi_all, ]
individuals <- as.character(rownames(pheno))
geno <- geno[individuals,]
indi_all <- rownames(geno)

##### taxa ####

traits <- c("otu", "genus", "family", "order", "class", "phylum")
for (trait in traits){
taxa <- read.csv(paste0("../../../Phenotyping_27.02.2020/tables_core/",trait, "_table_f2_core.csv"), header=T, row.names = 1)
taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),]
names_taxa<- substr(rownames(taxa), 1,nchar(rownames(taxa))-4)
taxa <- as.data.frame(sapply(taxa, function(x) inv.logit(x), USE.NAMES = T))
rownames(taxa)<- names_taxa
da <- merge(taxa[1], pheno, by="row.names")
individuals <-da$Row.names

##### kinship matrix ######
# genotype kinship
library(argyle)
F2 <- read.plink("../../../Cleaning_snps/clean_f2")
F2 <- recode(F2, "relative")

colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")  

#Select the individuals from the genotype file 
F2 <- as.matrix(as.data.frame(F2))
ind <- rownames(taxa)

individuals <- ncol(F2)
F2_part <- apply(F2[,7:individuals], 2, as.numeric)# make sure that the genotype information(0,1,2) is numeric

F2_part <- F2_part[,colnames(F2_part) %in% ind]

F2 <- cbind(snps[,c(1,5,6)],F2_part)

write.table(F2,"geno.txt",sep = " ",quote = FALSE,row.names = FALSE,
            col.names = FALSE)

write.table(data.frame(rownames(taxa),1),"pheno.txt", quote=F, row.names=F, col.names=F)

# # compute kinship matrix

system(paste0(gemma, " -g geno.txt -gk 1 -p pheno.txt -debug -o kinship_", DNAorRNA  ))




results <- data.frame()
for (tx in 1:ncol(taxa)){
taxa2 <- taxa[tx]
tax <- names(taxa2)
cat("Running model for ", tax, ".\n")
names(taxa2)<-"tax"
lmm.data <- merge(taxa2, pheno, by="row.names")
rownames(lmm.data)<- lmm.data$Row.names
lmm.data$Row.names<- NULL
individuals <-rownames(lmm.data)

kinship <- read.table(paste0("./output/kinship_", DNAorRNA,".cXX.txt"))
kinship<- as.matrix(kinship)
rownames(kinship)<- colnames(F2)[4:ncol(F2)]
colnames(kinship)<- colnames(F2)[4:ncol(F2)]

#snp <- sub[1]
df<-data.frame(lmm.data)

null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
# r1<- residuals(null_model)
# qqnorm(r1)
# qqline(r1)
# hist(r1, breaks = 30)

# heritability
h2 <- VarProp(null_model)

# null_model_reduced <- update(null_model, . ~ . - (1|mating.pair))
m1_null <- update(null_model, . ~ . - (1|id))

rlrt_h2 <- exactRLRT(
  mA = null_model, # the full model under the alternative
  m0 = m1_null, # the model under the null
  seed = 1
)
rlrt_h2

lrt_h2 <- anova(m1_null, null_model)
lrt_h2
scv <- cv(taxa2$tax)

results[tax,"h2_id"]<- h2[1,6]
results[tax,"h2_mating_pair"]<- h2[2,6]
results[tax,"h2_residual"]<- h2[3,6]
results[tax,"rlrt"]<- rlrt_h2$p.value
results[tax,"lrt"]<- lrt_h2$`Pr(>Chisq)`[2]
results[tax,"Sample_coef_variation"]<- scv
}
write.table(results, paste0(trait,"_",DNAorRNA,"_snp_heritabilitylme4_sample_coef_var.csv"), sep=";")
}
}

col_pal<- c("light grey","#8F2D56", "#FBB13C","#218380" )

lme4qtl_dna %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_DNA.pdf", height=15)


lme4qtl_rna %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("RNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_RNA.pdf", height=15)

