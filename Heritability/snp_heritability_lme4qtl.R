setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/snp_heritability/lme4qtl/")

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

null_model_reduced <- update(null_model, . ~ . - (1|mating.pair))
m1_null <- update(null_model, . ~ . - (1|id))

rlrt_h2 <- exactRLRT(
  null_model_reduced, # the reduced model with only the effect to be tested
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
write.csv2(results, paste0(trait,"_",DNAorRNA,"_snp_heritabilitylme4_sample_coef_var.csv"))
}
}


Sv_gemma <- readxl::read_excel("../cospeciation_heritability_11_2020.xlsx", sheet=1)
sv_lme4qtl_dna <- read.csv2("otu_DNA_snp_heritabilitylme4_sample_coef_var.csv")
sv_lme4qtl_rna <- read.csv2("otu_RNA_snp_heritabilitylme4_sample_coef_var.csv")
library(corrplot)
library(Hmisc)
library(ggrepel)
library(ggpmisc)


sv_lme4qtl <- merge(sv_lme4qtl_dna, sv_lme4qtl_rna, by="X")
colnames(sv_lme4qtl)<- c("SV", "h2.id.dna", "h2.mp.dna", "h2.resid.dna","rlrt.dna", "lrt.dna", "sam.covar.dna", "h2.id.rna", "h2.mp.rna", "h2.resid.rna", "rlrt.rna", "lrt.rna","sam.covar.rna")

sv_all <- merge(sv_lme4qtl, Sv_gemma, by.x="SV", by.y="Taxon")
rcorr(sv_all$h2.id.dna, sv_all$DNA, type="spearman")
rcorr(sv_all$h2.id.rna, sv_all$RNA, type="spearman")

# corrplot(cor.mat)
sv_all$h2.id.dna<-as.numeric(sv_all$h2.id.dna)
sv_all$CospeciationScoreGenus<- as.numeric(sv_all$CospeciationScoreGenus)
sv_dna<-ggplot(sv_all, aes(x=h2.id.dna, y=CospeciationScoreGenus))+geom_point() +   geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), 
                                                  col="pink"),label.x  = 'middle', label.y = 0.9, size = 3)+ 
  geom_text_repel(aes(label = SV), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.9, 
               parse = TRUE, size = 3)+theme_bw() + ggtitle("DNA")
sv_dna
sv_all$h2.id.rna<-as.numeric(sv_all$h2.id.rna)
sv_rna<-ggplot(sv_all, aes(x=h2.id.rna, y=CospeciationScoreGenus))+geom_point() +  
geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", 
                                                                signif(..p.value.., digits = 4), sep = "" ), col="pink"),
                  label.x  = 'middle', label.y = 0.9, size = 3)+ geom_text_repel(aes(label = SV), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.9, parse = TRUE, 
               size = 3)+theme_bw() + ggtitle("RNA")
sv_rna

sv_lme4qtl_dna_m <- sv_lme4qtl_dna %>% 
  mutate(X=forcats::fct_reorder(X,as.numeric(h2_id))) %>% 
  dplyr::select(X, h2_id, h2_mating_pair, h2_residual) %>% 
  pivot_longer(!X,names_to="h2_kind", values_to="h2_estimate") %>% 
  
  ggplot(aes(x=X, y=h2_estimate))+geom_bar(aes(fill=forcats::fct_rev(h2_kind)),position="stack",stat="identity")+coord_flip()+
  scale_fill_d3() +labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")
sv_lme4qtl_dna_m

sv_lme4qtl_rna_m <- sv_lme4qtl_rna %>% 
  mutate(X=forcats::fct_reorder(X,as.numeric(h2_id))) %>% 
  dplyr::select(X, h2_id, h2_mating_pair, h2_residual) %>% 
  pivot_longer(!X,names_to="h2_kind", values_to="h2_estimate") %>% 
  
  ggplot(aes(x=X, y=h2_estimate))+geom_bar(aes(fill=forcats::fct_rev(h2_kind)),position="stack",stat="identity")+coord_flip()+
  scale_fill_d3() +labs(x="", y="Heritability estimate", fill="")+ggtitle("RNA")
sv_lme4qtl_rna_m

(sv_lme4qtl_dna_m + sv_lme4qtl_rna_m +plot_layout(guides="collect"))/(sv_dna+sv_rna)+plot_annotation(tag_levels = "A")
ggsave("heritability_lme4qtl_sv.pdf", width=18,height=20)
#### Genus ####
genus_lme4qtl_dna <- read.csv2("genus_DNA_snp_heritabilitylme4_sample_coef_var.csv")
genus_lme4qtl_rna <- read.csv2("genus_RNA_snp_heritabilitylme4_sample_coef_var.csv")
genus_lme4qtl <- merge(genus_lme4qtl_dna, genus_lme4qtl_rna, by="X")
colnames(genus_lme4qtl)<- c("Genus", "h2.id.dna", "h2.mp.dna", "h2.resid.dna","rlrt.dna", "lrt.dna", "sam.covar.dna", "h2.id.rna", "h2.mp.rna", "h2.resid.rna", "rlrt.rna", "lrt.rna","sam.covar.rna")
genus_gemma <- readxl::read_excel("../cospeciation_heritability_11_2020.xlsx", sheet="Genus")
genus_all <- merge(genus_lme4qtl, genus_gemma, by.x="Genus", by.y="genus")
genus_all <- merge(genus_all, Sv_gemma[,c(9,11)], by="Genus")
genus_all<- genus_all %>% 
  distinct()
genus_all <- genus_all[-1,]
genus_all$h2.id.dna<-as.numeric(genus_all$h2.id.dna)
genus_all$CospeciationScoreGenus<- as.numeric(genus_all$CospeciationScoreGenus)
cospec_genus_dna<-ggplot(genus_all, aes(x=h2.id.dna, y=CospeciationScoreGenus))+geom_point() +   geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ),
                                                  col="pink"),label.x  = 'middle', label.y = 1, size = 3)+ 
  geom_text_repel(aes(label = Genus), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", 
               label.y = 0.85, parse = TRUE, size = 3)+theme_test() + ggtitle("DNA")
cospec_genus_dna
genus_all$h2.id.rna<-as.numeric(genus_all$h2.id.rna)
cospec_genus_rna<-ggplot(genus_all, aes(x=h2.id.rna, y=CospeciationScoreGenus))+geom_point() +  
  geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), 
                                                  col="pink"),label.x  = 'middle', label.y = 0.9, size = 3)+ 
  geom_text_repel(aes(label = Genus), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, 
               label.x = "middle", label.y = 0.85, parse = TRUE, size = 3)+theme_test() + ggtitle("RNA")
cospec_genus_rna

rcorr(genus_all$h2.id.dna, genus_all$`DNA pve`, type="spearman")
rcorr(genus_all$h2.id.rna, genus_all$`RNA pve`, type="spearman")

rcorr(genus_all$h2.id.rna, genus_all$h2.id.dna, type="spearman")

library(viridis)
library(hrbrthemes)
library(ggpubr)
library(ggsci)
library(patchwork)
genus_lme4qtl_dna_m <- genus_lme4qtl_dna %>% 
  mutate(X=forcats::fct_reorder(X,as.numeric(h2_id))) %>% 
  dplyr::select(X, h2_id, h2_mating_pair, h2_residual) %>% 
  pivot_longer(!X,names_to="h2_kind", values_to="h2_estimate") %>% 
  
  ggplot(aes(x=X, y=h2_estimate))+geom_bar(aes(fill=forcats::fct_rev(h2_kind)),position="stack",stat="identity")+coord_flip()+
  scale_fill_d3() +labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")
genus_lme4qtl_dna_m

genus_lme4qtl_rna_m <- genus_lme4qtl_rna %>% 
  mutate(X=forcats::fct_reorder(X,as.numeric(h2_id))) %>% 
  dplyr::select(X, h2_id, h2_mating_pair, h2_residual) %>% 
  pivot_longer(!X,names_to="h2_kind", values_to="h2_estimate") %>% 
  
  ggplot(aes(x=X, y=h2_estimate))+geom_bar(aes(fill=forcats::fct_rev(h2_kind)),position="stack",stat="identity")+coord_flip()+
  scale_fill_d3() +labs(x="", y="Heritability estimate", fill="")+ggtitle("RNA")
genus_lme4qtl_rna_m

(genus_lme4qtl_dna_m + genus_lme4qtl_rna_m +plot_layout(guides="collect"))/(cospec_genus_dna+cospec_genus_rna)+
  plot_annotation(tag_levels = "A")
ggsave("heritability_lme4qtl_genus.pdf", width=14, height = 10)
# cow_cospec<-cowplot::plot_grid(cospec_genus_dna,cospec_genus_rna)
# (genus_lme4qtl_dna_m + genus_lme4qtl_rna_m +plot_layout(guides="collect"))/cow_cospec
