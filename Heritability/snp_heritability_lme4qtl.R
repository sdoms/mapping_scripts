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

xx <- readxl::read_excel("h2estimates_lme4qtl.xlsx", sheet = "DNA")
tax_table<- read.delim("../../../Phenotyping_27.02.2020/tables_core/otu_tax_table_f2_core2.csv", sep=",")
xx$name <- NA
xx$name[33:47] <- xx$taxa[33:47]
xx2<- xx %>% 
  left_join(tax_table,by=c("taxa"= "SV")) 

xx2$name[1:32] <- paste0(xx2$taxa[1:32], " (", xx2$Genus[1:32], ")")
dna_xx <- xx2 %>% 
  mutate(taxa=forcats::fct_reorder(name,as.numeric(as.character(h2_id)))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=forcats::fct_rev(h2_kind)),position="stack",stat="identity")+coord_flip()+
  scale_fill_d3() +labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")

her_dna <- xx2
dna_xx
library(ggsci)
xx <- readxl::read_excel("h2estimates_lme4qtl.xlsx", sheet = "RNA")
# tax_table<- read.delim("../../../Phenotyping_27.02.2020/tables_core/otu_tax_table_f2_core2.csv", sep=",")
xx$name <- NA
xx$name[41:70] <- xx$taxa[41:70]
xx2<- xx %>% 
  left_join(tax_table,by=c("taxa"= "SV"))

xx2$name[1:40] <- paste0(xx2$taxa[1:40], " (", xx2$Genus[1:40], ")")
xx3<- xx2 %>% 
  mutate(taxa=forcats::fct_reorder(name,as.numeric(as.character(h2_id)))) %>% 

  dplyr::select(taxa, h2_id, h2_mating_pair, h2_residual, tax_level) %>% 
  pivot_longer(!c(taxa, tax_level),names_to="h2_kind", values_to="h2_estimate")
her_rna <- xx2
    
rna_xx <- xx3%>%  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=forcats::fct_rev(h2_kind)),
                                              position="stack",stat="identity")+coord_flip()+
  scale_fill_d3() +labs(x="", y="Heritability estimate", fill="")+ggtitle("RNA")
rna_xx
library(patchwork)
dna_xx+rna_xx+plot_layout(guides="collect") + plot_annotation(tag_levels = "A")
ggsave("h2_lme4ql.svg", width = 14, height=10)

# compare to the number of associations ####
## DNA ####
results <- readxl::read_excel("../../Bacterial traits/DNA/summary_tables/genome_wide_significant/gwscan/SV_results_DNA.xlsx", skip=3)
count_taxa_SV <- results %>% 
  group_by(P) %>% 
  count(SV) %>% 
  pivot_wider(id_cols=SV,names_from=P,values_from = n)
count_taxa_SV$tax<- count_taxa_SV$SV
count_taxa_SV$SV<- NULL

results <- readxl::read_excel("../../Bacterial traits/DNA/summary_tables/genome_wide_significant/gwscan/genus_results_DNA.xlsx", sheet="all")
results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P.*", "", results$tax)

count_taxa_gen <- results %>% 
  drop_na(name) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

results <- readxl::read_excel("../../Bacterial traits/DNA/summary_tables/genome_wide_significant/gwscan/family_results_DNA.xlsx", sheet="all")

results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P_.*", "", results$tax)

count_taxa_fam <- results %>% 
  drop_na(P) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

results <- readxl::read_excel("../../Bacterial traits/DNA/summary_tables/genome_wide_significant/gwscan/order_results_DNA.xlsx", sheet="all")

results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P_.*", "", results$tax)

count_taxa_ord <- results %>% 
  drop_na(P) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

results <- readxl::read_excel("../../Bacterial traits/DNA/summary_tables/genome_wide_significant/gwscan/class_results_DNA.xlsx", sheet="all")

results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P_.*", "", results$tax)

count_taxa_class <- results %>% 
  drop_na(P) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

results <- readxl::read_excel("../../Bacterial traits/DNA/summary_tables/genome_wide_significant/gwscan/phylum_results_DNA.xlsx", sheet="all")

results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P_.*", "", results$tax)

count_taxa_phy <- results %>% 
  drop_na(P) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

count_taxa_SV$tax_level <- "ASV"
count_taxa_gen$tax_level <- "genus"
count_taxa_fam$tax_level <- "family"
count_taxa_ord$tax_level <- "order"
count_taxa_class$tax_level <- "class"
count_taxa_phy$tax_level <- "phylum"

count_DNA_taxa <- rbind(count_taxa_class, count_taxa_phy)
count_DNA_taxa <- rbind(count_DNA_taxa, count_taxa_ord)
count_DNA_taxa <- rbind(count_DNA_taxa, count_taxa_fam)
count_DNA_taxa <- rbind(count_DNA_taxa, count_taxa_gen)
count_DNA_taxa <- rbind(count_DNA_taxa, count_taxa_SV)

count_DNA_taxa$`NA` <- NULL

library(ggpmisc)
library(ggrepel)
tt_dna<-count_DNA_taxa %>% 
  left_join(her_dna, by = c("tax" = "taxa")) %>% 
  drop_na(h2_id) %>% 
  dplyr::select(1:11) %>% 
  mutate(add.P=replace_na(add.P,0), dom.P=replace_na(dom.P,0),P=replace_na(P,0) ) %>% 
  mutate(allP=add.P+P+dom.P)

ggplot(tt_dna,aes(x=allP, y=h2_id))+geom_point()+geom_smooth(method = lm)+stat_fit_glance(method = 'lm',geom = 'text',
                                                                                    aes(label = paste("P = ", signif(..p.value.., digits = 4), sep = "" ), 
                                                                                        col="pink"),label.x  = "left", label.y = 0.25, size = 3)+ geom_text_repel(aes(label = tax), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "left", label.y = 0.9, parse = TRUE, size = 3)




cor.test(tt_dna$allP,tt_dna$h2_id, method="kendall") # not significant
cor.test(tt_dna$add.P,tt_dna$h2_id, method="kendall") 
cor.test(tt_dna$dom.P,tt_dna$h2_id, method="kendall") 
cor.test(tt_dna$P,tt_dna$h2_id, method="kendall") 

count_DNA_taxa<- count_DNA_taxa %>% 
  mutate(add.P=replace_na(add.P,0), dom.P=replace_na(dom.P,0),P=replace_na(P,0) ) %>% 
  mutate(allP=add.P+P+dom.P)
write.table(count_DNA_taxa, file="../../Bacterial traits/DNA/summary_tables/genome_wide_significant/Taxa_with_number_associations.tsv")
# all_dna<- readxl::read_excel("h2estimates_lme4qtl.xlsx", sheet = "allDNA")
# tt_dna_all<-count_DNA_taxa %>% 
#   left_join(all_dna, by = c("tax" = "taxa")) %>% 
#   drop_na(h2_id) %>% 
#   dplyr::select(1:11) %>% 
#   mutate(allP=add.P+P+dom.P)
# 
# ggplot(tt_dna_all,aes(x=allP, y=h2_id))+geom_point()+geom_smooth(method = lm)+stat_fit_glance(method = 'lm',geom = 'text',
#                                                                                           aes(label = paste("P = ", signif(..p.value.., digits = 4), sep = "" ), 
#                                                                                               col="pink"),label.x  = "left", label.y = 0.25, size = 3)+ geom_text_repel(aes(label = tax), size = 3)+
#   stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "left", label.y = 0.9, parse = TRUE, size = 3)


# RNA ####
results <- readxl::read_excel("../../Bacterial traits/rna/summary_tables/genome_wide_significant/gwscan/SV_results_RNA.xlsx", skip=1)
count_taxa_SV <- results %>% 
  group_by(P) %>% 
  count(SV) %>% 
  pivot_wider(id_cols=SV,names_from=P,values_from = n)
count_taxa_SV$tax<- count_taxa_SV$SV
count_taxa_SV$SV<- NULL

results <- readxl::read_excel("../../Bacterial traits/rna/summary_tables/genome_wide_significant/gwscan/genus_results_RNA.xlsx", sheet="all")
results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P.*", "", results$tax)

count_taxa_gen <- results %>% 
  drop_na(name) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

results <- readxl::read_excel("../../Bacterial traits/rna/summary_tables/genome_wide_significant/gwscan/family_results_RNA.xlsx", sheet="all")

results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P_.*", "", results$tax)

count_taxa_fam <- results %>% 
  drop_na(P) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

results <- readxl::read_excel("../../Bacterial traits/rna/summary_tables/genome_wide_significant/gwscan/order_results_RNA.xlsx", sheet="all")

results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P_.*", "", results$tax)

count_taxa_ord <- results %>% 
  drop_na(P) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

results <- readxl::read_excel("../../Bacterial traits/rna/summary_tables/genome_wide_significant/gwscan/class_results_RNA.xlsx", sheet="all")

results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P_.*", "", results$tax)

count_taxa_class <- results %>% 
  drop_na(P) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

results <- readxl::read_excel("../../Bacterial traits/rna/summary_tables/genome_wide_significant/gwscan/phylum_results_RNA.xlsx", sheet="all")

results$tax <- gsub("_....P_.*", "", results$name)
results$tax <- gsub("_P_.*", "", results$tax)

count_taxa_phy <- results %>% 
  drop_na(P) %>% 
  group_by(P) %>% 
  count(tax) %>% 
  pivot_wider(id_cols=tax,names_from=P,values_from = n)

count_taxa_SV$tax_level <- "ASV"
count_taxa_gen$tax_level <- "genus"
count_taxa_fam$tax_level <- "family"
count_taxa_ord$tax_level <- "order"
count_taxa_class$tax_level <- "class"
count_taxa_phy$tax_level <- "phylum"

count_rna_taxa <- rbind(count_taxa_class, count_taxa_phy)
count_rna_taxa <- rbind(count_rna_taxa, count_taxa_ord)
count_rna_taxa <- rbind(count_rna_taxa, count_taxa_fam)
count_rna_taxa <- rbind(count_rna_taxa, count_taxa_gen)
count_rna_taxa <- rbind(count_rna_taxa, count_taxa_SV)

count_rna_taxa$`NA` <- NULL

library(ggpmisc)
library(ggrepel)
tt_rna<-count_rna_taxa %>% 
  left_join(her_rna, by = c("tax" = "taxa")) %>% 
  drop_na(h2_id) %>% 
  dplyr::select(1:11) %>% 
  mutate(add.P=replace_na(add.P,0), dom.P=replace_na(dom.P,0),P=replace_na(P,0) ) %>% 
  mutate(allP=add.P+P+dom.P)

ggplot(tt_rna,aes(x=allP, y=h2_id))+geom_point()+geom_smooth(method = lm)+stat_fit_glance(method = 'lm',geom = 'text',
                                                                                          aes(label = paste("P = ", signif(..p.value.., digits = 4), sep = "" ), 
                                                                                              col="pink"),label.x  = "left", label.y = 0.25, size = 3)+ geom_text_repel(aes(label = tax), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "left", label.y = 0.9, parse = TRUE, size = 3)




cor.test(tt_rna$allP,tt_rna$h2_id, method="kendall") # not significant
cor.test(tt_rna$add.P,tt_rna$h2_id, method="kendall") 
cor.test(tt_rna$dom.P,tt_rna$h2_id, method="kendall") 
cor.test(tt_rna$P,tt_rna$h2_id, method="kendall") 

count_rna_taxa<- count_rna_taxa %>% 
  mutate(add.P=replace_na(add.P,0), dom.P=replace_na(dom.P,0),P=replace_na(P,0) ) %>% 
  mutate(allP=add.P+P+dom.P)
write.table(count_rna_taxa, file="../../Bacterial traits/RNA/summary_tables/genome_wide_significant/Taxa_with_number_associations.tsv")

write_excel_csv2(her_dna, file="sign_heritability_dna.csv")
write_excel_csv2(her_rna, file="sign_heritability_rna.csv")

her_all <- her_dna %>% 
  full_join(her_rna, by="name")
brks <- seq(-1, 1, 0.2)
lbls = as.character(c(seq(1, 0, -0.2), seq(0.2, 1, 0.2)))
d3_dub<- c( "#2CA02CFF","#2CA02CFF","#1F77B4FF","#1F77B4FF", "#FF7F0EFF","#FF7F0EFF")
her_barplot <- her_all%>% 
  mutate(her_dna=h2_id.x, h2_mating_pair_dna=h2_mating_pair.x,h2_residual_dna=h2_residual.x, 
         her_rna=h2_id.y*(-1), h2_mating_pair_rna=h2_mating_pair.y*(-1), h2_residual_rna=h2_residual.y*-1) %>% 
  mutate(name=forcats::fct_reorder(name,as.numeric(her_rna))) %>% 
  dplyr::select(name,her_dna, her_rna, h2_mating_pair_dna, h2_residual_dna, h2_mating_pair_rna, h2_residual_rna) %>% 
  pivot_longer(!name,names_to="h2_kind", values_to="h2_estimate") %>% 
  
  ggplot(aes(x=name, y=h2_estimate))+geom_bar(aes(fill=forcats::fct_relevel(h2_kind, 
                                                                            c( "h2_residual_rna",
                                                                               "h2_residual_dna",
                                                                               "h2_mating_pair_dna",
                                                                               "h2_mating_pair_rna",
                                                                               "her_rna", "her_dna"))),
                                              position="stack",stat="identity", width=0.9)+
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip()+
  scale_fill_manual(values=d3_dub) +labs(x="", y="Heritability estimate", fill="") +
  theme_bw() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks.y  = element_blank())+
  geom_hline(yintercept = 0, color="black")
her_barplot
ggsave("heritabilty_DNA_RNA_lme4qtl.pdf")
ggsave("heritabilty_DNA_RNA_lme4qtl.svg")
gen_melt <-melt(genus, id.vars = "genus", measure.vars = c("DNA", "RNA"))
# X Axis Breaks and Labels 

tot_gen <-ggplot(gen_melt, aes(x = reorder(genus, abs(value)), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(y="SNP heritability", x="", fill="") +
  theme_bw() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks.y  = element_blank()) +   # Centre plot title
  scale_fill_brewer(palette = "Dark2")  # Color palette
tot_gen
ggsave("genus_SNP_heritability.pdf", plot=tot_gen)
