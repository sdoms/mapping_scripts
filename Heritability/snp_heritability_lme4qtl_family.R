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
library(argyle)

gemma <-"/usr/local/bin/gemma"
DNAorRNA="DNA"
PCAs<- read.table("~/Documents/Research/Experiments/Final_QTL_mapping/Scripts/Revisions/output/PCA_projection.txt")
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Plots/Effect_plots/ref_alt_min_major.Rdata")
alpha <- read.csv("../../../Phenotyping_27.02.2020/alpha_values_otu.csv") %>% filter(Family=="F2") %>% pivot_wider(id_cols = "ID", names_from="variable", values_from="value")
beta <- read.csv("../../../Phenotyping_27.02.2020/bray_curtis_axis.csv")[,1:6]
alpha_beta <- left_join(alpha, beta, by=c("ID"="X")) 
for (DNAorRNA in c("DNA", "RNA")){
  ##### Phenotypes ####
  pheno<-read.csv("../../../Phenotypes_histology_new.csv", sep=";", header=T)
  rownames(pheno)<-pheno$Mouse_Name
  pheno$id <- pheno$Mouse_Name
  pheno_pca<- PCAs %>% rownames_to_column("id") %>% left_join(pheno, by="id")
  rownames(pheno_pca)<-pheno_pca$id
  #### Genotypes ####
  geno<- BEDMatrix("../../../Cleaning_snps/clean_f2")
  rownames(geno)<- substr(rownames(geno),5 ,19)
  rownames(geno) <- gsub(x = rownames(geno), pattern = "\\/", replacement = ".")  
  colnames(geno)<- substr(colnames(geno),1 ,nchar(colnames(geno))-2)
  indi_all <- rownames(geno)
  load("../../../Cleaning_snps/clean_snps.Rdata")
  pheno_pca <- pheno_pca[indi_all, ]
  individuals <- pheno_pca$id
  geno <- geno[individuals,]
  indi_all <- rownames(geno)
  
  # ## Alpha and beta ####
  # div <- alpha_beta[which(grepl(paste0("_",DNAorRNA),alpha_beta$ID)),]
  # div$ID<- substr(div$ID, 1,nchar(div$ID)-4)
  # 
  # da <- div %>% left_join( pheno_pca, by=c("ID"="id"))
  # individuals <-da$ID
  # 
  # ##### kinship matrix ######
  # # genotype kinship
  # F2 <- read.plink("../../../Cleaning_snps/clean_f2")
  # F2 <- recode(F2, "relative")
  # 
  # colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")  
  # 
  # #Select the individuals from the genotype file 
  # F2 <- as.matrix(as.data.frame(F2))
  # 
  # ind <- ncol(F2)
  # F2_part <- apply(F2[,7:ind], 2, as.numeric)# make sure that the genotype information(0,1,2) is numeric
  # 
  # F2_part <- F2_part[,individuals]
  # 
  # F2 <- cbind(genotypes_ref[,c(7,6,5)],F2_part)
  # 
  # write.table(F2,"geno.txt",sep = " ",quote = FALSE,row.names = FALSE,
  #             col.names = FALSE)
  # write.table(data.frame(rep(1, nrow(da))),"pheno.txt", quote=F, row.names=F, col.names=F)
  # 
  # # # compute kinship matrix
  # 
  # system(paste0(gemma, " -g geno.txt -gk 1 -p pheno.txt -debug -o kinship_", DNAorRNA  ))
  # div <- column_to_rownames(div, var = "ID")
  # results<-data.frame()
  # for (tx in 1:ncol(div)){
  #   taxa2 <- div[tx]
  #   tax <- names(taxa2)
  #   cat("Running model for ", tax, ".\n")
  #   names(taxa2)<-"tax"
  #   lmm.data <- merge(taxa2, pheno_pca, by="row.names")
  #   rownames(lmm.data)<- lmm.data$Row.names
  #   lmm.data$Row.names<- NULL
  #   individuals <-rownames(lmm.data)
  #   
  #   kinship <- read.table(paste0("./output/kinship_", DNAorRNA,".cXX.txt"))
  #   kinship<- as.matrix(kinship)
  #   rownames(kinship)<- colnames(F2)[4:ncol(F2)]
  #   colnames(kinship)<- colnames(F2)[4:ncol(F2)]
  #   
  #   #snp <- sub[1]
  #   df<-data.frame(lmm.data[rownames(kinship),])
  #   
  #   null_model <- relmatLmer(tax ~(1|mating.pair)+(1|family.x)+(1|id), df,relmat=list(id=kinship))
  #   # m1 <-update(null_model, . ~ . + (1|family.x))
  #   # r1<- residuals(null_model)
  #   # qqnorm(r1)
  #   # qqline(r1)
  #   # hist(r1, breaks = 30)
  #   
  #   # heritability
  #   h2 <- VarProp(null_model)
  #   
  #   null_model_reduced <- update(null_model, . ~ . - (1|mating.pair)-(1|family.x))
  #   m1_null <- update(null_model, . ~ . - (1|id))
  # 
  #   rlrt_h2 <- exactRLRT(
  #     null_model_reduced, # the reduced model with only the effect to be tested
  #     mA = null_model, # the full model under the alternative
  #     m0 = m1_null, # the model under the null
  #     seed = 1
  #   )
  #   rlrt_h2
  #   
  #   lrt_h2 <- anova(m1_null, null_model)
  #   lrt_h2
  #   scv <- cv(taxa2$tax)
  #   
  #   results[tax,"h2_id"]<- h2[1,6]
  #   results[tax,"h2_mating_pair"]<- h2[2,6]
  #   results[tax,"h2_family"]<- h2[3,6]
  #   results[tax,"h2_residual"]<- h2[3,6]
  #    results[tax,"rlrt"]<- rlrt_h2$p.value
  #   results[tax,"lrt"]<- lrt_h2$`Pr(>Chisq)`[2]
  #   results[tax,"Sample_coef_variation"]<- scv
  # }
  # write.table(results, paste0("Alpha_Beta_",DNAorRNA,"_snp_heritabilitylme4_mating_pair_family.csv"), sep=";")
  
  
  
  ##### taxa ####
  
  traits <- c("otu", "genus", "family", "order", "class", "phylum")
  for (trait in traits){
    taxa <- read.csv(paste0("../../../Phenotyping_27.02.2020/tables_core/",trait, "_table_f2_core.csv"), header=T, row.names = 1)
    taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),]
    names_taxa<- substr(rownames(taxa), 1,nchar(rownames(taxa))-4)
    taxa <- as.data.frame(sapply(taxa, function(x) inv.logit(x), USE.NAMES = T))
    rownames(taxa)<- names_taxa
    
    da <- taxa[1] %>% rownames_to_column("id") %>% left_join( pheno_pca, by="id")
    individuals <-da$id
    
    ##### kinship matrix ######
    # genotype kinship
    F2 <- read.plink("../../../Cleaning_snps/clean_f2")
    F2 <- recode(F2, "relative")
    
    colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")  
    
    #Select the individuals from the genotype file 
    F2 <- as.matrix(as.data.frame(F2))
    
    ind <- ncol(F2)
    F2_part <- apply(F2[,7:ind], 2, as.numeric)# make sure that the genotype information(0,1,2) is numeric
    
    F2_part <- F2_part[,individuals]
    
    F2 <- cbind(genotypes_ref[,c(7,6,5)],F2_part)
    
    write.table(F2,"geno.txt",sep = " ",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
    write.table(data.frame(rep(1, nrow(da))),"pheno.txt", quote=F, row.names=F, col.names=F)
    
    # # compute kinship matrix
    
    system(paste0(gemma, " -g geno.txt -gk 1 -p pheno.txt -debug -o kinship_", DNAorRNA  ))
    
    
    
    results <- data.frame()
    for (tx in 1:ncol(taxa)){
      if(DNAorRNA=="RNA"&trait=="otu"&tx==8|tx==74) next
      
      taxa2 <- taxa[tx]
      tax <- names(taxa2)
      cat("Running model for ", tax, ".\n")
      names(taxa2)<-"tax"
      lmm.data <- merge(taxa2, pheno_pca, by="row.names")
      rownames(lmm.data)<- lmm.data$Row.names
      lmm.data$Row.names<- NULL
      individuals <-rownames(lmm.data)
      
      kinship <- read.table(paste0("./output/kinship_", DNAorRNA,".cXX.txt"))
      kinship<- as.matrix(kinship)
      rownames(kinship)<- colnames(F2)[4:ncol(F2)]
      colnames(kinship)<- colnames(F2)[4:ncol(F2)]
      
      #snp <- sub[1]
      df<-data.frame(lmm.data)
      
      null_model <- relmatLmer(tax ~(1|mating.pair)+ (1|family.x)+(1|id), df,relmat=list(id=kinship))
     # m1 <-update(null_model, . ~ . + (1|family.x))
      # r1<- residuals(null_model)
      # qqnorm(r1)
      # qqline(r1)
      # hist(r1, breaks = 30)
      
      # heritability
      h2 <- VarProp(null_model)
      
      null_model_reduced <- update(null_model, . ~ . - (1|mating.pair)-(1|family.x))
      m1_null <- update(null_model, . ~ . - (1|id))
      
      rlrt_h2 <- exactRLRT(
       # null_model_reduced, # the reduced model with only the effect to be tested
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
      results[tax,"h2_family"]<- h2[3,6]
      results[tax,"h2_residual"]<- h2[4,6]
      results[tax,"rlrt"]<- rlrt_h2$p.value
      results[tax,"lrt"]<- lrt_h2$`Pr(>Chisq)`[2]
      results[tax,"Sample_coef_variation"]<- scv
    }
    write.table(results, paste0(trait,"_",DNAorRNA,"_snp_heritabilitylme4_mating_pair_family.csv"), sep=";")
  }
}

theme_set(theme_test())
lme4qtl_dna<- readxl::read_excel("lme4qtl_family.xlsx", sheet="DNA")

lme4qtl_rna<- readxl::read_excel("lme4qtl_family.xlsx", sheet="RNA")

cospec <- read.table("../cospeciation_rates.csv", sep=";", header=T)

lme4qtl_dna %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+
  stat_cor(label.x = 0.4, label.y=0.8, col="blue", method = "spearman")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("DNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_family_all_DNA.pdf")
lme4qtl_rna %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+
  stat_cor(label.x = 0.4, label.y=0.8, col="blue", method="spearman")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("RNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_family_all_RNA.pdf")


lme4qtl_dna %>% filter(rlrt<0.05) %>%  inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+stat_cor(label.x = 0.4, label.y=0.8, col="blue")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("DNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_family_sign_only_DNA.pdf")
lme4qtl_rna %>% filter(rlrt<0.05) %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+stat_cor(label.x = 0.4, label.y=0.8, col="blue")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("RNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_family_sign_only_RNA.pdf")

col_pal<- c("light grey","#8F2D56", "#FBB13C","#218380" )


lme4qtl_dna %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 

  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_DNA_family.pdf", height=15)

lme4qtl_dna %>% filter(rlrt<0.05) %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_DNA_family_sign.pdf", height=10)

lme4qtl_dna %>% filter(lrt<0.05) %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_DNA_family_sign_lrt.pdf", height=10)

lme4qtl_rna %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("RNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_RNA_family.pdf", height=15)

lme4qtl_rna %>% filter(rlrt<0.05) %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("RNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_RNA_family_sign.pdf", height=10)

lme4qtl_rna %>% filter(lrt<0.05) %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("RNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_RNA_family_sign_lrt.pdf", height=10)


# scatter dna vs rNa
inner_join(lme4qtl_dna, lme4qtl_rna, by="taxa") %>% 
ggplot(aes(x=h2_id.x, y=h2_id.y))+geom_point() +theme_test()+labs(x="h2 DNA", y="h2 RNA") +
geom_smooth(method="lm")+stat_cor()+
  geom_text_repel(aes(label = taxa), size = 3)+ 
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", alpha=0.7) 
ggsave("correlation_DNA_RNA_family.pdf")


inner_join(lme4qtl_dna, lme4qtl_rna, by="taxa") %>% filter(rlrt.x<0.05&rlrt.y<0.05) %>% 
  ggplot(aes(x=h2_id.x, y=h2_id.y))+geom_point() +theme_test()+labs(x="h2 DNA", y="h2 RNA") +
  geom_smooth(method="lm")+stat_cor()+
  geom_text_repel(aes(label = taxa), size = 3)+ 
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", alpha=0.7) 
ggsave("correlation_DNA_RNA_family_sig_only.pdf")

### standardised ####

lme4qtl_dna<- readxl::read_excel("lme4_mating_pair_family_standardised.xlsx", sheet="DNA")

lme4qtl_rna<- readxl::read_excel("lme4_mating_pair_family_standardised.xlsx", sheet="RNA")


cospec <- read.table("../cospeciation_rates.csv", sep=";", header=T)

lme4qtl_dna %>% rename("taxa"="...1") %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+
  stat_cor(label.x = 0.4, label.y=0.8, col="blue", method="spearman")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("DNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_family_all_stand_DNA.pdf")
lme4qtl_rna %>% rename("taxa"="...1") %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+
  stat_cor(label.x = 0.4, label.y=0.9, col="blue", method="spearman")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("RNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_family_all_stand_RNA.pdf")


lme4qtl_dna %>% rename("taxa"="...1")%>% filter(rlrt<0.05) %>%  inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+stat_cor(label.x = 0.4, label.y=0.8, col="blue")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("DNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_family_sign_only_DNA.pdf")
lme4qtl_rna %>% rename("taxa"="...1") %>% filter(rlrt<0.05) %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+stat_cor(label.x = 0.4, label.y=0.8, col="blue")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("RNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_family_sign_only_RNA.pdf")

col_pal<- c("light grey","#8F2D56", "#FBB13C","#218380" )


lme4qtl_dna %>% rename("taxa"="...1")%>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_DNA_family_stand.pdf", height=15)

lme4qtl_dna %>% rename("taxa"="...1")%>% filter(rlrt<0.05) %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("DNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_DNA_family_sign_stand.pdf", height=10)

lme4qtl_rna %>% rename("taxa"="...1")%>% filter(rlrt<0.05) %>% 
  mutate(taxa=forcats::fct_reorder(taxa,as.numeric(h2_id))) %>% 
  dplyr::select(taxa, h2_id, h2_mating_pair,h2_family, h2_residual) %>% 
  pivot_longer(!taxa,names_to="h2_kind", values_to="h2_estimate") %>% 
  mutate(h2_kind=factor(h2_kind, levels=c("h2_residual","h2_family","h2_mating_pair","h2_id" ))) %>% 
  
  ggplot(aes(x=taxa, y=h2_estimate))+geom_bar(aes(fill=h2_kind),position="stack",stat="identity")+coord_flip() +
  labs(x="", y="Heritability estimate", fill="")+ggtitle("RNA")+ scale_fill_manual(values = col_pal)
ggsave("all_heritability_estimates_RNA_family_sign_stand.pdf", height=10)
