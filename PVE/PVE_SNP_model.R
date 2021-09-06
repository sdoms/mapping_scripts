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
library(MuMIn)
#library(qqman)

# Parameters ####
DNAorRNA <- "RNA"
# read in combined sig.summaries file, the script will calculate for every marker-trait association
load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/all_RNA.Rdata")
# the loaded file is called all_files (just an rbind of sig.summaries results)
cat("Reading in phenotypes and covariates. \n")
pheno <- read.csv("./Phenotypes_new.csv", sep=";", header=T)
rownames(pheno)<-pheno$Mouse_Name
pheno$id <- pheno$Mouse_Name

# Genotypes ------------
cat("Reading in genotypes. \n")
geno_in<- BEDMatrix("Cleaning_snps/clean_f2")
rownames(geno_in)<- substr(rownames(geno_in),4 ,18)
rownames(geno_in) <- gsub(x = rownames(geno_in), pattern = "\\/", replacement = ".")  
colnames(geno_in)<- substr(colnames(geno_in),1 ,nchar(colnames(geno_in))-2)

load("Cleaning_snps/clean_snps.Rdata")

# recode from major allele count coding to -1(homo min), 0(hetero), 1(homo ref)
add.matrix <-ifelse(geno_in[,] == 0, 1, ifelse(geno_in[,] == 1, 0,  ifelse(geno_in[,] == 2, -1, NA)))

# Dominance
dom.matrix <- ifelse(abs(geno_in[,]) == 2, 0, ifelse(geno_in[,] == 1, 1,  ifelse(geno_in[,] == 0, 0, NA)))
rownames(dom.matrix)<- rownames(add.matrix)
colnames(dom.matrix)<- colnames(dom.matrix)

out<-data.frame(all_files[2:32],r2_glmm=NA,r2_lr=NA,r2_lm=NA)
# loops over all SNP-trait combo's
for (i in 1:nrow(all_files)){
  # parameters
  trait <- all_files$tax_level[i]
  snp<- all_files$peak.snp[i]
  tx <- all_files$trait[i]
  chr<- all_files$chr[i]
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
    kinship <- read.table(paste0("Kinship_DNA/kinship_chr",chr,".cXX.txt"))
  }else if (DNAorRNA=="RNA"){
    kinship <- read.table(paste0("./Kinship_RNA/kinship_chr",chr,".cXX.txt"))
  }
  kinship<- as.matrix(kinship)
  rownames(kinship)<- indi_all
  colnames(kinship)<-indi_all
  
  df<-data.frame(lmm.data,ad=add.mat[,snp],dom=dom.mat[,snp],gt=geno[,snp])
  df<-df[is.na(df$ad)==F & is.na(df$tax)==F &is.na(df$gt)==F & is.na(df$dom)==F,]
  null_model <- relmatLmer(tax ~ (1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  model <- relmatLmer(tax ~ ad+dom+(1|mating.pair)+ (1|id), df,relmat=list(id=kinship))
  
  out[i,"r2_glmm"]<-r.squaredGLMM(model,null_model)[1] ### this is the good one
  out[i,"r2_lr"]<-r.squaredLR(model,null = null_model)[1]
  out[i,"r2_lm"]<-summary(lm(tax~ad+dom, df))$adj.r.squared
  
  
  }


save(out, file="pve_snps_rna.Rdata")
  save(out, file="pve_snps_dna.Rdata")

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
