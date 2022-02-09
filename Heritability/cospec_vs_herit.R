setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/snp_heritability/cospeciation/")

library(readxl)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpmisc)
library(patchwork)
#color_palettes
cbbPalette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", '#6a3d9a', "#ffff99", "#b15928", "#40e0d0", "#ffff00", "#ee82ee", "#fb8072", "#191970", "#fdb462", "#CBD4D4", "#fccde5", "#d9d9d9", "#bc80bd", "#32cd32", "#ffed6f", "lightcoral" )
colors_palette = c("#a6cee3", "#1f78b4", "#f0f8ff", "#b2df8a", "#33a02c", "#fb9a99", "#ff69b4", "#e31a1c", "#ffe4e1", "#cd853f", "#f0bf6f", "#ff7f00", "#cab2d6", "#ba3d9a","#6a3d94","#ffff99", "#B15928", "#ffb90f", "#8b4513")

lme4qtl_dna<- readxl::read_excel("lme4qtl.xlsx", sheet="DNA")

lme4qtl_rna<- readxl::read_excel("lme4qtl.xlsx", sheet="RNA")

cospec <- read.table("../cospeciation_rates.csv", sep=";", header=T)

lme4qtl_dna %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+stat_cor(label.x = 0.4, label.y=0.8, col="blue")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("DNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_all_DNA.pdf")
lme4qtl_rna %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+stat_cor(label.x = 0.4, label.y=0.8, col="blue")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("RNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_all_RNA.pdf")

lme4qtl_dna %>% filter(rlrt<0.05) %>%  inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+stat_cor(label.x = 0.4, label.y=0.8, col="blue")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("DNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_sign_only_DNA.pdf")
lme4qtl_rna %>% filter(rlrt<0.05) %>% inner_join(cospec, by=c("taxa"="tax")) %>% 
  ggplot(aes(x=h2_id, y=cospeciation))+geom_point() +geom_smooth(method="lm")+stat_cor(label.x = 0.4, label.y=0.8, col="blue")+
  labs(x="h2SNP", y="Cospeciation rate")+ggtitle("RNA")+geom_text_repel(aes(label = taxa), size = 3)
ggsave("cospeciation_lme4qtl_sign_only_RNA.pdf")



