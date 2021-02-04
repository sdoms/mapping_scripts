setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")
library(tidyverse)
library(patchwork)
add_res <- read.table("ALL_sig_snps_add_DNA_P.txt")
dom_res <- read.table("ALL_sig_snps_dom_P_DNA.txt")
P_res <- read.table("ALL_sig_snps_P_DNA.txt")
allP <- read.table("ALL_sig_snps_allP_DNA.txt")

allP$degree_of_dominance <- allP$dom.T/allP$add.T
add_res$degree_of_dominance <- add_res$dom.T/add_res$add.T
dom_res$degree_of_dominance <- dom_res$dom.T/dom_res$add.T
P_res$degree_of_dominance <- P_res$dom.T/P_res$add.T

ggplot(allP, aes(x=add.T))+ geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = mean(add.T)), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07")
ggplot(allP, aes(x=degree_of_dominance))+ geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = 1), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07")
             
             
add <- ggplot(allP, aes(x=add.T)) + geom_histogram(colour="black", fill="white",breaks=seq(round(min(allP$add.T, na.rm=T),3)-0.2, 
                                                              round(max(allP$add.T, na.rm=T),3)+0.2, by=0.2)) + labs(x="Additive effect Z-score ")+ theme(axis.title = element_text(size=9))
dom <- ggplot(allP, aes(x=dom.T)) +geom_histogram(colour="black", fill="white",
                                                  breaks=seq(round(min(allP$dom.T, na.rm=T),3)-0.2, 
                                                             round(max(allP$dom.T, na.rm=T),3)+0.2, by=0.2))+ 
  labs(x="Dominance effect Z-score")+ 
  theme(axis.title = element_text(size=9))  
dom_degree <- ggplot(allP, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(round(min(allP$degree_of_dominance, na.rm=T),3)-0.03, round(max(allP$degree_of_dominance, na.rm=T),3)+0.03, by=0.03))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9))  

(add+dom) /dom_degree+ plot_annotation(title="All DNA", tag_levels = "A")
ggsave("DNA_zscores_degree_sep.pdf")
melted_allP <- reshape2::melt(allP, id.vars="marker",measure.vars=c("add.T", "dom.T"))
ad<-ggplot(melted_allP, aes(x=value))+ geom_histogram(aes(color = variable, fill = variable),
                                                    alpha = 0.2, position = "identity", bins=50) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + labs(x="z-score") +theme_test()
dom_dens <- ggplot(allP, aes(x=degree_of_dominance))+ geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = 1), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07") +theme_test()
ad/dom_dens
ggsave("DNA_zscores_degree.pdf")
allP_D<- allP
#### RNA
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")
# 
# add_res <- read.table("ALL_sig_snps_add_DNA_P.txt")
# dom_res <- read.table("ALL_sig_snps_dom_P_DNA.txt")
# P_res <- read.table("ALL_sig_snps_P_DNA.txt")
allP <- read.table("RNAALL_sig_snps_allP.txt")

allP$degree_of_dominance <- allP$dom.T/allP$add.T
add_res$degree_of_dominance <- add_res$dom.T/add_res$add.T
dom_res$degree_of_dominance <- dom_res$dom.T/dom_res$add.T
P_res$degree_of_dominance <- P_res$dom.T/P_res$add.T

ggplot(allP, aes(x=add.T))+ geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = mean(add.T)), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07")
ggplot(allP, aes(x=degree_of_dominance))+ geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = 1), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07")


add <- ggplot(allP, aes(x=add.T)) + geom_histogram(colour="black", fill="white",
                                                   breaks=seq(round(min(allP$add.T, na.rm=T),3)-0.2, 
                                                              round(max(allP$add.T, na.rm=T),3)+0.2, by=0.3)) + labs(x="Additive effect Z-score ")+ theme(axis.title = element_text(size=9))
dom <- ggplot(allP, aes(x=dom.T)) +geom_histogram(colour="black", fill="white",
                                                  breaks=seq(round(min(allP$dom.T, na.rm=T),3)-0.2, 
                                                             round(max(allP$dom.T, na.rm=T),3)+0.2, by=0.3))+ 
  labs(x="Dominance effect Z-score")+ 
  theme(axis.title = element_text(size=9))  
dom_degree <- ggplot(allP, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(round(min(allP$degree_of_dominance, na.rm=T),3)-0.05, round(max(allP$degree_of_dominance, na.rm=T),3)+0.05, by=0.05))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9))  

(add+dom) /dom_degree+ plot_annotation(title="All RNA", tag_levels = "A")
ggsave("RNA_zscores_degree_sep.pdf")
melted_allP <- reshape2::melt(allP, id.vars="marker",measure.vars=c("add.T", "dom.T"))
ad<-ggplot(melted_allP, aes(x=value))+ geom_histogram(aes(color = variable, fill = variable),
                                                      alpha = 0.2, position = "identity", bins=50) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + labs(x="z-score") +theme_test()
dom_dens <- ggplot(allP, aes(x=degree_of_dominance))+ geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = 1), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07") +theme_test()
ad/dom_dens
ggsave("RNA_zscores_degree.pdf")

allP_R <- allP
allP_R$DNAorRNA <- "RNA"
allP_D$DNAorRNA <- "DNA"
allP_all <- rbind(allP_D, allP_R)

melted_allP <- reshape2::melt(allP_all, id.vars="marker",measure.vars=c("add.T", "dom.T"))
ad<-ggplot(melted_allP, aes(x=value))+ geom_histogram(aes(color = variable, fill = variable),
                                                      alpha = 0.2, position = "identity", bins=80) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + labs(x="z-score") +theme_test()
dom_dens <- ggplot(allP_all, aes(x=degree_of_dominance))+ geom_histogram(colour="black", fill="white",
                                                                         breaks=seq(round(min(allP_all$degree_of_dominance, na.rm=T),3)-0.05, 
                                                                                    round(max(allP_all$degree_of_dominance, na.rm=T),3)+0.05, by=0.05))+ 
  labs(x="Degree of dominance (d/a)") + geom_vline(xintercept = 1, colour="red", linetype=2)+theme_test()
  
ad/dom_dens
ggsave("RNA_DNA_zscores_degree.pdf")
