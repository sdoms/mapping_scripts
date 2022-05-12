setwd("~/Documents/research/Experiments/Final_QTL_mapping/")
library(tidyverse)
library(viridis)
theme_set(theme_test())
load("pve_snps_rna.Rdata")
out_rna<-out
load("pve_snps_dna.Rdata")
out_dna<-out
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/otu/SV17_gwscan_allchr.Rdata")

load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/all_RNA.Rdata")

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
all<- input %>% inner_join(out_rna)

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
library(ggpubr)
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

#### RNA ####
load("Results/Bacterial traits/Revisions/pve_all_snps_rna.Rdata")
out_all_rna<-out
herit <- readxl::read_excel("Results/snp_heritability/lme4qtl/lme4qtl_family.xlsx", sheet="RNA")

herit <- readxl::read_excel("Results/snp_heritability/lme4qtl/h2estimates_lme4qtl.xlsx", sheet="allRNA")
herit.pve.rna.all <- herit %>% 
  inner_join(out_all_rna, by=c("taxa"="trait"))
herit.pve.rna.all %>% filter(rlrt<0.05) %>% 
  ggplot(aes(x=h2_id, y=r2_glmm)) + geom_point() +geom_smooth(method="lm" ) +geom_text_repel(aes(label=taxa))+
  stat_poly_eq( formula=y ~ x, aes(label = paste(..eq.label.., ..rr.label.., p.value.label, sep = "~~~")), 
                parse = TRUE)+theme_pubr()+geom_abline(slope=1, linetype="dashed")
ggsave("Results/Bacterial traits/Revisions/pve_all_vs_herit_rna.pdf")
ggsave("Results/Bacterial traits/Revisions/pve_all_vs_herit_rna_family.pdf")

herit.pve.rna.all %>% filter(rlrt<0.05) %>% mutate(h2_snp=h2_id-r2_glmm) %>% select(taxa,  r2_glmm,h2_snp)%>%
  pivot_longer(cols = !taxa, names_to="variable", values_to="value") %>% 
  ggplot(aes(x=taxa, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity") +coord_flip()

herit.pve.rna.all %>% 
  select(taxa,  r2_glmm,h2_id)%>%
  mutate(taxa = forcats::fct_reorder(taxa, h2_id)) %>% 
  pivot_longer(cols = !taxa, names_to="PVE", values_to="Estimate") %>% 
  ggplot(aes(x=taxa, y=Estimate, fill=PVE))+
  geom_bar(position="dodge", stat="identity") +coord_flip()+scale_fill_d3()
ggsave("Results/Bacterial traits/Revisions/PVE_all_snps_RNA_heritability_all.pdf")
ggsave("Results/Bacterial traits/Revisions/PVE_all_snps_RNA_heritability_all_family.pdf")

herit.pve.rna.all %>% 
  filter(rlrt<0.05) %>% 
  select(taxa,  r2_glmm,h2_id)%>%
  mutate(taxa = forcats::fct_reorder(taxa, h2_id)) %>% 
  pivot_longer(cols = !taxa, names_to="PVE", values_to="Estimate") %>% 
  ggplot(aes(x=taxa, y=Estimate, fill=PVE))+
  geom_bar(position="dodge", stat="identity") +coord_flip()+scale_fill_d3()
ggsave("Results/Bacterial traits/Revisions/PVE_all_snps_RNA_heritability_sig_only.pdf")
ggsave("Results/Bacterial traits/Revisions/PVE_all_snps_RNA_heritability_sig_only_family.pdf")

### DNA ####
load("Results/Bacterial traits/Revisions/pve_all_snps_dna.Rdata")

out_all_dna<-out
herit <- readxl::read_excel("Results/snp_heritability/lme4qtl/lme4qtl_family.xlsx", sheet="DNA")

herit <- readxl::read_excel("Results/snp_heritability/lme4qtl/h2estimates_lme4qtl.xlsx", sheet="allDNA")
herit.pve.dna.all <- herit %>% 
  inner_join(out_all_dna, by=c("taxa"="trait"))
herit.pve.dna.all %>% filter(rlrt<0.05) %>% 
  ggplot(aes(x=h2_id, y=r2_glmm)) + geom_point() +geom_smooth(method="lm" ) +geom_text_repel(aes(label=taxa))+
  stat_poly_eq( formula=y ~ x, aes(label = paste(..eq.label.., ..rr.label.., p.value.label, sep = "~~~")), 
                parse = TRUE)+theme_pubr()+geom_abline(slope=1, linetype="dashed")
ggsave("Results/Bacterial traits/Revisions/pve_all_vs_herit_dna.pdf")
ggsave("Results/Bacterial traits/Revisions/pve_all_vs_herit_dna_family.pdf")

herit.pve.dna.all %>% 
  select(taxa,  r2_glmm,h2_id)%>%
  mutate(taxa = forcats::fct_reorder(taxa, h2_id)) %>% 
  pivot_longer(cols = !taxa, names_to="PVE", values_to="Estimate") %>% 
  ggplot(aes(x=taxa, y=Estimate, fill=PVE))+
  geom_bar(position="dodge", stat="identity") +coord_flip()+scale_fill_d3()
ggsave("Results/Bacterial traits/Revisions/PVE_all_snps_dna_heritability_all.pdf")
ggsave("Results/Bacterial traits/Revisions/PVE_all_snps_dna_heritability_all_family.pdf")

herit.pve.dna.all %>% 
  filter(rlrt<0.05) %>% 
  select(taxa,  r2_glmm,h2_id)%>%
  mutate(taxa = forcats::fct_reorder(taxa, h2_id)) %>% 
  pivot_longer(cols = !taxa, names_to="PVE", values_to="Estimate") %>% 
  ggplot(aes(x=taxa, y=Estimate, fill=PVE))+
  geom_bar(position="dodge", stat="identity") +coord_flip()+scale_fill_d3()
ggsave("Results/Bacterial traits/Revisions/PVE_all_snps_dna_heritability_sig_only.pdf")
ggsave("Results/Bacterial traits/Revisions/PVE_all_snps_dna_heritability_sig_only_family.pdf")


library(hrbrthemes)
data <- herit.pve.dna.all %>% 
  rowwise() %>% 
  mutate( mymean = mean(c(r2_glmm,h2_id) )) %>% 
  arrange(mymean) %>% 
  mutate(taxa=factor(taxa, taxa))

# Lollipop Plot#####
ggplot(data) +
  geom_segment( aes(x=taxa, xend=taxa, y=r2_glmm, yend=h2_id), color="grey") +
  geom_point( aes(x=taxa, y=r2_glmm), color=rgb(0.2,0.7,0.1,0.5), size=3 ) + #green
  geom_point( aes(x=taxa, y=h2_id), color=rgb(0.7,0.2,0.1,0.5), size=3 ) + #red
  coord_flip()+theme_minimal()+
  xlab("") +
  ylab("Estimate")
ggsave("../../PVE_all_snps_dna_heritability_lollipop.pdf", height=10)
ggsave("../../PVE_all_snps_dna_heritability_lollipop.png", height=10)

data %>% filter(rlrt<0.05) %>% 
ggplot() +
  geom_segment( aes(x=taxa, xend=taxa, y=r2_glmm, yend=h2_id), color="grey") +
  geom_point( aes(x=taxa, y=r2_glmm), color=rgb(0.2,0.7,0.1,0.5), size=3 ) + #green
  geom_point( aes(x=taxa, y=h2_id), color=rgb(0.7,0.2,0.1,0.5), size=3 ) + #red
  coord_flip()+theme_minimal()+
  xlab("") +
  ylab("Estimate")
ggsave("../../PVE_all_snps_dna_heritability_lollipop_sig_only.pdf", height=6)
ggsave("../../PVE_all_snps_dna_heritability_lollipop_sig_only.png", height=6)
