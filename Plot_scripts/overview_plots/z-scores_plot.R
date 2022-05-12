setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")
library(tidyverse)
library(patchwork)
library(ggsci)
theme_set(theme_test())
load("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/DNA_all_sig_addP_snps_pval_corr_SW.Rdata")
load("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/DNA_all_sig_domP_snps_pval_corr_SW.Rdata")
load("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/DNA_all_sig_P_snps_pval_corr_SW.Rdata")
load("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/DNA_all_sig_snps_pval_corr_SW.Rdata")
# 
# add_res <- read.table("ALL_sig_snps_add_DNA_P.txt")
# dom_res <- read.table("ALL_sig_snps_dom_P_DNA.txt")
# P_res <- read.table("ALL_sig_snps_P_DNA.txt")
# sig.all <- read.table("ALL_sig_snps_allP_DNA.txt")

sig.all$degree_of_dominance <- sig.all$dom.T/abs(sig.all$add.T)
sig.add.p$degree_of_dominance <- sig.add.p$dom.T/abs(sig.add.p$add.T)
sig.dom.p$degree_of_dominance <- sig.dom.p$dom.T/abs(sig.dom.p$add.T)
sig.P$degree_of_dominance <- sig.P$dom.T/abs(sig.P$add.T)
sig.all$deg_beta <- sig.all$dom.Beta/abs(sig.all$add.Beta)

# ggplot(sig.all, aes(x=add.T))+ geom_density(aes(y = ..count..), fill = "lightgray") +
#   geom_vline(aes(xintercept = mean(add.T)), 
#              linetype = "dashed", size = 0.6,
#              color = "#FC4E07")
# ggplot(sig.all, aes(x=degree_of_dominance))+ geom_density(aes(y = ..count..), fill = "lightgray") +
#   geom_vline(aes(xintercept = 1), 
#              linetype = "dashed", size = 0.6,
#              color = "#FC4E07")
             
             
add <- ggplot(sig.all, aes(x=add.T)) + geom_histogram(colour="black", fill="white",breaks=seq(round(min(sig.all$add.T, na.rm=T),3)-0.2, 
                                                              round(max(sig.all$add.T, na.rm=T),3)+0.2, by=0.2)) + labs(x="Additive effect Z-score ")+ 
  theme(axis.title = element_text(size=9))
dom <- ggplot(sig.all, aes(x=dom.T)) +geom_histogram(colour="black", fill="white",
                                                  breaks=seq(round(min(sig.all$dom.T, na.rm=T),3)-0.2, 
                                                             round(max(sig.all$dom.T, na.rm=T),3)+0.2, by=0.2))+ 
  labs(x="Dominance effect Z-score")+ 
  theme(axis.title = element_text(size=9))  
dom_degree <- ggplot(sig.all, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(round(min(sig.all$degree_of_dominance, na.rm=T),3)-0.03, round(max(sig.all$degree_of_dominance, na.rm=T),3)+0.03, by=0.03))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9))  

(add+dom) /dom_degree+ plot_annotation(title="All DNA", tag_levels = "A")
ggsave("DNA_zscores_degree_sep_pval_corr.pdf")
melted_sig.all <- reshape2::melt(sig.all, id.vars="marker",measure.vars=c("add.T", "dom.T"))
ad<-ggplot(melted_sig.all, aes(x=value))+ geom_histogram(aes(color = variable, fill = variable),
                                                    alpha = 0.2, position = "identity", bins=50) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + labs(x="z-score") +theme_test()
dom_dens <- ggplot(sig.all, aes(x=degree_of_dominance))+ geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = 1), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07") +theme_test()
ad/dom_dens
ad/dom_degree
ggsave("DNA_zscores_degree_pval_corr.pdf")
sig.all_D<- sig.all

ggplot(melted_sig.all, aes(x=fct_reorder(marker, value), y=value))+geom_bar(aes(fill=variable), stat = "identity", position=position_dodge())+
  theme(axis.text.x = element_blank())
ggplot(sig.all_D, aes(x=add.T, y=dom.T))+geom_point()

#### RNA
setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")

load("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/RNA_all_sig_addP_snps_pval_corr_SW.Rdata")
load("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/RNA_all_sig_domP_snps_pval_corr_SW.Rdata")
load("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/RNA_all_sig_P_snps_pval_corr_SW.Rdata")
load("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/RNA_all_sig_snps_pval_corr_SW.Rdata")

sig.all$degree_of_dominance <- sig.all$dom.T/abs(sig.all$add.T)
sig.add.p$degree_of_dominance <- sig.add.p$dom.T/abs(sig.add.p$add.T)
sig.dom.p$degree_of_dominance <- sig.dom.p$dom.T/abs(sig.dom.p$add.T)
sig.P$degree_of_dominance <- sig.P$dom.T/abs(sig.P$add.T)

sig.all$deg_beta <- sig.all$dom.Beta/abs(sig.all$add.Beta)


ggplot(sig.all, aes(x=abs(add.Beta), y=deg_beta))+geom_point()
ggplot(sig.all, aes(x=dom.Beta, y=add.Beta))+geom_point() +geom_smooth(method="lm")
ggsave("dom.beta_add.beta_RNA_pval_corr.pdf")
ggplot(sig.all, aes(x=dom.Beta, y=abs(add.Beta)))+geom_point() +geom_smooth(method="lm")
ggsave("dom.beta_abs_add.beta_RNA_pval_corr.pdf")

# ggplot(sig.all, aes(x=add.T))+ geom_density(aes(y = ..count..), fill = "lightgray") +
#   geom_vline(aes(xintercept = mean(add.T)), 
#              linetype = "dashed", size = 0.6,
#              color = "#FC4E07")
# ggplot(sig.all, aes(x=degree_of_dominance))+ geom_density(aes(y = ..count..), fill = "lightgray") +
#   geom_vline(aes(xintercept = 1), 
#              linetype = "dashed", size = 0.6,
#              color = "#FC4E07")


add <- ggplot(sig.all, aes(x=add.T)) + geom_histogram(colour="black", fill="white",
                                                   breaks=seq(round(min(sig.all$add.T, na.rm=T),3)-0.2, 
                                                              round(max(sig.all$add.T, na.rm=T),3)+0.2, by=0.3)) + labs(x="Additive effect Z-score ")+ theme(axis.title = element_text(size=9))
dom <- ggplot(sig.all, aes(x=dom.T)) +geom_histogram(colour="black", fill="white",
                                                  breaks=seq(round(min(sig.all$dom.T, na.rm=T),3)-0.2, 
                                                             round(max(sig.all$dom.T, na.rm=T),3)+0.2, by=0.3))+ 
  labs(x="Dominance effect Z-score")+ 
  theme(axis.title = element_text(size=9))  
dom_degree <- ggplot(sig.all, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(round(min(sig.all$degree_of_dominance, na.rm=T),3)-0.05, round(max(sig.all$degree_of_dominance, na.rm=T),3)+0.05, by=0.05))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9))  

(add+dom) /dom_degree+ plot_annotation(title="All RNA", tag_levels = "A")
ggsave("RNA_zscores_degree_sep_pval_corr.pdf")
melted_sig.all <- reshape2::melt(sig.all, id.vars="marker",measure.vars=c("add.T", "dom.T"))
ad<-ggplot(melted_sig.all, aes(x=value))+ geom_histogram(aes(color = variable, fill = variable),
                                                      alpha = 0.2, position = "identity", bins=50) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + labs(x="z-score") +theme_test()
dom_dens <- ggplot(sig.all, aes(x=degree_of_dominance))+ geom_density(aes(y = ..count..), fill = "lightgray") +
  geom_vline(aes(xintercept = 1), 
             linetype = "dashed", size = 0.6,
             color = "#FC4E07") +theme_test()
ad/dom_dens
ggsave("RNA_zscores_degree_pval_corr.pdf")

sig.all_R <- sig.all
sig.all_R$DNAorRNA <- "RNA"
sig.all_D$DNAorRNA <- "DNA"
sig.all_all <- rbind(sig.all_D, sig.all_R)
save(sig.all_all, file="sig.all_all_pval_corr.Rdata")
melted_sig.all <- reshape2::melt(sig.all_all, id.vars="marker",measure.vars=c("add.T", "dom.T"))
ad<-ggplot(melted_sig.all, aes(x=value))+ geom_histogram(aes(color = variable, fill = variable),
                                                      alpha = 0.2, position = "identity", bins=80) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + labs(x="z-score") +theme_test()
dom_deg_all <- ggplot(sig.all_all, aes(x=degree_of_dominance))+ 
  geom_histogram(colour="black", fill="white",
                 breaks=seq(round(min(sig.all_all$degree_of_dominance, na.rm=T),3)-0.05, 
                                                                                    round(max(sig.all_all$degree_of_dominance, na.rm=T),3)+0.05, by=0.05))+ 
  labs(x="Degree of dominance (d/a)") + geom_vline(xintercept = 1, colour="red", linetype=2)+theme_test()
  
ad/dom_deg_all
ggsave("RNA_DNA_zscores_degree_pval_corr.pdf")

##### good plot #### 
sig.all_all$degree_of_dominance<- sig.all_all$dom.T/abs(sig.all_all$add.T)
sig.all_all$deg_beta <- sig.all_all$dom.Beta/abs(sig.all_all$add.Beta)
sig.all_all$cat<- cut(sig.all_all$degree_of_dominance, breaks=c(-Inf,-1.25,-0.75,-0.25,0.25,0.75,1.25,Inf), 
                   labels = c("underdominant","recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant"))
sig.all_cat<-sig.all_all %>% 
  group_by(cat) %>% 
  tally() %>% 
  mutate(percent=(n/1372)*100) %>% 
  drop_na%>% 
  ggplot(aes(x=cat, y=n))+geom_bar(stat="identity") + geom_text(aes(label=paste(round(percent,2), "%")), vjust=-0.5) +
  labs(x="", y="Frequency")
sig.all_cat
ggsave("dominance_categories_pval_corr.pdf")

ggplot(sig.all_all, aes(x=abs(add.Beta), y=deg_beta))+geom_point()
ggplot(sig.all_all, aes(x=dom.Beta, y=add.Beta))+geom_point() +geom_smooth(method="lm")
ggsave("dom.beta_add.beta_pval_corr.pdf")
ggplot(sig.all_all, aes(x=dom.Beta, y=abs(add.Beta)))+geom_point() +geom_smooth(method="lm")
ggsave("dom.beta_abs_add.beta_pval_corr.pdf")

ggplot(sig.all_all, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(-0.75,1.25, 0.5))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9)) +scale_x_continuous(breaks=c(-0.75,-0.25,0.25,0.75,1.25))
ggsave("dom_histogram_numbers_pval_corr.pdf")
ggplot(sig.all_all, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=c(-Inf, seq(-1.25,1.25, 0.5), Inf))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9)) +
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1), 
                                                               labels = c("recessive","partially recessive", "additive", "partially dominant", "dominant"))
ggsave("dom_histogram_labels_pval_corr.pdf")

ggplot(sig.all_all, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(-2,2, by=0.03))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9)) 

sig.all_all$bins<- cut(sig.all_all$degree_of_dominance, breaks=c(-Inf,seq(-2,2,0.05),Inf), 
                    labels = c("<-2", seq(-2,2,0.05)))
sig.all_all$cat<- cut(sig.all_all$degree_of_dominance, breaks=c(-Inf,-1.25,-0.75,-0.25,0.25,0.75,1.25,Inf), 
                   labels = c("underdominant","recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant"))

sig.all_cat<-sig.all_all %>% 
  group_by(cat) %>% 
  tally() %>% 
  mutate(percent=(n/1295)*100) %>% 
  drop_na %>% mutate(start=c(-Inf, -1.25,-0.75,0.75), stop=c(-1.25,-0.75,-0.25,1.25), 
                     middle=c(-1.4, -1,-0.5,1)) # there are no partially dominant and overdominant 

sig.all_sum<-sig.all_all %>% 
  group_by(bins) %>% 
  tally() %>% 
  drop_na 

# ggplot()+geom_point(data = sig.all_sum, mapping = aes(x=as.numeric(as.character(bins)), y=n))+
#   scale_x_continuous()+geom_rect(
#     aes(xmin = start, xmax = stop, fill = cat), 
#     ymin = -Inf, ymax = Inf, alpha = 0.2, 
#     data = sig.all_cat
#   )

ggplot()+
  geom_segment(data = sig.all_sum, aes(x=as.numeric(as.character(bins)), 
                                    xend=as.numeric(as.character(bins)), y=0, yend=n), size=5, alpha=0.9) +
  scale_x_continuous(breaks = c(-Inf,-2, -1.25,-0.75,-0.25,0.25,0.75,1.25,1.99, Inf), 
                     labels = c(-Inf,-2, -1.25,-0.75,-0.25,0.25,0.75,1.25,">2", Inf))+
  geom_rect(data = sig.all_cat,aes(xmin = start, xmax = stop, fill = cat), 
            ymin = -Inf, ymax = Inf, alpha = 0.2)+ 
  geom_text(data=sig.all_cat, mapping=aes(x=middle, y=305,label=paste(round(percent,2), "%")), vjust=-0.5, size=3.5) +
  geom_text(data=sig.all_cat, mapping=aes(x=middle, y=295,label=cat), vjust=-0.5, size=3)+
  labs(x=expression(paste("Dominance", italic( "(d/|a|)"))), y="Frequency", fill="")
ggsave("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/dominance_classified_pval_corr.pdf")

# dom - mus ####
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/snps.anc.inf.info.Rdata")
sig.all_all_minmaj<- read_delim("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Plots/Effect_plots/allP_all_minmaj.csv", delim=";")
sig.all_all_minmaj <- sig.all_all_minmaj[which(sig.all_all_minmaj$marker%in%sig.all_all$marker ),]
snps.anc.inf.info_without_X<- snps.anc.inf.info %>% 
  filter(chr.num<20)
tot_inf <- snps.anc.inf.info %>% 
  rownames_to_column("marker") %>% 
  inner_join(sig.all_all_minmaj, by="marker" ) %>% 
  mutate(dom.T=as.numeric(as.character(dom.Beta))/as.numeric(as.character(dom.StdErr))) %>% 
  mutate(degree_of_dominance=dom.T/abs(as.numeric(as.character(add.T)))) %>% 
  dplyr::select(marker, chr.num, anc.inf, dom.allele, mus.allele, A1.x, A2.x, 
                AA_min, BB_maj,tax, add.T, dom.T, degree_of_dominance, high )


tot_inf$cat<- cut(tot_inf$degree_of_dominance, breaks=c(-Inf,-1.25,-0.75,-0.25,0.25,0.75,1.25,Inf), 
                   labels = c("underdominant","recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant"))


tot_inf2 <- tot_inf %>% 
  mutate(min.allele=ifelse(AA_min==A1.x, "A", ifelse(AA_min==A2.x, "B", NA)), 
         maj.allele=ifelse(BB_maj==A1.x, "A", ifelse(BB_maj==A2.x, "B", NA))) %>% 
  mutate(min.allele.sp=ifelse(dom.allele==min.allele, "dom", ifelse(mus.allele==min.allele, "mus", NA))) %>% 
  mutate(maj.allele.sp=ifelse(dom.allele==maj.allele, "dom", ifelse(mus.allele==maj.allele, "mus", NA))) %>% 
  mutate(high.allele=ifelse(high=="minor_allele", min.allele.sp, ifelse(high=="major_allele", 
                                                                        maj.allele.sp, 
                                                                        ifelse(high=="H", cat, "NA"))))

# tot_inf2$high.allele<- c("underdominant", "overdominant", "dom", "mus", "other")

tot_inf_tax <- tot_inf2 %>% 
  group_by(tax, high.allele) %>% 
  drop_na() %>% 
  tally()

ggplot(tot_inf_tax, aes(x=fct_reorder(tax, n), y=n, fill=high.allele))+geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle=90))+coord_flip()+theme_test()+
  scale_y_continuous(limits = c(0,35), expand = c(0, 0))

tot_inf_cat <- tot_inf2 %>% 
  group_by(tax, high.allele, cat) %>% 
  drop_na() %>% 
  tally()
ggplot(tot_inf_cat, aes(x=fct_reorder(tax, n), y=n, fill=cat))+geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle=90))+coord_flip()+theme_test()+facet_wrap(~high.allele, scales = "free_y")


tot_inf_sum <- tot_inf2 %>% 
  group_by(high.allele) %>% 
  tally() 
tot_inf_sum$high.allele<- c("heterozygote overdominant", "dom", "mus", "unknown")

pie <-ggplot(tot_inf_sum, aes(x="", y=n, fill=factor(high.allele))) +
  geom_bar(stat="identity", width=1, color="white", alpha=0.8) +
  coord_polar("y", start=0) +
  
  theme_void() +# remove background, grid, numeric labels 
scale_fill_d3()

library(scales)
pie + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
                position = position_stack(vjust = 0.5)) +labs(fill="High allele")


ggsave("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/high.allele_mus_dom_pval_corr.pdf")


### all markers with non-informative 
tot_inf_all <- snps.anc.inf.info %>% 
  rownames_to_column("marker") %>% 
  right_join(sig.all_all_minmaj, by="marker" ) %>% 
  mutate(dom.T=as.numeric(as.character(dom.Beta))/as.numeric(as.character(dom.StdErr))) %>% 
  mutate(degree_of_dominance=dom.T/abs(as.numeric(as.character(add.T)))) %>% 
  dplyr::select(marker, chr.num, anc.inf, dom.allele, mus.allele, A1.x, A2.x, 
                AA_min, BB_maj,tax, add.T, dom.T, degree_of_dominance, high ) 
tot_inf_all$cat<- cut(tot_inf_all$degree_of_dominance, breaks=c(-Inf,-1.25,-0.75,-0.25,0.25,0.75,1.25,Inf), 
                  labels = c("underdominant","recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant"))


tot_inf_all2 <- tot_inf_all %>% 
  mutate(min.allele=ifelse(AA_min==A1.x, "A", ifelse(AA_min==A2.x, "B", NA)), 
         maj.allele=ifelse(BB_maj==A1.x, "A", ifelse(BB_maj==A2.x, "B", NA))) %>% 
  mutate(min.allele.sp=ifelse(is.na(dom.allele),"non-informative",
                              ifelse(dom.allele==min.allele, "dom", 
                                     ifelse(mus.allele==min.allele, "mus", NA)))) %>% 
  mutate(maj.allele.sp=ifelse(is.na(dom.allele),"non-informative",
                              ifelse(dom.allele==maj.allele, "dom", 
                                     ifelse(mus.allele==maj.allele, "mus", NA)))) %>% 
  mutate(high.allele=ifelse(is.na(dom.allele),as.character(cat),
                            ifelse(high=="minor_allele", min.allele.sp, 
                                   ifelse(high=="major_allele", maj.allele.sp,
                                          ifelse(high=="H", as.character(cat), "NA"))))) %>% 
  mutate(high.allele2=ifelse(is.na(dom.allele),"not-informative",
                            ifelse(high=="minor_allele", min.allele.sp, 
                                   ifelse(high=="major_allele", maj.allele.sp,
                                          ifelse(high=="H", as.character(cat), "NA")))))

tot_inf_tax_all <- tot_inf_all2 %>% 
  group_by(tax, high.allele2) %>% 
  tally()
tot_inf_tax_all_sum <- tot_inf_all2 %>% 
  group_by(tax) %>% 
  tally() %>% arrange(desc(n))

ggplot(tot_inf_tax_all, aes(x=factor(tax, levels=tot_inf_tax_all_sum$tax), y=n, fill=factor(high.allele2)))+
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle=90))+coord_flip()+theme_test()+
  scale_y_continuous( expand = c(0, 0))+scale_fill_d3()+labs(x="", y="Number of significant SNPs", fill="high allele")
ggsave("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/high.allele_mus_dom_all_tax_pval_corr.pdf", height = 15)

tot_inf_sum_all <- tot_inf_all2 %>% 
  group_by(high.allele2) %>% 
  tally() 
tot_inf_sum_all[5,1]<- "unknown"



pie_all <-ggplot(tot_inf_sum_all, aes(x="", y=n, fill=factor(high.allele2), color=factor(high.allele2))) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  
  theme_void() +# remove background, grid, numeric labels 
  scale_fill_d3()

library(scales)
pie_all + geom_text(aes(label = paste(round(n / sum(n) * 100, 1), "%")),
                position = position_stack(vjust = 0.5)) +labs(fill="High allele")


ggsave("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/high.allele_mus_dom_with_noninf_pval_corr.pdf")

# Let the musculus  allele be the major allele ==1
tot_inf2<- tot_inf2 %>% 
  mutate(deg_Dom_MD=ifelse(min.allele.sp=="mus", degree_of_dominance *-1, degree_of_dominance))
tot_inf2$cat_MD<- cut(tot_inf2$deg_Dom_MD, breaks=c(-Inf,-1.25,-0.75,-0.25,0.25,0.75,1.25,Inf), 
                  labels = c("underdominant","recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant"))
tot_inf2$bins_MD<- cut(tot_inf2$deg_Dom_MD, breaks=c(-Inf,seq(-2,2,0.05),Inf), 
                    labels = c("<-2", seq(-2,2,0.05)))

MD.allele<-tot_inf2 %>% 
  group_by(cat_MD) %>% 
  tally() %>% 
  mutate(percent=(n/690)*100) %>% 
  drop_na %>% mutate(start=c(-1.25,-0.75,-0.25,0.25,0.75,1.25), stop=c(-0.75,-0.25,0.25,0.75,1.25, Inf), 
                     middle=c(-1,-0.5,0,0.5,1,1.75))

MD.allele_sum<-tot_inf2 %>% 
  group_by(bins_MD) %>% 
  tally() %>% 
  drop_na 


ggplot()+
  geom_segment(data = MD.allele_sum, aes(x=as.numeric(as.character(bins_MD)), 
                                    xend=as.numeric(as.character(bins_MD)), y=0, yend=n), size=3, alpha=0.9) +
  scale_x_continuous(breaks = c(-Inf,-2, -1.25,-0.75,-0.25,0.25,0.75,1.25,1.99, Inf), 
                     labels = c(-Inf,-2, -1.25,-0.75,-0.25,0.25,0.75,1.25,">2", Inf))+
  geom_rect(data = MD.allele,aes(xmin = start, xmax = stop, fill = cat_MD), 
            ymin = -Inf, ymax = Inf, alpha = 0.2)+ 
  geom_text(data=MD.allele, mapping=aes(x=middle, y=305,label=paste(round(percent,2), "%")), vjust=-0.5, size=3.5) +
  geom_text(data=MD.allele, mapping=aes(x=middle, y=295,label=cat_MD), vjust=-0.5, size=3)+
  labs(x=expression(paste("Dominance", italic( "(d/|a|)"))), y="Frequency", fill="")+
  ggtitle(expression(paste(italic("M. m. musculus "), "allele"))) +theme_test()
ggsave("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/musculus_allele_dominance.pdf")

sig.all_all %>%  drop_na(add.T, dom.T) %>% 
ggplot( aes(x=add.T, y=dom.T, color=cat)) +
  geom_hline(yintercept = 0, color="black", alpha=0.8)+
  geom_vline(xintercept = 0, color="black", alpha=0.8) +
  geom_abline(slope=1, intercept=-1.25,color="grey", alpha=0.8) +
  geom_abline(slope=1, intercept=1.25, color="grey", alpha=0.8) +
  geom_point() +theme_test() +scale_color_d3()+
  annotate("text", x=-9, y=9, label="minor allele dominance", color="#218380")+
  annotate("text", x=5, y=9, label="major allele dominance", color="#218380")+
  annotate("text", x=-10, y=0.75, label="pure additive", color="#218380")+
  annotate("text", x=0, y=13, label="overdominance", color="#218380")+ 
  labs(x="Additive effect (a)", y="Dominance effect (d)", color="Major allele")


ggsave("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/scatter_ad_dom_pval_corr.pdf")

tot_inf_all2 %>%  drop_na(add.T, dom.T) %>% 
  ggplot( aes(x=add.T, y=dom.T, color=high.allele2)) +
  geom_hline(yintercept = 0, color="black", alpha=0.8)+
  geom_vline(xintercept = 0, color="black", alpha=0.8) +
  geom_abline(slope=1, intercept=-1.25,color="grey", alpha=0.8) +
  geom_abline(slope=1, intercept=1.25, color="grey", alpha=0.8) +
  geom_point() +theme_test() +scale_color_d3()+
  annotate("text", x=-9, y=9, label="minor allele dominance", color="#218380")+
  annotate("text", x=5, y=9, label="major allele dominance", color="#218380")+
  annotate("text", x=-10, y=0.75, label="pure additive", color="#218380")+
  annotate("text", x=0, y=13, label="overdominance", color="#218380")+ 
  labs(x="Additive effect (a)", y="Dominance effect (d)", color="Major allele")

ggsave("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/scatter_ad_dom_by_mus_dom_pval_corr.pdf")

#### fisher's exact test ####
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata")
snps.auto <- snps %>% 
  filter(chr %in% c(1:19))
informative_clean_snps <- snps.anc.inf.info_without_X %>% 
  rownames_to_column("marker") %>% 
  inner_join(snps, by="marker")

length(unique(tot_inf$marker))

conting_tab <- tibble("inform"=c(398,17757), "non-inform"=c(384,14086))
ft <-fisher.test(conting_tab)
ft$estimate

conting_tab <- tibble("sig"=c(398,384), "non-sig"=c(17757,14086))
ft <-fisher.test(conting_tab)
ft$estimate
ft
fixed_clean_snps <- snps.anc.inf.info_without_X %>% 
  rownames_to_column("marker") %>% 
  inner_join(snps, by="marker") %>% filter(anc.inf=="F")

fixed_signif_snps<-inner_join(fixed_clean_snps,sig.all_all, by="marker")

conting_tab <- tibble("sig"=c(50,732), "non-sig"=c(2765,29078))
ft <-fisher.test(conting_tab)
ft

# How many dom ar maj allele #####
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/otu/SV1_gwscan_allchr.Rdata")
snps_minmaj <- tot_inf2 %>% 
  inner_join(SV1_gwscan, by="marker") %>% 
  inner_join(tot_inf, by="marker") 
times_mus_maj <- snps_minmaj %>% 
  distinct(marker, maj.allele.sp) %>% 
  group_by(maj.allele.sp) %>% 
  tally()

times_mus_high.allele <- snps_minmaj %>% 
  distinct(marker, high.allele) %>% 
  group_by(high.allele) %>% 
  tally()
save.image("z-scores-plot-image.Rdata")
save(tot_inf_all2, file="/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Plots/Effect_plots/tot_inf_all2.Rdata")  


# core prevalence
core <- read.table("../../../../Phenotyping_27.02.2020/tables_core/data_tot_core.txt")


# RNA 
library(ggpubr)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
RNA_core <- core %>%  filter(core_rna=="yes") %>% select(taxa, med_abund_rna, mean_abund_rna,prev_01_rna)
sig.all_core_RNA <- sig.all %>% left_join(RNA_core, by=c("tax"="taxa"))
ggplot(sig.all_core_RNA, aes(x=deg_beta, y=med_abund_rna))+ geom_point()+stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))

ggplotRegression(lm(deg_beta ~ mean_abund_rna, data = sig.all_core_RNA))
ggsave("RNA_mean_abundVSdeg_beta.pdf")
ggplotRegression(lm(deg_beta ~ med_abund_rna, data = sig.all_core_RNA))
ggsave("RNA_med_abundVSdeg_beta.pdf")
ggplotRegression(lm(deg_beta ~ prev_01_rna, data = sig.all_core_RNA))
ggsave("RNA_mean_abundVSprev01.pdf")
