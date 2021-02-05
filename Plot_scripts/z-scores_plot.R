setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/")
library(tidyverse)
library(patchwork)
library(ggsci)
add_res <- read.table("ALL_sig_snps_add_DNA_P.txt")
dom_res <- read.table("ALL_sig_snps_dom_P_DNA.txt")
P_res <- read.table("ALL_sig_snps_P_DNA.txt")
allP <- read.table("ALL_sig_snps_allP_DNA.txt")

allP$degree_of_dominance <- allP$dom.T/abs(allP$add.T)
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

ggplot(melted_allP, aes(x=fct_reorder(marker, value), y=value))+geom_bar(aes(fill=variable), stat = "identity", position=position_dodge())+
  theme(axis.text.x = element_blank())


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

##### good plot #### 
allP_all$degree_of_dominance<- allP_all$dom.T/abs(allP_all$add.T)
allP_all$cat<- cut(allP_all$degree_of_dominance, breaks=c(-Inf,-1.25,-0.75,-0.25,0.25,0.75,1.25,Inf), 
                   labels = c("underdominant","recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant"))
allP_cat<-allP_all %>% 
  group_by(cat) %>% 
  tally() %>% 
  mutate(percent=(n/1372)*100) %>% 
  drop_na%>% 
  ggplot(aes(x=cat, y=n))+geom_bar(stat="identity") + geom_text(aes(label=paste(round(percent,2), "%")), vjust=-0.5) +
  labs(x="", y="Frequency")
ggsave("dominance_categories.pdf")

ggplot(allP_all, aes(x=abs(add.T), y=degree_of_dominance))+geom_point()
ggplot(allP_all, aes(x=dom.T, y=add.T))+geom_point() +geom_smooth(method="lm")
ggplot(allP_all, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(-0.75,12, 0.5))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9)) +scale_x_continuous(breaks=c(-0.75,-0.25,0.25,0.75,1.25,12))
ggsave("dom_histogram_numbers.pdf")
ggplot(allP_all, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(-0.75,1.5, 0.5, Inf))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9)) +scale_x_continuous(breaks=c(-0.5,0,0.5,1), labels = c("partially recessive", "additive", "partially dominant", "dominant"))
ggsave("dom_histogram_labels.pdf")

ggplot(allP_all, aes(x=degree_of_dominance)) +
  geom_histogram(colour="black", fill="grey",
                 breaks=seq(-2,2, by=0.03))+ 
  labs(x="Degree of dominance (d/a)")+ 
  theme(axis.title = element_text(size=9)) 

allP_all$bins<- cut(allP_all$degree_of_dominance, breaks=c(-Inf,seq(-2,2,0.05),Inf), 
                    labels = c("<-2", seq(-2,2,0.05)))
allP_all$cat<- cut(allP_all$degree_of_dominance, breaks=c(-Inf,-1.25,-0.75,-0.25,0.25,0.75,1.25,Inf), 
                   labels = c("underdominant","recessive", "partially recessive", "additive", "partially dominant", "dominant", "overdominant"))

allP_cat<-allP_all %>% 
  group_by(cat) %>% 
  tally() %>% 
  mutate(percent=(n/1372)*100) %>% 
  drop_na %>% mutate(start=c(-Inf, -1.25,-0.75,-0.25,0.75,1.25), stop=c(-1.25,-0.75,-0.25,0.25,1.25, Inf), 
                     middle=c(-1.4, -1,-0.5,0,1,1.75))

allP_sum<-allP_all %>% 
  group_by(bins) %>% 
  tally() %>% 
  drop_na 

ggplot()+geom_point(data = allP_sum, mapping = aes(x=as.numeric(as.character(bins)), y=n))+
  scale_x_continuous()+geom_rect(
    aes(xmin = start, xmax = stop, fill = cat), 
    ymin = -Inf, ymax = Inf, alpha = 0.2, 
    data = allP_cat
  )

ggplot()+
  geom_segment(data = allP_sum, aes(x=as.numeric(as.character(bins)), 
                                    xend=as.numeric(as.character(bins)), y=0, yend=n), size=3, alpha=0.9) +
  scale_x_continuous(breaks = c(-Inf,-2, -1.25,-0.75,-0.25,0.25,0.75,1.25,1.99, Inf), 
                     labels = c(-Inf,-2, -1.25,-0.75,-0.25,0.25,0.75,1.25,">2", Inf))+
  geom_rect(data = allP_cat,aes(xmin = start, xmax = stop, fill = cat), 
            ymin = -Inf, ymax = Inf, alpha = 0.2)+ 
  geom_text(data=allP_cat, mapping=aes(x=middle, y=305,label=paste(round(percent,2), "%")), vjust=-0.5, size=3.5) +
  geom_text(data=allP_cat, mapping=aes(x=middle, y=295,label=cat), vjust=-0.5, size=3)+
  labs(x=expression(paste("Dominance", italic( "(d/|a|)"))), y="Frequency", fill="")
ggsave("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/dominance_classified.pdf")

# dom - mus ####
tot_inf <- snps.anc.inf.info %>% 
  rownames_to_column("marker") %>% 
  inner_join(allP_all_minmaj, by="marker" ) %>% 
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
  mutate(high.allele=ifelse(high=="minor_allele", min.allele.sp, ifelse(high=="major_allele", maj.allele.sp, 
                                                                        ifelse(high=="H", cat, "NA"))))

tot_inf_sum <- tot_inf2 %>% 
  group_by(high.allele) %>% 
  tally() 
tot_inf_sum$high.allele<- c("dominant", "overdominant", "dom", "mus", "Other")
ggplot(tot_inf_sum, aes(x="", y=n, fill=high.allele)) +
  geom_bar(stat="identity", width=1, color="white", alpha=0.8) +
  coord_polar("y", start=0) +
  
  theme_void() +# remove background, grid, numeric labels 
scale_fill_d3()
ggsave("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/high.allele_mus_dom.pdf")

