setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/")
library(tidyverse)
load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/all_DNA.Rdata")
all_dna<- all_files
load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/all_RNA.Rdata")
all_rna <- all_files
all_rna$dna.rna<-"RNA"
all_dna$dna.rna <- "DNA"
all_dna_rna <- rbind(all_dna, all_rna)
save(all_dna_rna,file="/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_18032021_GW.Rdata")
 ####### GW #######
all_unique <- all_dna_rna %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_rna <- all_rna %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_dna <- all_dna %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
med_all<-median(all_unique$stop.LD.pos-all_unique$start.LD.pos)/1e6
med_rna<-median(all_uniq_rna$stop.LD.pos-all_uniq_rna$start.LD.pos)/1e6
med_dna<-median(all_uniq_dna$stop.LD.pos-all_uniq_dna$start.LD.pos)/1e6
all_uniq_dna$length <- (all_uniq_dna$stop.LD.pos-all_uniq_dna$start.LD.pos)/1e6
all_uniq_rna$length <- (all_uniq_rna$stop.LD.pos-all_uniq_rna$start.LD.pos)/1e6
all_unique$length <- (all_unique$stop.LD.pos-all_unique$start.LD.pos)/1e6

med0_all<- median(all_unique$length[all_unique$length>0])
med0_rna <- median(all_uniq_rna$length[all_uniq_rna$length>0])
med0_dna <- median(all_uniq_dna$length[all_uniq_dna$length>0])

minP_rna <- all_rna %>% 
  filter(P.type=="P") %>% 
  summarise(min(peak.P.type))
minP_dna <- all_dna %>% 
  filter(P.type=="P") %>% 
  summarise(min(peak.P.type))
minP_all <- all_dna_rna %>% 
  filter(P.type=="P") %>% 
  summarise(min(peak.P.type))


minaddP_rna <- all_rna %>% 
  filter(P.type=="add.P") %>% 
  summarise(min(peak.P.type))
minaddP_dna <- all_dna %>% 
  filter(P.type=="add.P") %>% 
  summarise(min(peak.P.type))
minaddP_all <- all_dna_rna %>% 
  filter(P.type=="add.P") %>% 
  summarise(min(peak.P.type))


mindomP_rna <- all_rna %>% 
  filter(P.type=="dom.P") %>% 
  summarise(min(peak.P.type))
mindomP_dna <- all_dna %>% 
  filter(P.type=="dom.P") %>% 
  summarise(min(peak.P.type))
mindomP_all <- all_dna_rna %>% 
  filter(P.type=="dom.P") %>% 
  summarise(min(peak.P.type))




all_unique.add <- all_dna_rna %>% 
  filter(P.type=="add.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_rna.add <- all_rna %>% 
  filter(P.type=="add.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_dna.add <- all_dna %>% 
  filter(P.type=="add.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)

all_unique.dom <- all_dna_rna %>% 
  filter(P.type=="dom.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_rna.dom <- all_rna %>% 
  filter(P.type=="dom.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_dna.dom <- all_dna %>% 
  filter(P.type=="dom.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)

all_unique.P <- all_dna_rna %>% 
  filter(P.type=="P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_rna.P <- all_rna %>% 
  filter(P.type=="P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_dna.P <- all_dna %>% 
  filter(P.type=="P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)

all_taxa_dna<- all_dna %>% 
  distinct(trait)
all_taxa_rna<- all_rna %>% 
  distinct(trait)
all_taxa<- all_dna_rna %>% 
  distinct(trait)
summary_numbers <- data.frame(matrix(ncol=12, nrow=3))
colnames(summary_numbers)<- c("Mapped.taxa", "Taxa.sig.hits", "Significant.hits", "Unique.sig.loci", "median.unique", "med.unique.0.removed","unique.P.loci", "unique.add.loci", "unique.dom.loci", "min.P", "min.addP", "min.domP")
rownames(summary_numbers)<- c("DNA", "RNA","Total")
summary_numbers$Mapped.taxa <- c(101,142,153)
summary_numbers$Taxa.sig.hits <- c(nrow(all_taxa_dna), nrow(all_taxa_rna), nrow(all_taxa))
summary_numbers$Significant.hits<- c(nrow(all_dna), nrow(all_rna), nrow(all_dna_rna))
summary_numbers$Unique.sig.loci<- c(nrow(all_uniq_dna), nrow(all_uniq_rna), nrow(all_unique))
summary_numbers$median.unique <- c(signif(med_dna,3), signif(med_rna, 3), signif(med_all, 3))
summary_numbers$med.unique.0.removed<- c(signif(med0_dna,3), signif(med0_rna, 3), signif(med0_all, 3))
summary_numbers$unique.add.loci <- c(nrow(all_uniq_dna.add), nrow(all_uniq_rna.add), nrow(all_unique.add))
summary_numbers$unique.dom.loci <- c(nrow(all_uniq_dna.dom), nrow(all_uniq_rna.dom), nrow(all_unique.dom))
summary_numbers$unique.P.loci <- c(nrow(all_uniq_dna.P), nrow(all_uniq_rna.P), nrow(all_unique.P))
summary_numbers$min.P <- c(as.character(signif(minP_dna,3)), as.character(signif(minP_rna, 3)), as.character(signif(minP_all,3)))
summary_numbers$min.addP <- c(as.character(signif(minaddP_dna,3)), as.character(signif(minaddP_rna, 3)), as.character(signif(minaddP_all,3)))
summary_numbers$min.domP <- c(as.character(signif(mindomP_dna,3)), as.character(signif(mindomP_rna, 3)), as.character(signif(mindomP_all,3)))
write_csv(summary_numbers, file="~/Documents/PhD/Experiments/Final_QTL_mapping/Results/summary_numbers_associations_GW.csv", quote=F)

uniq_GW <- all_dna_rna %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)

# number of QTL per taxon ####

loci_per_trait <- all_dna_rna %>% group_by(trait) %>%  tally()
tot_med_per_trait_all<-median(loci_per_trait$n)  
tot_mean_per_trait_all<-mean(loci_per_trait$n)

loci_per_trait_uniq <- all_dna_rna %>% group_by(trait) %>% distinct(chr,start.LD.pos,stop.LD.pos) %>%  tally()
tot_med_per_trait_all_uniq<-median(loci_per_trait_uniq$n)
tot_mean_per_trait_all_uniq<-mean(loci_per_trait_uniq$n)

loci_per_trait_rna <- all_rna %>% group_by(trait) %>%  tally()
tot_med_per_trait_rna<-median(loci_per_trait_rna$n)  
tot_mean_per_trait_rna<-mean(loci_per_trait_rna$n)

loci_per_trait_uniq_rna <- all_rna %>% group_by(trait) %>% distinct(chr,start.LD.pos,stop.LD.pos) %>%  tally()
tot_med_per_trait_all_uniq_rna<-median(loci_per_trait_uniq_rna$n)
tot_mean_per_trait_all_uniq_rna<-mean(loci_per_trait_uniq_rna$n)

loci_per_trait_dna <- all_dna %>% group_by(trait) %>%  tally()
tot_med_per_trait_dna<-median(loci_per_trait_dna$n)  
tot_mean_per_trait_dna<-mean(loci_per_trait_dna$n)

loci_per_trait_uniq_dna <- all_dna %>% group_by(trait) %>% distinct(chr,start.LD.pos,stop.LD.pos) %>%  tally()
tot_med_per_trait_all_uniq_dna<-median(loci_per_trait_uniq_dna$n)
tot_mean_per_trait_all_uniq_dna<-mean(loci_per_trait_uniq_dna$n)

summary_numbers$mean.loci.per.trait <- c(tot_mean_per_trait_dna, tot_mean_per_trait_rna, tot_mean_per_trait_all)
summary_numbers$median.loci.per.trait <- c(tot_med_per_trait_dna, tot_med_per_trait_rna, tot_med_per_trait_all)
summary_numbers$mean.loci.per.trait.unique <- c(tot_mean_per_trait_all_uniq_dna, tot_mean_per_trait_all_uniq_rna, tot_mean_per_trait_all_uniq)
summary_numbers$median.loci.per.trait.unique <- c(tot_med_per_trait_all_uniq_dna, tot_med_per_trait_all_uniq_rna, tot_med_per_trait_all_uniq)

# genes per locus ####

mean_genes_per_locus <- all_dna_rna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_genes))
med_genes_per_locus <- all_dna_rna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_genes))

mean_PC_genes_per_locus <- all_dna_rna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_protein_coding_genes))
med_PC_genes_per_locus <- all_dna_rna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_protein_coding_genes))

mean_genes_per_locus_dna <- all_dna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_genes))
med_genes_per_locus_dna <- all_dna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_genes))

mean_PC_genes_per_locus_dna <- all_dna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_protein_coding_genes))
med_PC_genes_per_locus_dna <- all_dna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_protein_coding_genes))

mean_genes_per_locus_rna <- all_rna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_genes))
med_genes_per_locus_rna <- all_rna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_genes))

mean_PC_genes_per_locus_rna <- all_rna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_protein_coding_genes))
med_PC_genes_per_locus_rna <- all_rna %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_protein_coding_genes))

summary_numbers$mean.genes.per.locus <- c(mean_genes_per_locus_dna$`mean(total_genes)`,mean_genes_per_locus_rna$`mean(total_genes)`,mean_genes_per_locus$`mean(total_genes)`)
summary_numbers$median.genes.per.locus <- c(med_genes_per_locus_dna$`median(total_genes)`,med_genes_per_locus_rna$`median(total_genes)`, med_genes_per_locus$`median(total_genes)`)

summary_numbers$mean.PC.genes.per.locus <- c(mean_PC_genes_per_locus_dna$`mean(total_protein_coding_genes)`,mean_PC_genes_per_locus_rna$`mean(total_protein_coding_genes)`, mean_PC_genes_per_locus$`mean(total_protein_coding_genes)`)
summary_numbers$median.PC.genes.per.locus<- c(med_PC_genes_per_locus_dna$`median(total_protein_coding_genes)`, med_PC_genes_per_locus_rna$`median(total_protein_coding_genes)`, med_PC_genes_per_locus$`median(total_protein_coding_genes)`)
sum_table_to_write<- as.data.frame(t(summary_numbers))
rownames(sum_table_to_write)<- colnames(summary_numbers)
sum_table_to_write <- sum_table_to_write %>% rownames_to_column("Variable")
write_delim(sum_table_to_write,"summary_numbers_associations_GW.csv", delim=";", col_names = T)

####### SW #######
SW_threshold =0.05/32625/118.2177

all_dna_SW <- all_dna %>% 
  filter(peak.P.type<SW_threshold)
all_rna_SW <- all_rna %>% 
  filter(peak.P.type<SW_threshold)
all_dna_rna_SW <- all_dna_rna %>% 
  filter(peak.P.type<SW_threshold)
save(all_dna_rna_SW,file="/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_18032021_SW.Rdata")

all_unique <- all_dna_rna_SW %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_rna <- all_rna_SW%>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_dna <- all_dna_SW%>% 
  distinct(chr, start.LD.pos, stop.LD.pos)

med_all<-median(all_dna_rna_SW$stop.LD.pos-all_dna_rna_SW$start.LD.pos)/1e6
med_rna<-median(all_rna_SW$stop.LD.pos-all_rna_SW$start.LD.pos)/1e6
med_dna<-median(all_dna_SW$stop.LD.pos-all_dna_SW$start.LD.pos)/1e6

med0_all<- median(all_dna_rna_SW$length[all_dna_rna_SW$length>0])
med0_rna <- median(all_rna_SW$length[all_rna_SW$length>0])
med0_dna <- median(all_dna_SW$length[all_dna_SW$length>0])

med_all.uniq<-median(all_unique$stop.LD.pos-all_unique$start.LD.pos)/1e6
med_rna.uniq<-median(all_uniq_rna$stop.LD.pos-all_uniq_rna$start.LD.pos)/1e6
med_dna.uniq<-median(all_uniq_dna$stop.LD.pos-all_uniq_dna$start.LD.pos)/1e6
all_uniq_dna$length <- (all_uniq_dna$stop.LD.pos-all_uniq_dna$start.LD.pos)/1e6
all_uniq_rna$length <- (all_uniq_rna$stop.LD.pos-all_uniq_rna$start.LD.pos)/1e6
all_unique$length <- (all_unique$stop.LD.pos-all_unique$start.LD.pos)/1e6

med0_all.uniq<- median(all_unique$length[all_unique$length>0])
med0_rna.uniq <- median(all_uniq_rna$length[all_uniq_rna$length>0])
med0_dna.uniq<- median(all_uniq_dna$length[all_uniq_dna$length>0])

minP_rna <- all_rna_SW%>% 
  filter(P.type=="P") %>% 
  summarise(min(peak.P.type))
minP_dna <- all_dna_SW%>% 
  filter(P.type=="P") %>% 
  summarise(min(peak.P.type))
minP_all <- all_dna_rna_SW %>% 
  filter(P.type=="P") %>% 
  summarise(min(peak.P.type))


minaddP_rna <- all_rna_SW%>% 
  filter(P.type=="add.P") %>% 
  summarise(min(peak.P.type))
minaddP_dna <- all_dna_SW%>% 
  filter(P.type=="add.P") %>% 
  summarise(min(peak.P.type))
minaddP_all <- all_dna_rna_SW %>% 
  filter(P.type=="add.P") %>% 
  summarise(min(peak.P.type))


mindomP_rna <- all_rna_SW%>% 
  filter(P.type=="dom.P") %>% 
  summarise(min(peak.P.type))
mindomP_dna <- all_dna_SW%>% 
  filter(P.type=="dom.P") %>% 
  summarise(min(peak.P.type))
mindomP_all <- all_dna_rna_SW %>% 
  filter(P.type=="dom.P") %>% 
  summarise(min(peak.P.type))




all_unique.add <- all_dna_rna_SW %>% 
  filter(P.type=="add.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_rna.add <- all_rna_SW%>% 
  filter(P.type=="add.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_dna.add <- all_dna_SW%>% 
  filter(P.type=="add.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)

all_unique.dom <- all_dna_rna_SW %>% 
  filter(P.type=="dom.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_rna.dom <- all_rna_SW%>% 
  filter(P.type=="dom.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_dna.dom <- all_dna_SW%>% 
  filter(P.type=="dom.P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)

all_unique.P <- all_dna_rna_SW %>% 
  filter(P.type=="P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_rna.P <- all_rna_SW%>% 
  filter(P.type=="P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)
all_uniq_dna.P <- all_dna_SW%>% 
  filter(P.type=="P") %>% 
  distinct(chr, start.LD.pos, stop.LD.pos)

all_taxa_dna<- all_dna_SW%>% 
  distinct(trait)
all_taxa_rna<- all_rna_SW%>% 
  distinct(trait)
all_taxa<- all_dna_rna_SW %>% 
  distinct(trait)
summary_numbers <- data.frame(matrix(ncol=12, nrow=3))
colnames(summary_numbers)<- c("Mapped.taxa", "Taxa.sig.hits", "Significant.hits", "Unique.sig.loci", "median.unique", "med.unique.0.removed","unique.P.loci", "unique.add.loci", "unique.dom.loci", "min.P", "min.addP", "min.domP")
rownames(summary_numbers)<- c("DNA", "RNA","Total")
summary_numbers$Mapped.taxa <- c(101,142,153)
summary_numbers$Taxa.sig.hits <- c(nrow(all_taxa_dna), nrow(all_taxa_rna), nrow(all_taxa))
summary_numbers$Significant.hits<- c(nrow(all_dna_SW), nrow(all_rna_SW), nrow(all_dna_rna_SW))
summary_numbers$Unique.sig.loci<- c(nrow(all_uniq_dna), nrow(all_uniq_rna), nrow(all_unique))
summary_numbers$median.all <- c(signif(med_dna,3), signif(med_rna, 3), signif(med_all, 3))
summary_numbers$med.all.0.removed<- c(signif(med0_dna,3), signif(med0_rna, 3), signif(med0_all, 3))
summary_numbers$median.unique <- c(signif(med_dna.uniq,3), signif(med_rna.uniq, 3), signif(med_all.uniq, 3))
summary_numbers$med.unique.0.removed<- c(signif(med0_dna.uniq,3), signif(med0_rna.uniq, 3), signif(med0_all.uniq, 3))
summary_numbers$unique.add.loci <- c(nrow(all_uniq_dna.add), nrow(all_uniq_rna.add), nrow(all_unique.add))
summary_numbers$unique.dom.loci <- c(nrow(all_uniq_dna.dom), nrow(all_uniq_rna.dom), nrow(all_unique.dom))
summary_numbers$unique.P.loci <- c(nrow(all_uniq_dna.P), nrow(all_uniq_rna.P), nrow(all_unique.P))
summary_numbers$min.P <- c(as.character(signif(minP_dna,3)), as.character(signif(minP_rna, 3)), as.character(signif(minP_all,3)))
summary_numbers$min.addP <- c(as.character(signif(minaddP_dna,3)), as.character(signif(minaddP_rna, 3)), as.character(signif(minaddP_all,3)))
summary_numbers$min.domP <- c(as.character(signif(mindomP_dna,3)), as.character(signif(mindomP_rna, 3)), as.character(signif(mindomP_all,3)))

# load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/genes_DNA_RNA.Rdata")
# uniq_GW <- genes_DNA_RNA %>% 
#   distinct(chr, start, stop)

# number of QTL per taxon ####
#save.image("summary_numbers_workspace_SW.Rdata")

loci_per_trait <- all_dna_rna_SW %>% group_by(trait) %>%  tally()
tot_med_per_trait_all<-median(loci_per_trait$n)  
tot_mean_per_trait_all<-mean(loci_per_trait$n)

loci_per_trait_uniq <- all_dna_rna_SW %>% group_by(trait) %>% distinct(chr,start.LD.pos,stop.LD.pos) %>%  tally()
tot_med_per_trait_all_uniq<-median(loci_per_trait_uniq$n)
tot_mean_per_trait_all_uniq<-mean(loci_per_trait_uniq$n)

loci_per_trait_rna <- all_rna_SW%>% group_by(trait) %>%  tally()
tot_med_per_trait_rna<-median(loci_per_trait_rna$n)  
tot_mean_per_trait_rna<-mean(loci_per_trait_rna$n)

loci_per_trait_uniq_rna <- all_rna_SW%>% group_by(trait) %>% distinct(chr,start.LD.pos,stop.LD.pos) %>%  tally()
tot_med_per_trait_all_uniq_rna<-median(loci_per_trait_uniq_rna$n)
tot_mean_per_trait_all_uniq_rna<-mean(loci_per_trait_uniq_rna$n)

loci_per_trait_dna <- all_dna_SW%>% group_by(trait) %>%  tally()
tot_med_per_trait_dna<-median(loci_per_trait_dna$n)  
tot_mean_per_trait_dna<-mean(loci_per_trait_dna$n)

loci_per_trait_uniq_dna <- all_dna_SW%>% group_by(trait) %>% distinct(chr,start.LD.pos,stop.LD.pos) %>%  tally()
tot_med_per_trait_all_uniq_dna<-median(loci_per_trait_uniq_dna$n)
tot_mean_per_trait_all_uniq_dna<-mean(loci_per_trait_uniq_dna$n)

write.table(loci_per_trait_uniq_dna, "sig_loci_per_trait_dna.txt")
write.table(loci_per_trait_uniq_rna, "sig_loci_per_trait_rna.txt")

summary_numbers$mean.loci.per.trait <- c(signif(tot_mean_per_trait_dna,2), signif(tot_mean_per_trait_rna,2), signif(tot_mean_per_trait_all,2))
summary_numbers$median.loci.per.trait <- c(tot_med_per_trait_dna, tot_med_per_trait_rna, tot_med_per_trait_all)
summary_numbers$mean.loci.per.trait.unique <- c(signif(tot_mean_per_trait_all_uniq_dna,2), signif(tot_mean_per_trait_all_uniq_rna,2), signif(tot_mean_per_trait_all_uniq,2))
summary_numbers$median.loci.per.trait.unique <- c(tot_med_per_trait_all_uniq_dna, tot_med_per_trait_all_uniq_rna, tot_med_per_trait_all_uniq)

# genes per locus ####

mean_genes_per_locus <- all_dna_rna_SW %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_genes))
med_genes_per_locus <- all_dna_rna_SW %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_genes))

mean_PC_genes_per_locus <- all_dna_rna_SW %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_protein_coding_genes))
med_PC_genes_per_locus <- all_dna_rna_SW %>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_protein_coding_genes))

mean_genes_per_locus_dna <- all_dna_SW%>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_genes))
med_genes_per_locus_dna <- all_dna_SW%>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_genes))

mean_PC_genes_per_locus_dna <- all_dna_SW%>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_protein_coding_genes))
med_PC_genes_per_locus_dna <- all_dna_SW%>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_protein_coding_genes))

mean_genes_per_locus_rna <- all_rna_SW%>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_genes))
med_genes_per_locus_rna <- all_rna_SW%>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_genes))

mean_PC_genes_per_locus_rna <- all_rna_SW%>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(mean(total_protein_coding_genes))
med_PC_genes_per_locus_rna <- all_rna_SW%>% distinct(chr,start.LD.pos,stop.LD.pos,.keep_all = T) %>% 
  summarise(median(total_protein_coding_genes))

summary_numbers$mean.genes.per.locus <- c(signif(mean_genes_per_locus_dna$`mean(total_genes)`,2),signif(mean_genes_per_locus_rna$`mean(total_genes)`,2),signif(mean_genes_per_locus$`mean(total_genes)`,2))
summary_numbers$median.genes.per.locus <- c(med_genes_per_locus_dna$`median(total_genes)`,med_genes_per_locus_rna$`median(total_genes)`, med_genes_per_locus$`median(total_genes)`)

summary_numbers$mean.PC.genes.per.locus <- c(signif(mean_PC_genes_per_locus_dna$`mean(total_protein_coding_genes)`,2),signif(mean_PC_genes_per_locus_rna$`mean(total_protein_coding_genes)`,2), 
                                             signif(mean_PC_genes_per_locus$`mean(total_protein_coding_genes)`,2))
summary_numbers$median.PC.genes.per.locus<- c(med_PC_genes_per_locus_dna$`median(total_protein_coding_genes)`, med_PC_genes_per_locus_rna$`median(total_protein_coding_genes)`, med_PC_genes_per_locus$`median(total_protein_coding_genes)`)
sum_table_to_write<- as.data.frame(t(summary_numbers))
rownames(sum_table_to_write)<- colnames(summary_numbers)
sum_table_to_write <- sum_table_to_write %>% rownames_to_column("Variable")
write_delim(sum_table_to_write,"summary_numbers_associations_SW.csv", delim=";", col_names = T)


