effect.plot <- function(taxa, trait, DNAorRNA, snp){
  require(readxl)
  require(reshape2)
  require(cowplot)
  require(ggpubr)
  require(tidyverse)
  kleurtjes <- c("#782ca8",
                 "#009e68",
                 "#ff415f")
  
  cbbPalette <- c("#1f78b4",  "#a6cee3","#b2df8a","#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", '#6a3d9a', "#ffff99", "#b15928", "#40e0d0", "#ffff00", "#ee82ee", "#fb8072", "#191970", "#fdb462", "#CBD4D4", "#fccde5", "#d9d9d9", "#bc80bd", "#32cd32", "#ffed6f", "lightcoral" )
  # load("~/Documents/PhD/Experiments/QTL_mapping_results/cleaning/clean_geno.Rdata")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_geno.Rdata")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata")
  load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Plots/Effect_plots/tot_inf_all2.Rdata")  
  theme_set(theme_test())
  clean_geno <- as.data.frame(clean_geno)
  colnames(clean_geno) <- gsub(x = colnames(clean_geno), pattern = "\\/", replacement = ".")  
  all_mark <-clean_geno[snp,]
  info_marker <- tot_inf_all2 %>% filter(marker==snp )
    a <- info_marker$AA_min
    b <- info_marker$BB_maj


  genotypes <- t(all_mark)
  
  consensus <-read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv", sep=",")
  allele_freq <-read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/allele_freq_consensus_mus_dom.csv", sep=",")
  rownames(consensus)<- consensus$X 
  consensus$X <- NULL
  mus_dom <- consensus[snp,]
  
  
  # abundances
  #DNAorRNA <- "DNA"
  # otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Miseq_hybrid_mice/DNA_and_RNA/tables_core/",trait, "_table_f2_core.csv"))
  otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/",trait, "_table_f2_core.csv"))
  rownames(otu_table)<- otu_table$X
  otu_table <- otu_table[which(grepl(paste0("_",DNAorRNA),rownames(otu_table))),] # only do the DNA
  rownames(otu_table) <- substr(rownames(otu_table), 1,15)
  
  SV <- as.data.frame(otu_table[,taxa])
  rownames(SV) <- rownames(otu_table)
  tot_with_SV <- merge(genotypes, SV, by="row.names")
  colnames(tot_with_SV)<- c("individuals", "marker", "abundance")
  tot_with_SV$counts <- tot_with_SV$abundance*10000
  tot_with_SV$log_counts <- log10(tot_with_SV$counts+1)
  tot_with_SV$marker<- as.factor(tot_with_SV$marker)
  
  # Summarizing data
  tot_with_SV_m <- tot_with_SV %>% 
    group_by(marker) %>% 
    summarise(mean = mean(abundance),
              se = psych::describe(abundance)$se, n=n())
  # mus dom alleles 
  if (info_marker$min.allele.sp!="non-informative"){
    min.allele <- info_marker$min.allele.sp
    maj.allele <- info_marker$maj.allele.sp
  } else {
    min.allele <- ""
    maj.allele <- ""
  }
  tot_with_SV_m$marker<- factor(tot_with_SV_m$marker, levels = c(a, "H", b), labels=c(paste(a, min.allele), "H", paste(b, maj.allele)))
  
  # Creating the plot
  fig <-tot_with_SV_m %>% 
    drop_na(marker) %>% 
    ggplot(aes(x = marker, 
               y = mean, 
               color = marker)) +

    geom_point(size=3) +
    geom_errorbar(aes(ymin = mean-1.96*se, 
                      ymax = mean+1.96*se),
                  width = .1)+ labs(caption = paste0("Mmm: ", mus_dom[,1], ", Mmd: ", mus_dom[,2]))+
    theme(legend.position = "none") + 
    geom_text(aes(label=paste("n=",n), y=max(mean+3*se)))+
    labs(x="Genotypes",y="Relative abundance")+ scale_color_manual(values = kleurtjes)+
    ggtitle(paste0(taxa," abundance for marker ", snp))+ scale_y_continuous(labels = scales::percent)
  
fig
  #annotate_figure(fig3, fig.lab = "GG: domesticus, AA: musculus", fig.lab.pos = "bottom.right", fig.lab.size = 12)
  ggsave(file=paste0(taxa,"_marker_", snp,"_",DNAorRNA, "_effect_plot.pdf"), plot=fig)

  fig2<-ggplot(tot_with_SV, aes(x=marker, y=log_counts)) + geom_boxplot()+ 
    # stat_compare_means(comparisons = my_comp, label = "p.signif")  
    theme(legend.position = "none") +  geom_jitter( aes(color=marker))+
    labs(x="Genotypes",y="Log10 transformed absolute abundance")+ scale_color_manual(values = cbbPalette)+
   
    ggtitle(paste0(taxa," abundance for marker ", snp)) + labs(caption = paste0("Mmm: ", mus_dom[,1], ", Mmd: ", mus_dom[,2]))
  ggsave(file=paste0(taxa,"_marker_", snp,"_",DNAorRNA, "log_box.pdf"), plot=fig2)
  return(fig)
  
}
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Plots/Effect_plots/")
# no AA
effect.plot("SV115", "otu", "RNA", "JAX00022875")


# H is high ###
effect.plot("SV264", "otu", "RNA", "UNCHS047355")
effect.plot("G_Odoribacter", "genus", "RNA", "UNC6535447")

effect.plot("SV91", "otu", "RNA", "UNCHS011150")

# Major allele high

effect.plot("SV91", "otu", "RNA", "UNCHS043099")
effect.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "JAX00499850")
effect.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNC6440308")

# min allele is high and H is the lowest

effect.plot("SV145", "otu", "DNA", "JAX00214148")
effect.plot("SV97", "otu", "RNA", "JAX00150093")
effect.plot("SV36", "otu", "RNA", "UNC21196986")
effect.plot("SV36", "otu", "RNA", "UNC24886001")

# min - H - maj
## H is closer to min
effect.plot("G_Odoribacter", "genus", "RNA", "UNCHS001411")
effect.plot("SV52", "otu", "RNA", "JAX00279640")
effect.plot("G_Odoribacter", "genus", "RNA", "UNCHS001330")
effect.plot("SV21", "otu", "DNA", "JAX00277190")

# H is closer to maj
effect.plot("SV7", "otu", "DNA", "UNCHS014192")
effect.plot("SV35", "otu", "DNA", "JAX00193278")
effect.plot("SV64", "otu", "RNA", "UNC28015845")

#### effect plots for top genes ####
#Adcyap1
effect.plot("SV97", "otu", "RNA", "UNCHS045289")

#gng12
effect.plot("G_Roseburia", "genus", "DNA", "UNCHS017593")
effect.plot("SV40", "otu", "RNA", "UNCHS017593")
effect.plot("SV63", "otu", "RNA", "UNCHS017593")
effect.plot("SV151", "otu", "RNA", "UNCHS017593")

#Cbr3
effect.plot("SV3", "otu", "DNA", "UNCHS043371")
effect.plot("SV3", "otu", "DNA", "UNC27424333")

#Serpina3b
effect.plot("SV158", "otu", "RNA", "UNC21788960")
effect.plot("SV370", "otu", "RNA", "UNC21788960")
effect.plot("G_Marvinbryantia", "genus", "RNA", "UNC21788960")

#Oxgr1
effect.plot("SV22", "otu", "RNA", "UNC24886982")
effect.plot("SV36", "otu", "RNA", "UNC24886982")
effect.plot("SV97", "otu", "RNA", "UNC24886982")

#Prok2
effect.plot("SV35", "otu", "DNA", "UNC11691749")
effect.plot("SV20", "otu", "RNA", "UNC11691749")
effect.plot("SV31", "otu", "RNA", "UNC11691749")
effect.plot("SV35", "otu", "RNA", "UNC11691749")
effect.plot("SV181", "otu", "RNA", "UNC11691749")
effect.plot("G_Lactobacillus", "genus", "RNA", "UNC11691749")
effect.plot("F_Lactobacillaceae", "family", "RNA", "UNC11691749")
effect.plot("O_Lactobacillales", "order", "RNA", "UNC11691749")
effect.plot("C_Bacilli", "class", "RNA", "UNC11691749")

#Tacr3
effect.plot("SV45", "otu", "DNA", "JAX00113489")
effect.plot("SV45", "otu", "DNA", "UNCHS010166")

effect.plot("SV4", "otu", "RNA", "JAX00113489")
effect.plot("SV4", "otu", "RNA", "UNCHS010166")

effect.plot("SV82", "otu", "RNA", "JAX00113489")
effect.plot("SV82", "otu", "RNA", "UNCHS010166")

effect.plot("SV101", "otu", "RNA", "JAX00113489")
effect.plot("SV101", "otu", "RNA", "UNCHS010166")

effect.plot("SV151", "otu", "RNA", "JAX00113489")
effect.plot("SV151", "otu", "RNA", "UNCHS010166")

# Chordc1
effect.plot("SV26", "otu", "RNA", "UNC15924078")
effect.plot("SV281", "otu", "RNA", "UNC15924078")
effect.plot("SV26", "otu", "RNA", "JAX00168681")
effect.plot("SV281", "otu", "RNA", "JAX00168681")
effect.plot("SV26", "otu", "RNA", "UNC15929728")
effect.plot("SV281", "otu", "RNA", "UNC15929728")

#Mchr1
effect.plot("SV36", "otu", "RNA", "UNC25970881")
effect.plot("SV97", "otu", "RNA", "UNC25970881")
effect.plot("SV155", "otu", "RNA", "UNC25970881")
effect.plot("SV170", "otu", "RNA", "UNC25970881")
effect.plot("SV180", "otu", "RNA", "UNC25970881")

#Cyp24a1, Bcas1, Tshz2, Pfdn4, Dok5
effect.plot("G_Eisenbergiella", "genus", "DNA", "UNC4492815")
effect.plot("G_Eisenbergiella", "genus", "DNA", "UNCHS007808")
effect.plot("SV23", "otu", "DNA", "UNC4365491")
effect.plot("SV32", "otu", "RNA", "UNC4600348")

# Tshz2
effect.plot("SV23", "otu", "DNA", "UNC4468510")
effect.plot("SV32", "otu", "RNA", "UNC4473045")

# Adcy2
effect.plot("SV36", "otu", "DNA", "UNCHS036333")
# Chrm3
effect.plot("SV293", "otu", "RNA", "UNC22111265")
effect.plot("SV133", "otu", "RNA", "UNC22119063")
effect.plot("G_Fusicatenibacter", "genus", "RNA", "UNC22119063")

#Rplp2
effect.plot("SV70", "otu", "DNA", "UNC14020115")
effect.plot("G_Odoribacter", "genus", "DNA", "UNC14020115")
effect.plot("SV52", "otu", "RNA", "UNC14020115")
effect.plot("G_Odoribacter", "genus", "RNA", "UNC14020115")
effect.plot("SV3", "otu", "RNA", "JAX00659104")

# GRik3
effect.plot("SV97", "otu", "RNA", "UNCHS012711")
effect.plot("SV98", "otu", "RNA", "UNCHS012711")
effect.plot("unclassified_O_Clostridiales", "genus", "RNA", "UNCHS012711")
effect.plot("SV234", "otu", "RNA", "UNCHS012723")

# Myom2
effect.plot("SV64", "otu", "DNA", "JAX00661779")
effect.plot("SV98", "otu", "DNA", "JAX00661779")
effect.plot("SV119", "otu", "DNA", "JAX00661779")
effect.plot("SV64", "otu", "RNA", "JAX00661779")
effect.plot("SV98", "otu", "RNA", "JAX00661779")
effect.plot("SV396", "otu", "RNA", "JAX00661779")

# Nhlrc3
effect.plot("SV85", "otu", "RNA", "UNC5277799")

# Nebl
effect.plot("SV50", "otu", "DNA", "JAX00484113")
effect.plot("SV113", "otu", "DNA", "JAX00484113")
effect.plot("SV49", "otu", "RNA", "JAX00484113")
effect.plot("SV60", "otu", "RNA", "JAX00484113")

# Slc26a3
effect.plot("SV29", "otu", "DNA", "UNC20836314")
effect.plot("unclassified_C_Deltaproteobacteria", "genus", "DNA", "UNC20836314")
effect.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNC20836314")
effect.plot("G_Acetatifactor", "genus", "RNA", "UNC20836314")

# Ptgfr
effect.plot("SV23", "otu", "DNA", "UNC6545950")
effect.plot("SV23", "otu", "RNA", "UNC6545950")

effect.plot("SV23", "otu", "DNA", "UNC6548766")
effect.plot("SV23", "otu", "RNA", "UNC6548766")

effect.plot("SV99", "otu", "RNA", "UNC6545950")
effect.plot("SV99", "otu", "RNA", "UNC6548766")

effect.plot("SV101", "otu", "RNA", "UNC6545950")
effect.plot("SV101", "otu", "RNA", "UNC6548766")

effect.plot("G_Odoribacter", "genus", "RNA", "UNC6545950")
effect.plot("G_Odoribacter", "genus", "RNA", "UNC6548766")

# Aacs
effect.plot("SV6", "otu", "DNA", "UNCHS015800")
effect.plot("SV10", "otu", "DNA", "UNCHS015800")
effect.plot("SV13", "otu", "DNA", "UNCHS015800")

effect.plot("G_Mucispirillum", "genus", "DNA", "UNCHS015800")
effect.plot("F_Deferribacteraceae", "family", "DNA", "UNCHS015800")
effect.plot("O_Deferribacterales", "order", "DNA", "UNCHS015800")
effect.plot("C_Deferribacteres", "class", "DNA", "UNCHS015800")
effect.plot("P_Deferribacteres", "phylum", "DNA", "UNCHS015800")

# Fxr1
effect.plot("SV22", "otu", "DNA", "JAX00519091")
effect.plot("SV97", "otu", "RNA", "JAX00519091")

# Alpi & Spp2
effect.plot("SV17", "otu", "DNA", "JAX00192085")

effect.plot("SV26", "otu", "RNA", "UNC1098588")

effect.plot("SV370", "otu", "RNA", "UNC1098588")
effect.plot("SV370", "otu", "RNA", "JAX00192085")

effect.plot("G_Marvinbryantia", "genus", "RNA", "UNC1098588")
effect.plot("G_Marvinbryantia", "genus", "RNA", "JAX00192085")

# Tomm22
effect.plot("SV155", "otu", "RNA", "UNCHS040615")
effect.plot("SV170", "otu", "RNA", "UNCHS040615")
effect.plot("SV180", "otu", "RNA", "UNCHS040615")

effect.plot("SV36", "otu", "RNA", "UNC25970881")
effect.plot("SV36", "otu", "RNA", "UNC25960296")

# Grm8
effect.plot("SV113", "otu", "DNA", "UNCHS016696")
effect.plot("SV113", "otu", "DNA", "UNC10780304")

# Grm8
effect.plot("SV6", "otu", "DNA", "UNCHS019799")
effect.plot("SV13", "otu", "DNA", "UNCHS019799")
effect.plot("SV25", "otu", "RNA", "UNCHS019873")
effect.plot("SV49", "otu", "RNA", "JAX00150794")
effect.plot("G_Odoribacter", "genus", "RNA", "UNC12506197")







