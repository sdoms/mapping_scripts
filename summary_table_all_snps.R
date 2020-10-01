setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/")
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/")
library(ggplot2)
library(patchwork)
musdom <- read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv", sep=",")
load("../../../Cleaning_snps/consensus_F0_all.Rdata")
parents <- as.data.frame(reduced_geno)
parents$X <- rownames(parents)
musdom <- merge(musdom, parents, by="X", )
#trait <- "class"
for (trait in c("phylum", "class", "order", "family", "genus", "otu")){
  tax_table <- read.csv(paste0("../../../Phenotyping_27.02.2020/tables_core/",trait,"_tax_table_f2_core.csv" ))
  
  
  #gws <- "C_Bacilli"
  Xtrue <- TRUE
  chrom_range <- c(1:19)
  sig.P <- data.frame()
  sig.add.P <- data.frame()
  sig.dom.P <- data.frame()
  colnames(tax_table)<- tolower(colnames(tax_table))
  if (trait =="otu"){
    colnames(tax_table) <- c("otu", "kingdom", "phylum", "class", "order", "family", "genus")
  }
  for (gws in tax_table[,trait]){
    gwscan <- data.frame()
    for (chr in 1:19){
      tx <- readRDS(paste0(trait, "/",gws, "_chr_", chr, "with_add_dom.rds"))
      gwscan <- rbind(gwscan, tx)
    }
    if (file.exists(paste0(trait,"/",gws, "_chrX.rds"))){
      X_chrom <- readRDS(paste0(trait,"/",gws, "_chrX.rds"))
      #X_chrom$chr <- "20"
      gwscan <- dplyr::bind_rows(gwscan, X_chrom)
      chrom_range <- c(1:19, "X")
      
    }
    # Bonferonni threshold 
    threshold <- 0.05/nrow(gwscan)
    
    sig_gwscan_P <- gwscan[which(gwscan$P<threshold),]
    sig_gwscan_add.P <- gwscan[which(gwscan$add.P<threshold),]
    sig_gwscan_dom.P <- gwscan[which(gwscan$dom.P<threshold),]
    sig.P <- rbind(sig.P, sig_gwscan_P)
    sig.add.P <- rbind(sig.add.P, sig_gwscan_add.P)
    sig.dom.P <- rbind(sig.dom.P, sig_gwscan_dom.P)
    
  }
  # add mus dom
  sig.P <- merge(sig.P, musdom, by.x = "marker", by.y="X", all.x = T)
  sig.add.P <- merge(sig.add.P, musdom, by.x = "marker", by.y="X", all.x = T)
  sig.dom.P <- merge(sig.dom.P, musdom, by.x = "marker", by.y="X", all.x = T)
  
  write.table(x = unique(sig.P), file=paste0(trait, "_all_sig_snps_P.txt"))
  write.table(x = unique(sig.add.P), file=paste0(trait, "_all_sig_snps_add_P.txt"))
  write.table(x = unique(sig.dom.P), file=paste0(trait, "_all_sig_snps_dom_P.txt"))
  
  all_sig1 <- rbind(sig.P, sig.add.P)
  all_sig <- rbind(all_sig1, sig.dom.P)
  deduped.data <- unique( all_sig )
  if (trait == "family" | trait == "genus"){
    deduped.data <- deduped.data[which(deduped.data$marker!="UNCHS040366"),]
  }

  add <- ggplot(deduped.data, aes(x=add.T)) + geom_histogram() + labs(x="Additive effect Z-score ")+ theme(axis.title = element_text(size=9))
  dom <- ggplot(deduped.data, aes(x=dom.T)) +geom_histogram(colour="black", fill="white")+ labs(x="Dominance effect Z-score")+ theme(axis.title = element_text(size=9))  
  
  
  add+dom + plot_annotation(title=trait, tag_levels = "A")
  ggsave(paste0(trait, "_all_sig_snps-zscore.pdf"))
  
} 

###### all snps taxa ###

#musdom <- read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv", sep=";")
#trait <- "class"
sig.P <- data.frame()
sig.add.P <- data.frame()
sig.dom.P <- data.frame()
for (trait in c("phylum", "class", "order", "family", "genus", "otu")){
  tax_table <- read.csv(paste0("../../../Phenotyping_27.02.2020/tables_core/",trait,"_tax_table_f2_core.csv" ))
  #gws <- "C_Bacilli"
  Xtrue <- TRUE

  chrom_range <- c(1:19)

  colnames(tax_table)<- tolower(colnames(tax_table))
  if (trait =="otu"){
    colnames(tax_table) <- c("otu", "kingdom", "phylum", "class", "order", "family", "genus")
  }
  
  for (gws in tax_table[,trait]){
    gwscan <- data.frame()
    for (chr in 1:19){
      tx <- readRDS(paste0(trait, "/",gws, "_chr_", chr, "with_add_dom.rds"))
      gwscan <- rbind(gwscan, tx)
    }
    if (file.exists(paste0(trait,"/",gws, "_chrX.rds"))){
      X_chrom <- readRDS(paste0(trait,"/",gws, "_chrX.rds"))
      #X_chrom$chr <- "20"
      gwscan <- dplyr::bind_rows(gwscan, X_chrom)
      chrom_range <- c(1:19, "X")
      
    }
    # Bonferonni threshold 
    threshold <- 0.05/nrow(gwscan)
    
    sig_gwscan_P <- gwscan[which(gwscan$P<threshold),]
    sig_gwscan_add.P <- gwscan[which(gwscan$add.P<threshold),]
    sig_gwscan_dom.P <- gwscan[which(gwscan$dom.P<threshold),]
    sig.P <- rbind(sig.P, sig_gwscan_P)
    sig.add.P <- rbind(sig.add.P, sig_gwscan_add.P)
    sig.dom.P <- rbind(sig.dom.P, sig_gwscan_dom.P)
    
  }
 
} 

# add mus dom
sig.P <- merge(sig.P, musdom, by.x = "marker", by.y="X", all.x = T)
sig.add.P <- merge(sig.add.P, musdom, by.x = "marker", by.y="X", all.x = T)
sig.dom.P <- merge(sig.dom.P, musdom, by.x = "marker", by.y="X", all.x = T)

write.table(x = unique(sig.P), file="ALL_sig_snps_P.txt")
write.table(x = unique(sig.add.P), file= "ALL_sig_snps_add_P.txt")
write.table(x = unique(sig.dom.P), file="ALL_sig_snps_dom_P.txt")

all_sig1 <- rbind(sig.P, sig.add.P)
all_sig <- rbind(all_sig1, sig.dom.P)
deduped.data <- unique( all_sig )

  deduped.data <- deduped.data[which(deduped.data$marker!="UNCHS040366"),]

add <- ggplot(deduped.data, aes(x=add.T)) + geom_histogram() + labs(x="Additive effect Z-score ")+ theme(axis.title = element_text(size=9))
dom <- ggplot(deduped.data, aes(x=dom.T)) +geom_histogram(colour="black", fill="white")+ labs(x="Dominance effect Z-score")+ theme(axis.title = element_text(size=9))  


add+dom + plot_annotation(title="All RNA", tag_levels = "A")
ggsave("ALL_sig_snps-zscore.pdf")

write.table(x = deduped.data, file="ALL_sig_snps_allP.txt")

library(dplyr)
ex <- deduped.data %>%
  group_by(marker) %>%
  mutate(all_taxa = paste(tax, collapse = " | "))  %>% 
  add_count(marker, name="count") %>% 
  mutate(chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>%
  select(marker,chr, pos,A=A2,B=A1,mus,dom, AA,AB,BB,count, all_taxa, FSP3, FSP5, HAP1, HAP2, HOP3, HOP6, TUP3,TUP4) %>%  
  distinct() %>% 
  arrange(chr,pos)
write.table(ex,"summary_table_DNA_markers_count.txt")  

write.table(ex,"summary_table_RNA_markers_count.txt")  




# pretty table ##### for paper? 
devtools::install_github("glin/reactable")
library(reactable)
reactable(ex)

devtools::install_github("haozhu233/kableExtra")
library(kableExtra)
ex %>%
  kbl() %>%
  kable_classic(full_width = T)
