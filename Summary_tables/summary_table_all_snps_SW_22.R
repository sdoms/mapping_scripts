setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/")

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


library(stringr)
library(patchwork)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db) 
library(tidyverse)

threshold <- 0.05/32625/118.2177

musdom <- read.csv("~/Documents/Research/Experiments/Final_QTL_mapping/Results/allele_freq_consensus_mus_dom.csv", sep=",")
musdom$dd_frew <- NULL
load("../../Cleaning_snps/consensus_F0_all.Rdata")
parents <- as.data.frame(reduced_geno)
parents$X <- rownames(parents)
musdom <- merge(musdom, parents, by="X", )
DNAorRNA <- "DNA"
outputdir <- "Genes/all_genes_snps/"

#####################################################################
##                        Summary table 1:                         ##
##  a table that lists all significant SNPs for a taxonomic level  ##
#####################################################################

# genome-wide 
DNAorRNA<- "DNA"
for (DNAorRNA in c("DNA", "RNA")){
  lambdaGC <- read.table(paste0("Revisions/ld_score_regression/lambdaGC_", DNAorRNA, ".csv"), sep=";")
  lambdaGC$taxa <- rownames(lambdaGC)
  
  taxa_for_correction <- filter(lambdaGC, lambdaGC_P>1.05) 
  all_taxa <- read.table(paste0("../../Scripts/Data/trait.list.", tolower(DNAorRNA), ".csv"), sep=";", header=T)
  tot_files <-tibble()
  # trait <- "otu"
  for (trait in c("otu", "genus", "family", "order", "class", "phylum")){
    all_files <- tibble()
    for (tax in all_taxa[,trait]){
      if (tax %in% taxa_for_correction$taxa){
        if (file.exists(paste0("./",DNAorRNA, "/sig.summaries_pval_corr/", trait, "/", tax,"_gwscan_sigSNPS.Rdata"))){
          infile<-loadRData(paste0("./",DNAorRNA, "/sig.summaries_pval_corr/", trait, "/", tax,"_gwscan_sigSNPS.Rdata"))
          if (nrow(infile)==0) next
          infile$tax <- tax
          infile$Pval_corr <- "Y"
          all_files <- rbind(all_files, infile)
        }
        
        
      } else {
        if (file.exists(paste0("./",DNAorRNA, "/sig.summaries/", trait, "/", tax,"_gwscan_sigSNPS.Rdata"))){
          infile<-loadRData(paste0("./",DNAorRNA, "/sig.summaries/", trait, "/", tax,"_gwscan_sigSNPS.Rdata"))
          if (nrow(infile)==0) next
          infile$tax <- tax
          infile$Pval_corr <- "N"
          all_files <- rbind(all_files, infile)
        }
        
        
        
      }
    }
    tot_files <- rbind(tot_files, all_files)
    write_delim(all_files, file=paste0(outputdir,  DNAorRNA, "_", trait, "_all_sig_snps_pval_corr.csv"), delim=";")
    save(all_files, file=paste0(outputdir, DNAorRNA, "_", trait, "_all_sig_snps_pval_corr.Rdata"))
  }
  
  write_delim(tot_files, file=paste0(outputdir,  DNAorRNA, "_all_sig_snps_pval_corr.csv"), delim=";")
  save(tot_files, file=paste0(outputdir, DNAorRNA, "_all_sig_snps_pval_corr.Rdata"))
  
}
# study-wide 
for (DNAorRNA in c("DNA", "RNA")){
load(paste0(outputdir, DNAorRNA, "_all_sig_snps_pval_corr.Rdata"))
sig.P <- tot_files %>% filter(P<threshold)
sig.add.p  <- tot_files %>% filter(add.P<threshold)
sig.dom.p  <- tot_files %>% filter(dom.P<threshold)

sig.all <- tot_files %>% filter(P<threshold|add.P<threshold|dom.P<threshold)

write_delim(sig.P, file=paste0(outputdir,  DNAorRNA, "_all_sig_P_snps_pval_corr_SW.csv"), delim=";")
save(sig.P, file=paste0(outputdir, DNAorRNA, "_all_sig_P_snps_pval_corr_SW.Rdata"))
write_delim(sig.add.p, file=paste0(outputdir,  DNAorRNA, "_all_sig_addP_snps_pval_corr_SW.csv"), delim=";")
save(sig.add.p, file=paste0(outputdir, DNAorRNA, "_all_sig_addP_snps_pval_corr_SW.Rdata"))
write_delim(sig.dom.p, file=paste0(outputdir,  DNAorRNA, "_all_sig_domP_snps_pval_corr_SW.csv"), delim=";")
save(sig.dom.p, file=paste0(outputdir, DNAorRNA, "_all_sig_domP_snps_pval_corr_SW.Rdata"))

write_delim(sig.all, file=paste0(outputdir,  DNAorRNA, "_all_sig_snps_pval_corr_SW.csv"), delim=";")
save(sig.all, file=paste0(outputdir, DNAorRNA, "_all_sig_snps_pval_corr_SW.Rdata"))
}
################################################################################################
##                                      Summary table 2:                                      ##
##  a table that lists all uniaue significant SNPs for all P values for all taxonomic levels  ##
##              with a count of the times it was significant in different traits              ##
################################################################################################

all_sig_snps_dna <- loadRData(paste0(outputdir,"DNA_all_sig_snps_pval_corr_SW.Rdata"))
all_sig_snps_dna$DNAorRNA <- "DNA"
all_sig_snps_rna <- loadRData(paste0(outputdir,"RNA_all_sig_snps_pval_corr_SW.Rdata"))
all_sig_snps_rna$DNAorRNA <- "RNA"

all_sig_snps <- rbind(all_sig_snps_dna, all_sig_snps_rna)

ex_dna <- all_sig_snps_dna %>% left_join(musdom, by=c("marker"="X")) %>%
  group_by(marker) %>%
  mutate(all_taxa = paste(tax, collapse = " | "))  %>% 
  add_count(marker, name="count") %>% ungroup() %>% 
  mutate(chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>%
  dplyr::select(marker,chr, pos,A=A2,B=A1,mus,dom, AA,AB,BB,count, all_taxa, FSP3, FSP5, HAP1, HAP2, HOP3, HOP6, TUP3,TUP4) %>%  
  distinct() %>% 
  arrange(chr,pos) %>%
  mutate(chr = replace_na(as.character(chr), "X"))
write.table(ex_dna,paste0(outputdir,"summary_table_DNA_markers_count_pval_cor_SW.txt")  )

ex_rna <- all_sig_snps_rna %>% left_join(musdom, by=c("marker"="X")) %>%
  group_by(marker) %>%
  mutate(all_taxa = paste(tax, collapse = " | "))  %>% 
  add_count(marker, name="count") %>% ungroup() %>% 
  mutate(chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>%
  dplyr::select(marker,chr, pos,A=A2,B=A1,mus,dom, AA,AB,BB,count, all_taxa, FSP3, FSP5, HAP1, HAP2, HOP3, HOP6, TUP3,TUP4) %>%  
  distinct() %>% 
  arrange(chr,pos) %>%
  mutate(chr = replace_na(as.character(chr), "X"))
write.table(ex_rna,paste0(outputdir,"summary_table_RNA_markers_count_pval_cor_SW.txt")  )

ex_all <- all_sig_snps %>% left_join(musdom, by=c("marker"="X")) %>%
  group_by(marker) %>%
  mutate(all_taxa = paste(tax, collapse = " | "))  %>% 
  add_count(marker, name="count") %>% ungroup() %>% 
  mutate(chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>%
  dplyr::select(marker,chr, pos,A=A2,B=A1,mus,dom, AA,AB,BB,count, all_taxa, FSP3, FSP5, HAP1, HAP2, HOP3, HOP6, TUP3,TUP4) %>%  
  distinct() %>% 
  arrange(chr,pos) %>%
  mutate(chr = replace_na(as.character(chr), "X"))
write.table(ex_rna,paste0(outputdir,"summary_table_all_markers_count_pval_cor_SW.txt")  )

################################################################################################
##                                      Summary table 3:                                      ##
##  a table that lists all unique significant SNPs for all P values for all taxonomic levels  ##
##              with a count of the times it was significant in different traits              ##
##                            with the closest genes to the marker                            ##
################################################################################################

# Add closest genes to al snps
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db) 


ex <- ex_rna
ex <- ex_dna
input <- ex %>%
  drop_na(chr) %>% 
  mutate(chr=paste0("chr", chr))%>% 
  dplyr::select(rsid=marker, chr, pos) 
#input <- input[1:3,]
input
#         rsid  chr       pos
# 1  rs3753344 chr1   1142150
# 2 rs12191877 chr6  31252925
# 3   rs881375 chr9 123652898

target <- with(input,
               GRanges( seqnames = Rle(chr),
                        ranges   = IRanges(pos, end=pos, names=rsid),
                        strand   = Rle(strand("*")) ) )
target







loc <- locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants())
loc
names(loc) <- NULL

p_ids <- unlist(loc$PRECEDEID, use.names=FALSE) 
exp_ranges <- rep(loc,  elementNROWS(loc$PRECEDEID))
p_dist <- GenomicRanges::distance(exp_ranges, TxDb.Mmusculus.UCSC.mm10.knownGene, id=p_ids, type="gene")
head(p_dist)
exp_ranges$PRECEDE_DIST <- p_dist
exp_ranges

## Collapsed view of ranges, gene id and distance:
loc$PRECEDE_DIST <- relist(p_dist, loc$PRECEDEID)
loc


f_ids <- unlist(loc$FOLLOWID, use.names=FALSE) 
exp_ranges_f <- rep(loc,  elementNROWS(loc$FOLLOWID))

f_dist <- GenomicRanges::distance(exp_ranges_f, TxDb.Mmusculus.UCSC.mm10.knownGene, id=f_ids, type="gene")
head(f_dist)
exp_ranges_f$FOLLOW_DIST <- f_dist
exp_ranges_f

## Collapsed view of ranges, gene id and distance:
loc$FOLLOW_DIST <- relist(f_dist, loc$FOLLOWID)
loc

# find closest
# find closest
loc$PRECEDEID_place <-lapply(loc$PRECEDE_DIST, function(x) which.min(x))
loc$PRECEDEID_place <- as.numeric(as.character(loc$PRECEDEID_place))

#loc$PRECEDEID_closest <- lapply(loc$PRECEDEID, function(x) x[unlist(loc$PRECEDEID_place)])

loc$FOLLOWID_place <-lapply(loc$FOLLOW_DIST, function(x) which.min(x))
loc$FOLLOWID_place <- as.numeric(as.character(loc$FOLLOWID_place))

#out$FOLLOWID_closest <- apply(out, 1, function(x) dplyr::nth(x[14],unlist(x[20]), default=NULL))


out <- as.data.frame(loc)
out$names <- names(target)[ out$QUERYID ]

precede_ID <- data.frame(out$PRECEDEID)
colnames(precede_ID)<- "ID"
precede_ID$place <- out$PRECEDEID_place
for (x in 1:nrow(precede_ID)){
  n <- precede_ID[x,"place"]
  precede_ID[x,"closest"] <- as.vector(precede_ID[x,1])[[1]][n]
}

follow_ID <- data.frame(out$FOLLOWID)
colnames(follow_ID)<- "ID"
follow_ID$place <- out$FOLLOWID_place
for (x in 1:nrow(follow_ID)){
  n <- follow_ID[x,"place"]
  follow_ID[x,"closest"] <- as.vector(follow_ID[x,1])[[1]][n]
}


out <- out[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID")]
out <- cbind(out, precede_ID$closest)
out <- cbind(out, follow_ID$closest)

out <- unique(out)
out

colnames(out)<- c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "PRECEDEID", "FOLLOWID")
Symbol2id <- as.list( org.Mm.egSYMBOL2EG )
id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
names(id2Symbol) <- unlist(Symbol2id)

x <- unique( with(out, c(levels(as.factor(GENEID)), levels(as.factor(PRECEDEID)), levels(as.factor(FOLLOWID)))))
table( x %in% names(id2Symbol)) # good, all found

out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
out$PRECEDESYMBOL <- id2Symbol[ as.character(out$PRECEDEID) ]
out$FOLLOWSYMBOL <- id2Symbol[ as.character(out$FOLLOWID) ]
out
final_out <-merge(ex, out, by.x="marker", by.y="names", all.x=T)

ex_dist <- final_out %>%
  group_by(marker) %>%
  mutate(locations= paste(LOCATION, collapse = " | "),all_genes = paste(GENESYMBOL, collapse = " | "), all_preceding_genes = paste(PRECEDESYMBOL, collapse = " | "), all_following_genes = paste(FOLLOWSYMBOL, collapse = " | ")) %>% 
  distinct(marker, locations,  all_genes, all_preceding_genes, all_following_genes) %>% 
  left_join(final_out[,c(1:20)], by="marker") %>% 
  distinct() %>% 
  mutate(chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>%
  arrange(chr,pos) %>% 
  mutate(chr=replace_na(as.character(chr),"X"))

write.csv(ex_dist, paste0(outputdir,"markers_with_genes_DNA.csv"), row.names = F)

write.csv(ex_dist, paste0(outputdir,"markers_with_genes_RNA.csv"), row.names = F)


final_out_dna <- read.csv(paste0(outputdir, "markers_with_genes_DNA.csv"))
final_out_rna <- read.csv(paste0(outputdir, "markers_with_genes_RNA.csv"))

final_out_dna$DNA_RNA <- "DNA"
final_out_rna$DNA_RNA <- "RNA"

all_final <- rbind(final_out_dna, final_out_rna)

tot_count <-all_final %>%
  group_by(marker, DNA_RNA) %>% 
  mutate(taxa=ifelse(DNA_RNA == "DNA", paste("DNA:", all_taxa), paste("RNA:", all_taxa))) %>% 
  group_by(marker) %>% 
  summarise(total_count = sum(count))
ex_all <- all_final %>%
  group_by(marker, DNA_RNA) %>% 
  mutate(taxa=ifelse(DNA_RNA == "DNA", paste("DNA:", all_taxa), paste("RNA:", all_taxa))) %>% 
  group_by(marker) %>% 
  mutate(all_taxa = paste(taxa, collapse = " | ")) %>% 
  dplyr::select(-c(taxa, DNA_RNA, count)) %>% 
  distinct(marker, .keep_all = T) %>% 
  ungroup() %>% 
  left_join(tot_count)
write.csv(ex_all, file=paste0(outputdir,"all_markers_RNA_DNA-with-genes_SW_pval_corr.csv"))

################################################################################################
##                                      Summary table 4:                                      ##
##  a table that calculates the amount of significant regions per taxon                       ##
################################################################################################
library(ggsci)
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_22022022_SW_pval_corr.Rdata")
tax_table<- read.table("~/Documents/Research/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/otu_tax_table_f2_core.csv", sep=";", header = T)[,c(3,10)]
count_taxa <- all_dna_rna_SW %>% group_by(trait, P.type, dna.rna) %>% summarise(n=n()) %>% ungroup() %>% left_join(tax_table, by=c("trait"="SV")) %>% 
  mutate(taxon=ifelse(!is.na(Genus), paste0("A",trait, " (", Genus, ")"), trait), n_min= ifelse(dna.rna=="DNA", n*-1, n)) 

ggplot(count_taxa, aes(x=reorder(taxon, n, sum), y=n_min, fill=P.type))+geom_bar(stat = "identity", position="stack")+coord_flip()+theme_minimal()+scale_fill_d3()+
  labs(x="", y="Number of siginificant loci", fill="P value")+geom_hline(yintercept = 0)
ggsave("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/sig_associations_per_taxa_dna_rna.pdf")

ggplot(count_taxa, aes(x=reorder(taxon, n, sum), y=n, fill=P.type))+geom_bar(stat = "identity", position="stack")+coord_flip()+theme_test()+scale_fill_d3()+
  labs(x="", y="Number of siginificant loci", fill="P value")
ggsave("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/sig_associations_per_taxa.pdf")

