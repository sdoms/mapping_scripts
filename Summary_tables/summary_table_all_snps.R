setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/")
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/")
library(ggplot2)
library(patchwork)
library(dplyr)
musdom <- read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv", sep=",")
load("../../../Cleaning_snps/consensus_F0_all.Rdata")
parents <- as.data.frame(reduced_geno)
parents$X <- rownames(parents)
load("../../../Genotype_data/GM_snps.Rdata")
parents2 <- merge(GM_snps[,c("marker","rsID")],parents, by.x= "marker", by.y="X")
musdom <- merge(musdom, parents2,by.x="X", by.y="marker", )

outputdir <- "summary_tables/genome_wide_significant/all_signif_snps/"
#####################################################################
##                        Summary table 1:                         ##
##  a table that lists all significant SNPs for a taxonomix level  ##
#####################################################################


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
    tax_table <- tax_table[,c(3,5:10)]
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
  
  write.table(x = unique(sig.P), file=paste0(outputdir,trait, "_all_sig_snps_P.txt"))
  write.table(x = unique(sig.add.P), file=paste0(outputdir,trait, "_all_sig_snps_add_P.txt"))
  write.table(x = unique(sig.dom.P), file=paste0(outputdir, trait, "_all_sig_snps_dom_P.txt"))
  
  ###########################################################################################
  ##                                   Summary table 2:                                    ##
  ##  a table that lists all uniaue significant SNPs for all P values per taxonomix level  ##
  ###########################################################################################
  
  
  all_sig1 <- rbind(sig.P, sig.add.P)
  all_sig <- rbind(all_sig1, sig.dom.P)
  deduped.data <- unique( all_sig )
  if (trait == "family" | trait == "genus"){
    deduped.data <- deduped.data[which(deduped.data$marker!="UNCHS040366"),]
  }
  
  ##################################################################
  ##           Plot zscores all unique significant snps           ##
  ##################################################################
  
  
  add <- ggplot(deduped.data, aes(x=add.T)) + geom_histogram() + labs(x="Additive effect Z-score ")+ theme(axis.title = element_text(size=9))
  dom <- ggplot(deduped.data, aes(x=dom.T)) +geom_histogram(colour="black", fill="white")+ labs(x="Dominance effect Z-score")+ theme(axis.title = element_text(size=9))  
  
  
  add+dom + plot_annotation(title=trait, tag_levels = "A")
  ggsave(paste0(outputdir, trait, "_all_sig_snps-zscore.pdf"))
  
} 

###### all snps taxa #####

###########################################################################################
##                                   Summary table 3:                                    ##
##  a table that lists all uniaue significant SNPs per P value for all taxonomic levels  ##
###########################################################################################


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
    tax_table <- tax_table[,c(3,5:10)]
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

# write summary table 3
write.table(x = unique(sig.P), file=paste0(outputdir,"ALL_sig_snps_P.txt"))
write.table(x = unique(sig.add.P), file= paste0(outputdir,"ALL_sig_snps_add_P.txt"))
write.table(x = unique(sig.dom.P), file=paste0(outputdir,"ALL_sig_snps_dom_P.txt"))

all_sig1 <- rbind(sig.P, sig.add.P)
all_sig <- rbind(all_sig1, sig.dom.P)
deduped.data <- unique( all_sig )

  deduped.data <- deduped.data[which(deduped.data$marker!="UNCHS040366"),]

add <- ggplot(deduped.data, aes(x=add.T)) + geom_histogram() + labs(x="Additive effect Z-score ")+ theme(axis.title = element_text(size=9))
dom <- ggplot(deduped.data, aes(x=dom.T)) +geom_histogram(colour="black", fill="white")+ labs(x="Dominance effect Z-score")+ theme(axis.title = element_text(size=9))  


add+dom + plot_annotation(title="All", tag_levels = "A")
ggsave(paste0(outputdir,"ALL_sig_snps-zscore.pdf"))

################################################################################################
##                                      Summary table 4:                                      ##
##  a table that lists all uniaue significant SNPs for all P values for all taxonomic levels  ##
################################################################################################


tt <- deduped.data %>% 
  mutate(D.A=dom.T/add.T, MAF=(AA/n+AA/n+AB/n)/2,chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>% 
  arrange(chr, pos)
write.table(x = tt, file=paste0(outputdir,"ALL_sig_snps_allP.txt"))


################################################################################################
##                                      Summary table 5:                                      ##
##  a table that lists all uniaue significant SNPs for all P values for all taxonomic levels  ##
##              with a count of the times it was significant in different traits              ##
################################################################################################


ex <- tt %>%
  group_by(marker) %>%
  mutate(all_taxa = paste(tax, collapse = " | "))  %>% 
  add_count(marker, name="count") %>% 
  mutate(chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>%
  dplyr::select(marker,chr, pos,A=A2,B=A1,mus,dom, AA,AB,BB,count, all_taxa, FSP3, FSP5, HAP1, HAP2, HOP3, HOP6, TUP3,TUP4) %>%  
  distinct() %>% 
  arrange(chr,pos) %>%
  mutate(chr = replace_na(chr, "X"))
write.table(ex,paste0(outputdir,"summary_table_DNA_markers_count.txt")  )

write.table(ex,paste0(outputdir,"summary_table_RNA_markers_count.txt")  )

# Add closest genes to al snps
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db) 
library(tidyr)


################################################################################################
##                                      Summary table 6:                                      ##
##  a table that lists all uniaue significant SNPs for all P values for all taxonomic levels  ##
##              with a count of the times it was significant in different traits              ##
##                            with the closest genes to the marker                            ##
################################################################################################


input <- ex %>%
  drop_na(chr) %>% 
  mutate(chr=paste0("chr", chr), pos=pos*1e6)%>% 
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
  mutate(chr=replace_na(chr,"X"))

write.csv(ex_dist, paste0(outputdir,"markers_with_genes_DNA.csv"), row.names = F)

write.csv(ex_dist, paste0(outputdir,"markers_with_genes_RNA.csv"), row.names = F)

final_out_dna <- read.csv(paste0("../DNA/", outputdir, "markers_with_genes_DNA.csv"))
final_out_rna <- read.csv(paste0("../RNA/",outputdir, "markers_with_genes_RNA.csv"))

final_out_dna$DNA_RNA <- "DNA"
final_out_rna$DNA_RNA <- "RNA"

all_final <- rbind(final_out_dna, final_out_rna)

ex_all <- all_final %>%
  group_by(marker, DNA_RNA) %>% 
  mutate(taxa=ifelse(DNA_RNA == "DNA", paste("DNA:", all_taxa), paste("RNA:", all_taxa))) %>% 
  group_by(marker) %>% 
  mutate(all_taxa = paste(taxa, collapse = " | ")) %>% 
  dplyr::select(-c(taxa, DNA_RNA,AA,AB, BB, count)) %>% 
  distinct(marker, .keep_all = T) 
write.csv(ex_all, file="../Shared/all_markers_RNA_DNA-with-genes_GW.csv")
