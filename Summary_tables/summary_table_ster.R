
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacula/")
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(qvalue)
musdom <- read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv", sep=",")
load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/consensus_F0_all.Rdata")
parents <- as.data.frame(reduced_geno)
parents$X <- rownames(parents)
load("~/Documents/PhD/Experiments/Final_QTL_mapping/Genotype_data/GM_snps.Rdata")
parents2 <- merge(GM_snps[,c("marker","rsID")],parents, by.x= "marker", by.y="X")
musdom <- merge(musdom, parents2,by.x="X", by.y="marker", )


#trait <- "class"
sig.P <- data.frame()
sig.add.P <- data.frame()
sig.dom.P <- data.frame()
for (trait in c("bac_PC01","bac_PC02","bac_PC03","bac_PC04","bac_PC05","bac_PC06",
                "bac_PC07","bac_PC08","bac_PC09","bac_PC10","bac_PC11","bac_PC12" )){

  #tax_table <- read.csv(paste0("../../../Phenotyping_27.02.2020/tables_core/",trait,"_tax_table_f2_core.csv" ))

  #gws <- "C_Bacilli"
  Xtrue <- TRUE
  chrom_range <- c(1:19)

  
    gwscan <- data.frame()
    for (chr in 1:19){
      tx <- readRDS(paste0(trait, "/",trait, "_chr_", chr, "with_add_dom.rds"))
      gwscan <- rbind(gwscan, tx)
    }
    if (file.exists(paste0(trait,"/",trait, "_chrX.rds"))){
      X_chrom <- readRDS(paste0(trait,"/",trait, "_chrX.rds"))
      #X_chrom$chr <- "20"
      gwscan <- dplyr::bind_rows(gwscan, X_chrom)
      chrom_range <- c(1:19, "X")
      
    }
    # Bonferonni threshold 
    threshold_P <- 0.05/nrow(gwscan)
    threshold_addP <- threshold_P
    threshold_domP <- threshold_P
    
    # q_P <-qvalue(gwscan[,"P"])
    # threshold_P <- (max(na.omit(q_P$pvalues[q_P$qvalues <=0.1])))
    # 
    # q_add <-qvalue(gwscan[,"add.P"])
    # threshold_addP <- (max(na.omit(q_add$pvalues[q_add$qvalues <=0.1])))
    # 
    # q_dom <-qvalue(gwscan[,"dom.P"])
    # threshold_domP <- (max(na.omit(q_dom$pvalues[q_dom$qvalues <=0.1])))
    # 
    sig_gwscan_P <- gwscan[which(gwscan$P<threshold_P),]
    sig_gwscan_add.P <- gwscan[which(gwscan$add.P<threshold_addP),]
    sig_gwscan_dom.P <- gwscan[which(gwscan$dom.P<threshold_domP),]
    
    sig.P <- rbind(sig.P, sig_gwscan_P)
    sig.add.P <- rbind(sig.add.P, sig_gwscan_add.P)
    sig.dom.P <- rbind(sig.dom.P, sig_gwscan_dom.P)
    
}
  # add mus dom
  sig.P <- merge(sig.P, musdom, by.x = "marker", by.y="X", all.x = T)
  sig.add.P <- merge(sig.add.P, musdom, by.x = "marker", by.y="X", all.x = T)
  sig.dom.P <- merge(sig.dom.P, musdom, by.x = "marker", by.y="X", all.x = T)
  
  write.table(x = unique(sig.P), file="PC_bacula_all_sig_snps_P.txt")
  write.table(x = unique(sig.add.P), file="PC_bacula_all_sig_snps_add_P.txt")
  write.table(x = unique(sig.dom.P), file= "PC_bacula_all_sig_snps_dom_P.txt")
  
  all_sig1 <- rbind(sig.P, sig.add.P)
  all_sig <- rbind(all_sig1, sig.dom.P)
  deduped.data <- unique( all_sig )
 
  
  add <- ggplot(deduped.data, aes(x=add.T)) + geom_histogram() + labs(x="Additive effect Z-score ")+ theme(axis.title = element_text(size=9))
  dom <- ggplot(deduped.data, aes(x=dom.T)) +geom_histogram(colour="black", fill="white")+ labs(x="Dominance effect Z-score")+ theme(axis.title = element_text(size=9))  
  
  
  add+dom + plot_annotation(title=trait, tag_levels = "A")
  ggsave("bacula_PC_all_sig_snps-zscore.pdf")
  

tt <- deduped.data %>% 
  mutate(D.A=dom.T/add.T, MAF=(AA/n+AA/n+AB/n)/2,chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>% 
  arrange(chr, pos) %>% 
  mutate(chr=replace_na(chr, "X"))
write.table(x = tt, file="bacula_PC_ALL_sig_snps_allP.txt")


ex <- tt %>%
  group_by(marker) %>%
  mutate(all_taxa = paste(tax, collapse = " | "))  %>% 
  add_count(marker, name="count") %>% 
  mutate(chr = as.numeric(as.character(chr)), pos=as.numeric(as.character(pos))) %>%
  dplyr::select(marker,chr, pos,A=A2,B=A1,mus,dom, AA,AB,BB,count, all_taxa, FSP3, FSP5, HAP1, HAP2, HOP3, HOP6, TUP3,TUP4) %>%  
  distinct() %>% 
  arrange(chr,pos) %>% 
  mutate(chr=replace_na(chr, "X"))
write.table(ex,"bacula_PC_summary_table_markers_count.txt") 

 

# Add closest genes to al snps
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db) 
library(tidyr)





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
write.csv(final_out, "bacula_PC_markers_with_genes.csv", row.names = F)





# pretty table ##### for paper? 
devtools::install_github("glin/reactable")
library(reactable)
reactable(ex)

devtools::install_github("haozhu233/kableExtra")
library(kableExtra)
ex %>%
  kbl() %>%
  kable_classic(full_width = T)

