setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/")
# load("../all_DNA_RNA_18032021_SW.Rdata")
load("../all_DNA_RNA_22022022_SW_pval_corr.Rdata")
library(tidyverse)
piece <- all_dna_rna_SW %>% 
  mutate(length=stop.LD.pos-start.LD.pos) %>% 
  filter(length<1e6&length>0) %>% 
  arrange(length)

ex_dist <- piece %>%
  group_by(chr, peak.pos) %>%
  mutate(taxa= paste(dna.rna, ":",trait, collapse = " | ")) %>% 
  distinct(chr, peak.pos,  .keep_all=T) %>% arrange(length) %>% 
  dplyr::select(taxa, length, chr, peak.pos,markers, list_protein_coding_genes, list_total_genes, closest_gene, closest_gene_type, closest_gene_description) %>% 
  ungroup()
write_delim(ex_dist, "genes_1Mb_intervals_pval_corr.csv", delim=";" )

library(VariantAnnotation)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
input <- ex_dist %>%
  drop_na(chr) %>% 
  mutate(chr=paste0("chr", chr), pos=peak.pos)%>% 
  dplyr::select(rsid=markers, chr, pos) 
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
final_out <-merge(ex_dist, out, by.x="markers", by.y="names", all.x=T) %>% 
  distinct(markers, .keep_all=T) %>% arrange(length) %>% 
  mutate(all_genes=paste(PRECEDESYMBOL, GENESYMBOL, FOLLOWSYMBOL) ) %>% 
  dplyr::select(markers, taxa, length, chr, peak.pos, all_genes, LOCATION, closest_gene, 
                closest_gene_type, list_total_genes, closest_gene_description) %>% 
  mutate(all_genes=gsub(x=all_genes,pattern = "NA", replacement = "")) %>% 
mutate(all_genes=gsub(x=trimws(all_genes),pattern = "  ", replacement = ", ")) 

for (n in 1:nrow(final_out)){
  all_tax<- unique(trimws(unlist(str_split(final_out$taxa[n], "\\|"))))
  final_out$taxa[n] <- paste(all_tax, collapse = " | ")
}
write_delim(final_out, "genes_1Mb_intervals_closest_genes_pval_corr.csv", delim=";" )

