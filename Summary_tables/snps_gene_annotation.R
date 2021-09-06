load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata")
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/")
library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db) 
library(tidyverse)
snps <- snps %>%  filter(chr %in% c(1:19, "X")) %>%  mutate(chr=paste0("chr", chr))
target <- with(snps,
               GRanges( seqnames = Rle(chr),
                        ranges   = IRanges(start=pos*1e6, end=(pos*1e6+1), names = marker ),
                        strand   = NULL) )

loc <- locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants(intergenic =IntergenicVariants(upstream=1.5e6, downstream=1.5e6)),ignore.strand=T)
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
for (i in 1:nrow(out)){
  out$gene[i] <- paste(out$PRECEDESYMBOL[i],out$GENESYMBOL[i], out$FOLLOWSYMBOL[i], collapse = ",")
  
}
final_out <-
  snps %>% mutate(start=pos*1e6) %>% 
  left_join(out,by=c("marker"= "names"))
final_out <- final_out %>% 
  dplyr::select(marker, chr, pos, gene, LOCATION)
write_delim(final_out, "snps_gene_annotation.csv", delim=";")

