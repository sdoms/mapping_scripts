get_closest_genes <- function (df){
  require(VariantAnnotation)
  require(org.Mm.eg.db)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  #require(BSgenome.Mmusculus.UCSC.mm10)
  require(tidyverse)
  # test
  # df<- sig.intervals.comb
  # df_with_ref<- left_join(df, preLDpruneSNPs, by=c("peak.snp"="marker"))

  # closest gene to peak snp
  target2 <- with(df,
                  GRanges( seqnames = Rle(paste0("chr",chr)),
                           ranges   = IRanges(peak.pos, end=peak.pos+1, names=names),
                           strand   = NULL ))
  #values(target2) <- DataFrame(varAllele = df_with_ref$A1)
    loc2 <- locateVariants(target2, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants())
  names(loc2)<- NULL
  
  # preceding
  p_ids <- unlist(loc2$PRECEDEID, use.names=FALSE) 
  exp_ranges <- rep(loc2,  elementNROWS(loc2$PRECEDEID))
  p_dist <- GenomicRanges::distance(exp_ranges, TxDb.Mmusculus.UCSC.mm10.knownGene, id=p_ids, type="gene")
  exp_ranges$PRECEDE_DIST <- p_dist
  ## Collapsed view of ranges, gene id and distance:
  loc2$PRECEDE_DIST <- relist(p_dist, loc2$PRECEDEID)
  # following
  f_ids <- unlist(loc2$FOLLOWID, use.names=FALSE) 
  exp_ranges_f <- rep(loc2,  elementNROWS(loc2$FOLLOWID))
  f_dist <- GenomicRanges::distance(exp_ranges_f, TxDb.Mmusculus.UCSC.mm10.knownGene, id=f_ids, type="gene")
  exp_ranges_f$FOLLOW_DIST <- f_dist
  ## Collapsed view of ranges, gene id and distance:
  loc2$FOLLOW_DIST <- relist(f_dist, loc2$FOLLOWID)
  # find closest
  loc2$PRECEDEID_place <-lapply(loc2$PRECEDE_DIST, function(x) which.min(x))
  loc2$PRECEDEID_place <- as.numeric(as.character(loc2$PRECEDEID_place))
  loc2$FOLLOWID_place <-lapply(loc2$FOLLOW_DIST, function(x) which.min(x))
  loc2$FOLLOWID_place <- as.numeric(as.character(loc2$FOLLOWID_place))
  out2 <- as.data.frame(loc2)
  out2$names <- names(target2)[ out2$QUERYID ]
  precede_ID <- data.frame(out2$PRECEDEID)
  colnames(precede_ID)<- "ID"
  precede_ID$place <- out2$PRECEDEID_place
  for (x in 1:nrow(precede_ID)){
    n <- precede_ID[x,"place"]
    precede_ID[x,"closest"] <- as.vector(precede_ID[x,1])[[1]][n]
  }
  follow_ID <- data.frame(out2$FOLLOWID)
  colnames(follow_ID)<- "ID"
  follow_ID$place <- out2$FOLLOWID_place
  for (x in 1:nrow(follow_ID)){
    n <- follow_ID[x,"place"]
    follow_ID[x,"closest"] <- as.vector(follow_ID[x,1])[[1]][n]
  }
  out2 <- out2[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID")]
  out2 <- cbind(out2, precede_ID$closest)
  out2 <- cbind(out2, follow_ID$closest)
  out2 <- unique(out2)
  colnames(out2)<- c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "PRECEDEID", "FOLLOWID")
  out2$GENESYMBOL <- id2Symbol[ as.character(out2$GENEID) ]
  out2$PRECEDESYMBOL <- id2Symbol[ as.character(out2$PRECEDEID) ]
  out2$FOLLOWSYMBOL <- id2Symbol[ as.character(out2$FOLLOWID) ]
  closest_gene_out <- out2 %>% 
    mutate(closest_genes=paste(PRECEDESYMBOL, GENESYMBOL, FOLLOWSYMBOL) ) %>% 
    mutate(closest_genes=gsub(x=closest_genes,pattern = "NA", replacement = "")) %>% 
    mutate(closest_genes=gsub(x=trimws(closest_genes),pattern = "  ", replacement = ", ")) %>% 
    dplyr::select(-GENEID, -PRECEDEID, -FOLLOWID, -FOLLOWSYMBOL, -PRECEDESYMBOL) 
  return(closest_gene_out)

  }
  
get_genes_interval <- function(df){
  df <- df %>% 
    mutate(new.start=ifelse(start.pos==stop.pos, start.pos-1e6, start.pos),
           new.stop=ifelse(start.pos==stop.pos, stop.pos+1e6, stop.pos))
  target <- with(df,
                 GRanges( seqnames = Rle(paste0("chr",chr)),
                          ranges   = IRanges(new.start, end=new.stop, names=names),
                          strand   = NULL ))
  
    loc <- locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants())
    
 
  names(loc)<- NULL
  out <- as.data.frame(loc)
  out$names <- names(target)[ out$QUERYID ]
  #out<- out[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "QUERYID")]
  
  # sig regions
  out<- out %>% distinct(GENEID, QUERYID, .keep_all=T) %>% drop_na(GENEID)
  
  
  
  x <- unique( with(out, c(levels(as.factor(GENEID)))))
  table( x %in% names(id2Symbol)) # good, all found
  out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
  all_genes_out <- out %>% 
    dplyr::group_by(QUERYID) %>% 
    mutate(number_of_genes=length(unique(GENESYMBOL)), all_genes= paste0(unique(GENESYMBOL), collapse= ", ")) %>% 
    ungroup() %>% 
    distinct(names,  number_of_genes, all_genes)
  return(all_genes_out)
}
