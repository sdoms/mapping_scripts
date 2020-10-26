summary_tables_histo <-
  function(trait,
           covar = F,
           method = "fdr",
           Xtrue = T) {
    require(ggplot2)
    require(patchwork)
    require(dplyr)
    require(qvalue)
    require(VariantAnnotation)
    require(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
    require(org.Mm.eg.db)
    require(tidyr)
    
    dir.create("./summary_tables", showWarnings = F)
    musdom <-
      read.csv(
        "~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv",
        sep = ","
      )
    load(
      "~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/consensus_F0_all.Rdata"
    )
    parents <- as.data.frame(reduced_geno)
    parents$X <- rownames(parents)
    load("~/Documents/PhD/Experiments/Final_QTL_mapping/Genotype_data/GM_snps.Rdata")
    parents2 <-
      merge(GM_snps[, c("marker", "rsID")], parents, by.x = "marker", by.y =
              "X")
    musdom <- merge(musdom, parents2, by.x = "X", by.y = "marker",)
    
    chrom_range <- c(1:19)
    gwscan <- data.frame()
    # Load in gwscan
    if (covar != F) {
      output <- paste0("./summary_tables/",trait, "_", covar,"_", method)
      
      for (chr in 1:19) {
        tx <-
          readRDS(
            paste0(
              trait,
              "/",
              trait,
              "_with_covar_",
              covar,
              "_chr_",
              chr,
              "with_add_dom.rds"
            )
          )
        gwscan <- rbind(gwscan, tx)
      }
      if (file.exists(paste0(trait, "/", trait, "_with_covar_", covar, "_chrX.rds"))) {
        X_chrom <-
          readRDS(paste0(trait, "/", trait, "_with_covar_", covar, "_chrX.rds"))
        #X_chrom$chr <- "20"
        gwscan <- dplyr::bind_rows(gwscan, X_chrom)
        chrom_range <- c(1:19, "X")
        
      }
    } else {
      output <- paste0("./summary_tables/",trait,"_", method)
      for (chr in 1:19) {
        tx <-
          readRDS(paste0(trait, "/", trait, "_chr_", chr, "with_add_dom.rds"))
        gwscan <- rbind(gwscan, tx)
      }
      if (file.exists(paste0(trait, "/", trait, "_chrX.rds"))) {
        X_chrom <- readRDS(paste0(trait, "/", trait, "_chrX.rds"))
        #X_chrom$chr <- "20"
        gwscan <- dplyr::bind_rows(gwscan, X_chrom)
        chrom_range <- c(1:19, "X")
        
      }
      
    }
    
    sig.P <- data.frame()
    sig.add.P <- data.frame()
    sig.dom.P <- data.frame()
    
    # Bonferonni threshold
    # threshold <- 0.05/nrow(gwscan)
    if (method == "fdr") {
      q_P <- qvalue(gwscan[, "P"])
      threshold_P <- (max(na.omit(q_P$pvalues[q_P$qvalues <= 0.1])))
      
      q_add <- qvalue(gwscan[, "add.P"])
      threshold_addP <-
        (max(na.omit(q_add$pvalues[q_add$qvalues <= 0.1])))
      
      q_dom <- qvalue(gwscan[, "dom.P"])
      threshold_domP <-
        (max(na.omit(q_dom$pvalues[q_dom$qvalues <= 0.1])))
      
      
    } else if (method == "bon") {
      # Bonferonni threshold
      threshold_P <- 0.05 / nrow(gwscan)
      threshold_addP <- threshold_P
      threshold_domP <- threshold_P
    }
    
    # select significant markers
    sig_gwscan_P <- gwscan[which(gwscan$P < threshold_P), ]
    sig_gwscan_add.P <- gwscan[which(gwscan$add.P < threshold_addP), ]
    sig_gwscan_dom.P <- gwscan[which(gwscan$dom.P < threshold_domP), ]
    
    # add mus dom
    sig.P <-
      merge(
        sig_gwscan_P,
        musdom,
        by.x = "marker",
        by.y = "X",
        all.x = T
      )
    sig.add.P <-
      merge(
        sig_gwscan_add.P,
        musdom,
        by.x = "marker",
        by.y = "X",
        all.x = T
      )
    sig.dom.P <-
      merge(
        sig_gwscan_dom.P,
        musdom,
        by.x = "marker",
        by.y = "X",
        all.x = T
      )
    
    # all unique significant snsp per trait
    write.table(x = unique(sig.P),
                file = paste0(output, "_all_sig_snps_P.txt"))
    write.table(x = unique(sig.add.P),
                file = paste0(output, "_all_sig_snps_add_P.txt"))
    write.table(x = unique(sig.dom.P),
                file = paste0(output, "_all_sig_snps_dom_P.txt"))
    
    # add markers together anf find unique markers
    all_sig1 <- rbind(sig.P, sig.add.P)
    all_sig <- rbind(all_sig1, sig.dom.P)
    deduped.data <- unique(all_sig)
    
    add <-
      ggplot(deduped.data, aes(x = add.T)) + geom_histogram() + labs(x = "Additive effect Z-score ") + theme(axis.title = element_text(size =
                                                                                                                                         9))
    dom <-
      ggplot(deduped.data, aes(x = dom.T)) + geom_histogram(colour = "black", fill =
                                                              "white") + labs(x = "Dominance effect Z-score") + theme(axis.title = element_text(size =
                                                                                                                                                  9))
    
    add + dom + plot_annotation(title = trait, tag_levels = "A")
    ggsave(paste0(output, "_all_sig_snps-zscore.pdf"))
    
    # write all unique snps
    tt <- deduped.data %>%
      mutate(
        D.A = dom.T / add.T,
        MAF = (AA / n + AA / n + AB / n) / 2,
        chr = as.numeric(as.character(chr)),
        pos = as.numeric(as.character(pos))
      ) %>%
      arrange(chr, pos) %>%
      mutate(chr = replace_na(chr, "X"))
    
    
    write.table(x = tt,
                file = paste0(output, "_all_unique_sig_snps_ALL_P.txt"))
    
    
    input <- tt %>%
      group_by(marker) %>%
      drop_na(chr) %>%
      mutate(chr = paste0("chr", chr), pos = pos * 1e6) %>%
      dplyr::select(rsid = marker, chr, pos) %>%
      distinct()
    target <- with(input,
                   GRanges(
                     seqnames = Rle(chr),
                     ranges   = IRanges(pos, end = pos, names = rsid),
                     strand   = Rle(strand("*"))
                   ))
    
    
    loc <- locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants())
    names(loc) <- NULL
    
    p_ids <- unlist(loc$PRECEDEID, use.names=FALSE) 
    exp_ranges <- rep(loc,  elementNROWS(loc$PRECEDEID))
    p_dist <- GenomicRanges::distance(exp_ranges, TxDb.Mmusculus.UCSC.mm10.knownGene, id=p_ids, type="gene")

    exp_ranges$PRECEDE_DIST <- p_dist
    
    ## Collapsed view of ranges, gene id and distance:
    loc$PRECEDE_DIST <- relist(p_dist, loc$PRECEDEID)

    f_ids <- unlist(loc$FOLLOWID, use.names=FALSE) 
    exp_ranges_f <- rep(loc,  elementNROWS(loc$FOLLOWID))
    
    f_dist <- GenomicRanges::distance(exp_ranges_f, TxDb.Mmusculus.UCSC.mm10.knownGene, id=f_ids, type="gene")

    exp_ranges_f$FOLLOW_DIST <- f_dist
    
    ## Collapsed view of ranges, gene id and distance:
    loc$FOLLOW_DIST <- relist(f_dist, loc$FOLLOWID)

    # find closest
    loc$PRECEDEID_place <-lapply(loc$PRECEDE_DIST, function(x) which.min(x))
    loc$PRECEDEID_place <- as.numeric(as.character(loc$PRECEDEID_place))
    
    loc$FOLLOWID_place <-lapply(loc$FOLLOW_DIST, function(x) which.min(x))
    loc$FOLLOWID_place <- as.numeric(as.character(loc$FOLLOWID_place))
    
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
    colnames(out)<- c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "PRECEDEID", "FOLLOWID")
    Symbol2id <- as.list( org.Mm.egSYMBOL2EG )
    id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
    names(id2Symbol) <- unlist(Symbol2id)
    
    x <- unique( with(out, c(levels(as.factor(GENEID)), levels(as.factor(PRECEDEID)), levels(as.factor(FOLLOWID)))))
    table( x %in% names(id2Symbol)) # good, all found
    
    out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
    out$PRECEDESYMBOL <- id2Symbol[ as.character(out$PRECEDEID) ]
    out$FOLLOWSYMBOL <- id2Symbol[ as.character(out$FOLLOWID) ]

    final_out <-merge(tt, out, by.x="marker", by.y="names", all.x=T)
    write.csv(final_out, paste0(output, "_markers_with_genes.csv"), row.names = F)
    
    
    
  }
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Sterility/")

summary_tables_histo("apo_per_tubule")
summary_tables_histo("arb_per_tubule")
summary_tables_histo("deg_per_tubule")
summary_tables_histo("exf_per_tubule")
summary_tables_histo("mgc_per_tubule")
summary_tables_histo("relative_testis_weight")
summary_tables_histo("sperm_count")
summary_tables_histo("spermatids_per_perimeter")
summary_tables_histo("spermatocytes_per_perimeter")
summary_tables_histo("spermatids_per_tubule")
summary_tables_histo("spermatocytes_per_tubule")
summary_tables_histo("tubule_area")
summary_tables_histo("vac_per_tubule")

summary_tables_histo("combined_testis_weight", covar="body_weight")
summary_tables_histo("sperm_count", covar="combined_testis_weight")
summary_tables_histo("spermatids_per_perimeter", covar="spermatocytes_per_perimeter")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacula")
summary_tables_histo("bac_PC01")
summary_tables_histo("bac_PC02")
summary_tables_histo("bac_PC03")
summary_tables_histo("bac_PC04")
summary_tables_histo("bac_PC05")
summary_tables_histo("bac_PC06")
summary_tables_histo("bac_PC07")
summary_tables_histo("bac_PC08")
summary_tables_histo("bac_PC09")
summary_tables_histo("bac_PC10")
summary_tables_histo("bac_PC11")
summary_tables_histo("bac_PC12")
summary_tables_histo("bac2d_area", covar="body_length")
summary_tables_histo("bac2d_base_width", covar="body_length")
summary_tables_histo("bac2d_length", covar="body_length")
summary_tables_histo("bac2d_shaft_width", covar="body_length")
summary_tables_histo("bac_surf_area", covar="body_length")
summary_tables_histo("bac_vol", covar="body_length")



summary_tables_histo("bac_PC01", method="bon")
summary_tables_histo("bac_PC02", method="bon")
summary_tables_histo("bac_PC03", method="bon")
summary_tables_histo("bac_PC04", method="bon")
summary_tables_histo("bac_PC05", method="bon")
summary_tables_histo("bac_PC06", method="bon")
summary_tables_histo("bac_PC07", method="bon")
summary_tables_histo("bac_PC08", method="bon")
summary_tables_histo("bac_PC09", method= "bon")
summary_tables_histo("bac_PC10", method= "bon")
summary_tables_histo("bac_PC11", method="bon")
summary_tables_histo("bac_PC12", method="bon")










