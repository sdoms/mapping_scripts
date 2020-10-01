plot.pretty.region <- function(gws, Xtrue, P="P", title=""){
  require(biomaRt)
  require(stringr)
  require(forcats)
  require(RColorBrewer)
  require(viridis)
  require(plyr)
  require(patchwork)
  require(tidyverse)

  
  
  # Biomarts
  snp.ensembl <- useEnsembl(biomart = "snp", dataset="mmusculus_snp", mirror="www")
  gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = 'www') # we will need an additional mart for genes
  gwscan <- data.frame()
  chrom_range <- c(1:19)
  for (chr in 1:19){
    tx <- readRDS(paste0(gws, "_chr_", chr, "with_add_dom.rds"))
    gwscan <- rbind(gwscan, tx)
  }
  if (Xtrue){
    X_chrom <- readRDS(paste0(gws, "_chrX.rds"))
    #X_chrom$chr <- "20"
    gwscan <- dplyr::bind_rows(gwscan, X_chrom)
    chrom_range <- c(1:19, "X")
    
}

  

  rownames(gwscan)<- gwscan$marker
  
  # Make the result directories
  dir.create("../figures", showWarnings = FALSE)
  dir.create("../genes", showWarnings = FALSE)
  dir.create("../results", showWarnings = FALSE)
  
  

  # Bonferonni threshold 
  # threshold <- -log10(0.05/nrow(gwscan))
  threshold <- -log10((0.05/nrow(gwscan))/118.2177)
  
  #Initializing parameters 
  intervals <- c()
  min_max_positions <- c(paste("chr", "bin", "start", "end", "\n"))
  out <- data.frame()
  #names(out)<- c("name", "chr", "bin", "start", "stop", "length", "significant_snps","peak_snp","position_peak_snp","markers", "total_genes", "list_total_genes", "refseq_genes", "list_refseq_genes" )
  
  # for each chromosome 
  for (chr in chrom_range){
    # only for this chromosome
    cat("Working on chromosome ", chr, ".\n")
    gwscan_chr <- gwscan[gwscan$chr==chr,]
    gwscan_chr <-na.omit(gwscan_chr, cols=P)
    gwscan_chr$log10p <- -log10(gwscan_chr[,P])

    if (length(which(gwscan_chr$log10p>=threshold))>=1){
      significant_snps <- gwscan_chr[which(gwscan_chr$log10p>=threshold),] 
      #0.9 ld -> interesting snps 1 -> pairwise ld for before, interesting snps last -> pairwise ld of after. 
      # then use the ld 0.9 threshold to determine the plotting region. --> eLIFE 
      load("~/Documents/PhD/Experiments/QTL_mapping_results/Genotype_data/GM_snps.Rdata")
      
      # combine significant regions < 10 Mb apart into a single region
      bins <- split(significant_snps,cut_interval(significant_snps$pos, length=10))
      
      j=0
      for (i in 1:length(bins)){
        if (nrow(bins[[i]])>0){
          j=j+1
          min_pos <-max(gwscan_chr$pos)
          max_pos <- min(gwscan_chr$pos)
          for (snp in rownames(bins[[i]])){
            system(paste0("plink --bfile ~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/hybrid_f2 -r2 --ld-snp ",snp," --ld-window-kb 10000 --ld-window-r2 0.9 --ld-window 99999 --allow-extra-chr"))
            ld_set <- read.table("plink.ld", header=T)
            rownames(ld_set)<- ld_set$SNP_B
            ld_set$marker <- ld_set$SNP_B
            ld_set<-match_df(GM_snps, ld_set, on="marker")
            min_pos <- min(c(ld_set$pos, min_pos)) 
            max_pos <- max(c(ld_set$pos, max_pos))
          }  
          intervals<- c(intervals,(max_pos-min_pos))
          # if (nrow(bins[[i]])==1){
          #   min_pos <- min_pos -0.5
          #   max_pos <- max_pos +0.5
          # }
          gwscan_region <- gwscan_chr[which(gwscan_chr$pos >= min_pos-0.5 & gwscan_chr$pos <= max_pos+0.5 ),]
          
          #max snp
          max_snp <- rownames(gwscan_region[which(gwscan_region$log10p==max(gwscan_region$log10p)),])[1]
        
         
          min_pos <- min_pos*1e6
          max_pos <- max_pos*1e6
          
          min_max_positions<-c(min_max_positions, paste(chr,j, min_pos, max_pos, "\n"))
          #get rsID
          rsId <- droplevels(GM_snps[max_snp, "rsID"][[1]])
          if (is.na(rsId)){
            rsId <- getBM(attributes = c('refsnp_id'), filters = c('chr_name', 'start', 'end'), values = list(chr, GM_snps[max_snp,"pos"]*1e6 , GM_snps[max_snp,"pos"]*1e6+1), mart = snp.ensembl)[1,1]
          }
          
          if (rsId=="logical(0)" | is.na(rsId)){ # when there is no rsId
            out.bm.snp2gene <- getBM(attributes = c('ensembl_gene_id',"chromosome_name", "start_position", "end_position"), filters = c("start", "end", "chromosome_name"), values=list(min_pos, max_pos, chr), mart=gene.ensembl)
            if (dim(out.bm.snp2gene)[1]==0){
              sel.chr <- chr
              sel.pos<- GM_snps[max_snp,"pos"]*1e6
            } else{
              out.bm.gene <- tryCatch(getBM(attributes = c('external_gene_name', 'gene_biotype', 'description'), 
                                            filters = c('ensembl_gene_id'),
                                            values = unique(out.bm.snp2gene$ensembl_gene_id), 
                                            mart = gene.ensembl), error=function(x) return(data.frame(external_gene_name=NA, gene_biotype=NA, description=NA)))
              sel.chr <- out.bm.snp2gene$chromosome_name
              sel.pos <- out.bm.snp2gene$start_position
              }
            
            

          } else{
            out.bm <- getBM(
              attributes = c("ensembl_gene_stable_id", 
                             "refsnp_id", 
                             "chr_name", 
                             "chrom_start", 
                             "chrom_end", 
                             "minor_allele", 
                             "minor_allele_freq"),
              # "ensembl_transcript_stable_id",
              # "consequence_type_tv"),
              filters = "snp_filter", 
              values = rsId,#dat.bmi.sel$SNP, 
              mart = snp.ensembl)
            if ((dim(out.bm)[1] == 0)){
              rsId <- as.character(getBM(attributes = c('refsnp_id'), filters = c('chr_name', 'start', 'end'), values = list(chr, GM_snps[max_snp,"pos"]*1e6 , GM_snps[max_snp,"pos"]*1e6+1), mart = snp.ensembl))
              out.bm <- getBM(
                attributes = c("ensembl_gene_stable_id", 
                               "refsnp_id", 
                               "chr_name", 
                               "chrom_start", 
                               "chrom_end", 
                               "minor_allele", 
                               "minor_allele_freq"),
                # "ensembl_transcript_stable_id",
                # "consequence_type_tv"),
                filters = "snp_filter", 
                values = rsId,#dat.bmi.sel$SNP, 
                mart = snp.ensembl)
            }
            
            #listAttributes(snp.ensembl) %>% slice(str_which(name, "gene")) 
            out.bm.snp2gene <- getBM(
              attributes = c('refsnp_id', 'allele', 'chrom_start', 'chr_name', 'ensembl_gene_stable_id'), 
              filters = c('snp_filter'), 
              values = rsId, 
              mart = snp.ensembl)
            
            out.bm.gene <- tryCatch(getBM(attributes = c('external_gene_name', 'gene_biotype', 'description'), 
                                 filters = c('ensembl_gene_id'),
                                 values = unique(out.bm.snp2gene$ensembl_gene_stable_id), 
                                 mart = gene.ensembl), error=function(e) return(data.frame(external_gene_name=NA, gene_biotype=NA, description=NA)))
            sel.chr <- out.bm.snp2gene$chr_name
            sel.pos <- out.bm.snp2gene$chrom_start
          }
          
          
          out.bm.gene
          
          range <- 1e5

          
          para.region <- gwscan_chr %>% filter(between(pos*1e6, min_pos-range, max_pos+range))
    
          # ld to max snp
          system(paste0("plink --bfile ~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/hybrid_f2 --r2 --ld-snp ",max_snp," --ld-window-kb 10000 --ld-window-r2 0 --ld-window 99999 --allow-extra-chr --out max_snp"))
          max_ld <- read.table("max_snp.ld", header=T)
          ss <- GM_snps  %>% filter(between(pos*1e6, min_pos-range, max_pos+range))
          max_ld2 <- merge(ss, max_ld, by.x="marker", by.y="SNP_B")
          para.region2 <- merge(para.region,max_ld2, by="marker",all.y=2 )
          para.region2$R2 <- as.numeric(para.region2$R2)
          
          # highlight lead snp
          para.region2$highlight <- ifelse(para.region2$marker == max_snp, "highlight", ifelse(is.na(para.region2[,P]), "No P value", "normal"))
          para.region2$pos<- para.region2$pos.y*1e6
          para.region2[is.na(para.region2[,P]), P]<- 1
          para.region2$log10p <- -log10(para.region2[,P])
         
          
          p1 <- ggplot() + 
            geom_point(data = para.region2,aes(pos, log10p,shape=highlight, colour = cut(R2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf)))) + 
            labs(title = title) + scale_color_brewer("R2", palette = "Set1") + scale_shape_manual(values=c(8,1, 16)) +labs(y="-log10(P)")
          p1
          
          out.bm.genes.region <- getBM(
            attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', "description"), 
            filters = c('chromosome_name','start','end'), 
            values = list(sel.chr, min_pos - range, max_pos + range), 
            mart = gene.ensembl)
          head(out.bm.genes.region)
          
          if ((dim(out.bm.genes.region)[1] == 0)){
            res <- data.frame(name=paste(gws,P,chr, j, sep = "_"), 
                              chr=chr, bin=j,P=P, start=min_pos, stop=max_pos, length=(max_pos -min_pos)/1e6, 
                              significant_snps=nrow(bins[[i]]),peak_snp=max_snp, rsId= rsId, 
                              position_peak_snp=significant_snps[max_snp, "pos"], P_peak_snp=gwscan_region[max_snp, "P"], 
                              add.Beta_peak_snp=gwscan_region[max_snp, "add.Beta"],add.StdErr_peak_snp=gwscan_region[max_snp, "add.StdErr"],add.Zscore_peak_snp=gwscan_region[max_snp, "add.T"],add.P_peak_snp=gwscan_region[max_snp, "add.P"],
                              dom.Beta_peak_snp=gwscan_region[max_snp, "dom.Beta"],dom.StdErr_peak_snp=gwscan_region[max_snp, "dom.StdErr"],dom.Zscore_peak_snp=gwscan_region[max_snp, "dom.T"],dom.P_peak_snp=gwscan_region[max_snp, "dom.P"],
                              markers=paste(bins[[i]]$marker, collapse = ", "),
                              closest_gene="", closest_gene_type= "", 
                              closest_gene_description= "", 
                              total_genes=0, 
                              list_total_genes="", 
                              total_protein_coding_genes= 0,
                              list_protein_coding_genes="")
            out <- rbind(out, res)
            next()
          }
          ## rank gene names according to start position
          out.bm.genes.region <- out.bm.genes.region %>% mutate(external_gene_name = fct_reorder(external_gene_name, 
                                                                                                 start_position, .desc = TRUE))
          
          ## plot
          ggplot(data = out.bm.genes.region) + 
            geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position)) + 
            coord_flip() + 
            ylab("")
          
          ## define plot range for x-axis
          plot.range <- c(min(sel.pos - range, out.bm.genes.region$start_position), 
                          max(sel.pos + range, out.bm.genes.region$end_position))
          
          ## rank gene_biotype label
          out.bm.genes.region <- out.bm.genes.region %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype), 
                                                                                               "protein_coding"), external_gene_name = fct_reorder2(external_gene_name, 
                                                                                                                                                    start_position, gene_biotype_fac, .desc = TRUE))
          
          ## plot
          p2 <- ggplot(data = out.bm.genes.region) + 
            geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
            coord_flip() + ylab("") +
            ylim(plot.range) + 
            geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
            labs(title = "", subtitle = paste0("Genes"), caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(), 
                  strip.text.y = element_text(angle = 0),
                  legend.position="bottom", 
                  panel.grid.major.y = element_blank()) + 
            expand_limits(y=c(-1, 1)) +
            scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1))) 
          ## hack to have 11 colors, but probably merging some gene biotype groups could make sense too. 
          ## consider colorblindr::palette_OkabeIto_black
          print(p2)
          p1b <- p1 + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlim(plot.range) +geom_vline(xintercept = min_pos, linetype="dotted") + geom_vline(xintercept = max_pos, linetype="dotted") +theme(legend.direction = "vertical", legend.box = "horizontal")
          
          p1b + p2 + plot_layout(ncol = 1, heights = c(6, 6)) 
          ggsave(paste0("../figures/", gws, "-", P,"-chr", chr, "_bin",  j ,"_region_plot.pdf"), width = 9, height =6)
          
          pos_max_snp <- significant_snps[max_snp, "pos"] *1e6
          if ((dim(out.bm.gene)[1] == 0)){
            absol_diff <- matrix(c(abs(out.bm.genes.region$start_position-pos_max_snp), 
                                     abs(out.bm.genes.region$end_position-pos_max_snp)), ncol=2, byrow = F)
            closest_gene <- out.bm.genes.region[which(absol_diff==min(absol_diff), arr.ind = T)[1],]
            
          } else {
            closest_gene <- out.bm.gene
          }
          
          # write to result file
          res <- data.frame(name=paste(gws,P,chr, j, sep = "_"), 
                            chr=chr, bin=j,P=P, start=min_pos, stop=max_pos, length=(max_pos -min_pos)/1e6, 
                            significant_snps=nrow(bins[[i]]),peak_snp=max_snp, rsId= rsId, 
                            position_peak_snp=significant_snps[max_snp, "pos"], P_peak_snp=gwscan_region[max_snp, "P"], 
                            add.Beta_peak_snp=gwscan_region[max_snp, "add.Beta"],add.StdErr_peak_snp=gwscan_region[max_snp, "add.StdErr"],add.Zscore_peak_snp=gwscan_region[max_snp, "add.T"],add.P_peak_snp=gwscan_region[max_snp, "add.P"],
                            dom.Beta_peak_snp=gwscan_region[max_snp, "dom.Beta"],dom.StdErr_peak_snp=gwscan_region[max_snp, "dom.StdErr"],dom.Zscore_peak_snp=gwscan_region[max_snp, "dom.T"],dom.P_peak_snp=gwscan_region[max_snp, "dom.P"],
                            markers=paste(bins[[i]]$marker, collapse = ", "),
                            closest_gene=paste(closest_gene$external_gene_name, collapse = ", "), closest_gene_type= paste(closest_gene$gene_biotype, collapse=", "), 
                            closest_gene_description= paste(closest_gene$description, collapse = ", "), 
                            total_genes=nrow(out.bm.genes.region), 
                            list_total_genes=paste(out.bm.genes.region$external_gene_name, collapse = ", "), 
                            total_protein_coding_genes= nrow(out.bm.genes.region[which(
                              out.bm.genes.region$gene_biotype=="protein_coding"),]),
                              list_protein_coding_genes=paste(out.bm.genes.region[which(out.bm.genes.region$gene_biotype=="protein_coding"),"external_gene_name"], collapse = ", "))
          out <- rbind(out, res)
          
          }
        
       
       cat("Ready with chromosome", chr, "bin", i, "j", j, "\n.") 
      }
    }
    cat("Ready with chromosome", chr, ".\n")
  }
  write(intervals, file=paste0("../results/",gws,"-", P,"_intervals.txt"))
  write(min_max_positions, file=paste0("../results/",gws, "-", P,"_min_max_positions.txt"))
  #write.table(out, file=paste0("../results/",gws, "-", P,"_results.txt"))
  return(out)
}


