# 
# trait <- "G_Roseburia"
# chr <- 15
# start <- 33956343
# stop <- 35275078
# peak_snp <- "UNCHS039968"
# dnaorrna <- "DNA"
# tax_level<- "genus"
plot.region <- function(trait,tax_level, chr, start, stop, peak_snp,P="P",dnaorrna, outputdir=".", gene.ensembl){
  require(tidyverse)
  require(biomaRt)
  require(patchwork)
  require(ggsci)
  begin <- as.numeric(as.character(start))
  end<- as.numeric(as.character(stop))
  load("~/Documents/PhD/Experiments/QTL_mapping_results/Genotype_data/GM_snps.Rdata")
  gwscan<- readRDS(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/", dnaorrna, 
                          "/",tax_level,"/",trait, "_chr_", chr, "with_add_dom.rds"))
  system(paste0("plink --bfile ~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Summary_tables/hybrid_f2_5 -r2 --ld-snp ",
                peak_snp," --ld-window-kb 10000 --ld-window-r2 0 --ld-window 99999 --allow-extra-chr"))
  ld_set <- read.table("plink.ld", header=T)
  ss <- GM_snps  %>% dplyr::filter(between(pos*1e6, start-10, stop+10)) # GM snps is sometimes rounded to 1e5
  max_ld <- merge(ss, ld_set, by.x="marker", by.y="SNP_B")
  
  para.region <- gwscan %>% filter(between(pos*1e6, start-10, stop+10))
  para.region2 <- merge(para.region,max_ld, by="marker",all.y=2 )
  para.region2$R2 <- as.numeric(para.region2$R2)
  
  # highlight lead snp
  para.region2$highlight <- ifelse(para.region2$marker == peak_snp, "highlight", 
                                   ifelse(is.na(para.region2[,P]), "No P value", "normal"))
  para.region2$pos<- para.region2$pos.y*1e6
  para.region2[is.na(para.region2[,P]), P]<- 1
  para.region2$log10p <- -log10(para.region2[,P])
  title <- paste(trait, dnaorrna, peak_snp)
  p1 <- ggplot() + 
    geom_point(data = para.region2,aes(pos, log10p,shape=highlight, colour = R2, fill=highlight), size=2.5, stroke=1) + 
    labs(title = title) + scale_color_viridis_c(option="inferno", direction = -1, limits=c(0,1)) + scale_shape_manual(values=c(23,1, 16)) +
    labs(y="-log10(P)")+theme_bw() +scale_fill_manual(values= c("#7fffd7", "white", "white"))
  p1
  if (stop-start>2e6){
    out.bm.genes.region <- getBM(
      attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', "description"), 
      filters = c('chromosome_name','start','end', 'biotype'), 
      values = list(chr,start, stop, "protein_coding"), 
      mart = gene.ensembl)
  }else if (start==stop){
    start= start -5e5
    stop=stop+5e5
    out.bm.genes.region <- getBM(
      attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', "description"), 
      filters = c('chromosome_name','start','end', 'biotype'), 
      values = list(chr,start, stop, "protein_coding"), 
      mart = gene.ensembl)
  }else{
    out.bm.genes.region <- getBM(
      attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', "description"), 
      filters = c('chromosome_name','start','end'), 
      values = list(chr,start, stop), 
      mart = gene.ensembl)
  }
    if (nrow(out.bm.genes.region)<1){
      repeat {
        start <- start -1e5
        stop <- stop +1e5
        out.bm.genes.region <- getBM(
          attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', "description"), 
          filters = c('chromosome_name','start','end', 'biotype'), 
          values = list(chr,start, stop, "protein_coding"), 
          mart = gene.ensembl)
        if (nrow(out.bm.genes.region)>0){
          break
        }
      }
      
    }else {
      out.bm.genes.region<- out.bm.genes.region %>% filter(gene_biotype%in%c("protein_coding", "lincRNA", "miRNA", "rRNA", 
                                                                             "scRNA", "snRNA", "snoRNA"))
    }
  
  
  
  ## rank gene names according to start position
  out.bm.genes.region <- out.bm.genes.region %>% mutate(external_gene_name = fct_reorder(external_gene_name, 
                                                                                         start_position, .desc = TRUE))


  
 
  
  ## rank gene_biotype label
  out.bm.genes.region <- out.bm.genes.region %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype), 
                                                                                       "protein_coding"), external_gene_name = fct_reorder2(external_gene_name, 
                                                                                                                                            start_position, gene_biotype_fac, .desc = TRUE))
  
  ## plot
  p2 <- ggplot(data = out.bm.genes.region) + 
    geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, 
                       group = gene_biotype_fac)) +theme_bw()+ ylim(start, stop)+
    coord_flip() + ylab("") +
    geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), 
              fontface = 2, alpha = I(0.7), hjust = "right", size= 3) + 
    labs(title = "", subtitle = "", caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          strip.text.y = element_text(angle = 0),
          legend.position="bottom", 
          panel.grid.major.y = element_blank()) + 
    expand_limits(y=c(-1, 1)) + scale_color_d3("category20c")
   
p2
  p1b <- p1 + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
    xlim(start-1,stop+1) +geom_vline(xintercept = begin, linetype="dotted") + geom_vline(xintercept = end, linetype="dotted") +
    theme(legend.direction = "vertical", legend.box = "horizontal")
  
  p1b / p2 
  ggsave(paste0(outputdir,"/" ,trait, "-", P,"-chr", chr, "_",  peak_snp ,"_", dnaorrna,"_region_plot.pdf"), 
         width = 9, height =9, dpi=300)
  
}

