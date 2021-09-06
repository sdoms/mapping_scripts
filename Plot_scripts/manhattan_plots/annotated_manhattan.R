# setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/otu/")
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/otu/")
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/genus/")
trait <- "SV184"
taxa <- data.frame()
threshold<- NA
Xtrue<-T
hlight<-NA
for (chr in 1:19){
  tx <- readRDS(paste0(trait, "_chr_", chr, "with_add_dom.rds"))
  taxa <- rbind(taxa, tx)
}
if (Xtrue){
  X_chrom <- readRDS(paste0(trait, "_chrX.rds"))
  X_chrom$chr <- "20"
  taxa_X <- dplyr::bind_rows(taxa, X_chrom)
} else{
  taxa_X <- taxa
}
if (!is.na(threshold)){
  threshold <- 0.05/nrow(taxa_X)/118.2177
}

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
require(patchwork)
require(gg.gap)

p1 <-gg.manhattan(gwscan = taxa_X, SNP="marker", CHR = "chr", BP="pos", P="P",threshold=threshold, col = mypalette, title = "", method="bon")

p2<- gg.manhattan(taxa,SNP="marker", CHR = "chr", BP="pos", P="add.P",threshold=threshold, col = mypalette, title = "", method="bon")


p3<-gg.manhattan(taxa,SNP="marker", CHR = "chr", BP="pos", P="dom.P", threshold=threshold,col = mypalette, title ="", method="bon")
p <- p1 /(p2+p3)
p
ggsave("ASV293_annotation.pdf")
ggsave("G_Dorea_annotation.pdf")
ggsave("ASV184.pdf")
# gapped Dorea
p1b<-gg.gap(p1, segments=list(c(23,30), c(35,50)), ylim = c(0,60), tick_width = c(1,1,2), rel_heights = c(1,0,0.2,0,0.2))
 
p1b
ggsave("G_Dorea_p1.ggrepel.pdf", plot=p1b, height = 15, width=16)
p2b<-gg.gap(p2, segments=list(c(27,42), c(50,80)), ylim = c(0,90), tick_width = c(1,2,2), rel_heights = c(1,0,0.2,0,0.2))
p2b
ggsave("G_Dorea_p2.ggap.pdf", plot=p2b, height = 15, width=16)
ggsave("G_Dorea_p2.ggap.png", plot=p2b, dpi=300, height = 15, width=16)

p3
p3b<-gg.gap(p3, segments=list(c(27,42), c(50,78)), ylim = c(0,90), tick_width = c(1,2,2), rel_heights = c(1,0,0.2,0,0.2))
p3b
ggsave("G_Dorea_p3.ggap.pdf", plot=p3b, height = 15, width=16)
ggsave("G_Dorea_p3.ggap.png", plot=p3b, dpi=300, height = 15, width=16)

ggsave("../manhattan_plots/annotated_ASV36_RNA_P.pdf", plot=p1, height = 15, width=16)

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
# gapped ASV293
p1b<-gg.gap(p1, segments=list(c(25,34), c(37,50)), ylim = c(0,60), tick_width = c(1,2,2), rel_heights = c(1,0,0.2,0,0.2))

p1b
ggsave("ASV293_p1.ggrepel.pdf", plot=p1b, height = 15, width=16)
p2b<-gg.gap(p2, segments=list(c(27,42), c(50,80)), ylim = c(0,90), tick_width = c(1,3,2), rel_heights = c(1,0,0.2,0,0.2))
p2b
ggsave("ASV293_p2.ggap.pdf", plot=p2b, height = 15, width=16)
ggsave("ASV293_p2.ggap.png", plot=p2b, dpi=300, height = 15, width=16)

p3
p3b<-gg.gap(p3, segments=list(c(27,42), c(50,84)), ylim = c(0,95), tick_width = c(1,2,2), rel_heights = c(1,0,0.2,0,0.2))
p3b
ggsave("ASV293_p3.ggap.pdf", plot=p3b, height = 15, width=16)
ggsave("ASV293_p3.ggap.png", plot=p3b, dpi=300, height = 15, width=16)

# gapped ASV184
p1
p1b<-gg.gap(p1, segments=list(c(27,39), c(41,62)), ylim = c(0,70), tick_width = c(1,2,2), rel_heights = c(1,0,0.2,0,0.2))

p1b
ggsave("ASV184_p1.ggrepel.pdf", plot=p1b, height = 15, width=16)
p2
p2b<-gg.gap(p2, segments=list(c(20,31),c(32,55), c(56,112)), ylim = c(0,120), tick_width = c(1,1,1,2), rel_heights = c(1,0,0.1,0,0.1,0,0.2)) 
p2b
ggsave("ASV184_p2.ggap.pdf", plot=p2b, height = 15, width=16)
ggsave("ASV184_p2.ggap.png", plot=p2b, dpi=300, height = 15, width=16)
p3b<-gg.gap(p3, segments=list(c(22,28),c(30,53), c(55,104)), ylim = c(0,115), tick_width = c(1,1,1,2), rel_heights = c(1,0,0.2,0,0.2,0,0.2)) 
p3b

ggsave("ASV184_p3.ggap.pdf", plot=p3b, height = 15, width=16)
ggsave("ASV184_p3.ggap.png", plot=p3b, dpi=300, height = 15, width=16)


# ggsave("../manhattan_plots/annotated_ASV36_RNA_P.pdf", plot=p1, height = 15, width=16)

gg.manhattan <- function(gwscan, SNP="marker", CHR="chr", BP="pos", P="P", threshold=NA, hlight=NA, col=mypalette, title="", method="bon"){
  require(readr)
  require(ggrepel)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(RColorBrewer)
  require(qvalue)
  snp_annotations <- read_delim("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/snps_gene_annotation.csv", delim=";")
  goodChrOrder<- c(1:19, "X")
  # make the df 
  gwscan<- gwscan %>% drop_na(P)
  df <- data.frame(SNP = gwscan[,SNP],CHR=gwscan[,CHR], BP=gwscan[,BP] * 1e6, P=gwscan[,P])
  colnames(df)<- c("SNP", "CHR", "BP", "P")
  df$CHR <- factor(df$CHR, levels=goodChrOrder)
  # df <- df[order(df$CHR),]
  
  if (method=="bon"){
    sig <- 0.05/nrow(df)
    #sugg = 1/nrow(df) 
  } else if (method =="fdr"){
    #threshold
    q <-qvalue(df$P)
    sig <- (max(na.omit(q$pvalues[q$qvalues <=0.1])))
  }
  
  #sig = 0.05/nrow(df) # significant threshold line
  sugg = 1/nrow(df)  # suggestive threshold line
  if (is.na(threshold)){
    threshold <- sig
    sugg = 1/nrow(df)  # suggestive threshold line
  } else{
    sugg=0.05/nrow(df)
    
  }
  
  topP <- df %>%  filter(P<threshold) %>% group_by(CHR) %>%  slice_min(order_by = P, n = 5)
  onetopP <- df %>%  filter(P<threshold) %>% group_by(CHR) %>%  slice_min(order_by = P, n = 2)
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% topP$SNP, "yes", "no")) %>%
    mutate( is_tophighlight=ifelse(SNP %in% onetopP$SNP, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no")) %>% 
    mutate(is_sugg = ifelse(P<sugg, "yes", "no")) %>% 
    left_join(snp_annotations, by=c("SNP"="marker"))
  
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) %>%  arrange(CHR)
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR), alpha=-log10(P)), size=1) +
    
    scale_alpha(range = c(0.4,1))+
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0)) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(threshold)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points
    geom_text_repel(data=df.tmp[df.tmp$is_tophighlight=="yes",], aes(label=gene), alpha=1, size=3, 
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines"),
                     min.segment.length = unit(0, 'lines'),
                     max.overlaps = 10) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # highlight significant points 
    # geom_point(data=subset(df.tmp, is_sugg=="yes"), color="#70d3e2", size=1) +
    # geom_point(data=subset(df.tmp, is_annotate=="yes"), color="#70e2b8", size=1) +
    # geom_point(data=subset(df.tmp, is_highlight=="yes"), color="#b870e2", size=1) +
    geom_point(data=subset(df.tmp, is_sugg=="yes"), color="#F6DD66", size=1) +
    geom_point(data=subset(df.tmp, is_annotate=="yes"), color="#70e27f", size=1) +
    # geom_point(data=subset(df.tmp, is_highlight=="yes"), color="#b870e2", size=1) +
    coord_cartesian(clip = 'off')+
    # Custom the theme:
    theme_bw(base_size = 12) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}
