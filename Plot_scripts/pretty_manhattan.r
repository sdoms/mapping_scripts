### test ####
# df <- phylum2
# SNP <- "marker"
# CHR <- "chr"
# BP <- "pos"
# P <- "P_Proteobacteria.gt.P"
# col <- mypalette
# title <- "P_Proteobacteria * gt "
# hlight <- c("JAX00682224", "UNCHS004753")
# threshold <- 1e-6
# sig = 0.05/nrow(df) # significant threshold line
# sugg = 1e-6 # suggestive threshold line
# q <-qvalue(as.numeric(gwscan[,P]))
# threshold <- -log10(max(na.omit(q$pvalues[q$qvalues <=0.1])))
# if is.na(threshold){
#   threshold <-NA
# }

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette


gg.manhattan <- function(gwscan, SNP="marker", CHR="chr", BP="pos", P="P", threshold=NA, hlight=NA, col=mypalette, title=""){
  require(readr)
  require(ggrepel)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(RColorBrewer)
  #require(qvalue)
  # make the df 
  gwscan<- gwscan %>% drop_na(P)
  df <- data.frame(SNP = gwscan[,SNP],CHR=as.numeric(gwscan[,CHR]), BP=as.numeric(gwscan[,BP]) * 1e6, P=as.numeric(gwscan[,P]))
  df <- df[order(df$CHR),]
  
  
  sig = 0.05/nrow(df) # significant threshold line
  sugg = 1/nrow(df)  # suggestive threshold line
  if (is.na(threshold)){
      threshold <- sig
    }


  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0)) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    
    # Add highlighted points

    
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
    
    # highlight significant points 
    geom_point(data=subset(df.tmp, is_annotate=="yes"), color="#F6DD66", size=2) +
    geom_point(data=subset(df.tmp, is_highlight=="yes"), color="blue", size=2) +
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
