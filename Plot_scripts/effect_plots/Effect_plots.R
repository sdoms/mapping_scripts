
abundance.plot <- function(taxa, trait, DNAorRNA, marker){
  require(readxl)
  require(ggplot2)
  require(reshape2)
  require(cowplot)
  require(ggpubr)
  
  cbbPalette <- c("#1f78b4",  "#a6cee3","#b2df8a","#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", '#6a3d9a', "#ffff99", "#b15928", "#40e0d0", "#ffff00", "#ee82ee", "#fb8072", "#191970", "#fdb462", "#CBD4D4", "#fccde5", "#d9d9d9", "#bc80bd", "#32cd32", "#ffed6f", "lightcoral" )
  # load("~/Documents/PhD/Experiments/QTL_mapping_results/cleaning/clean_geno.Rdata")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_geno.Rdata")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata")
  theme_set(theme_test())
  
  clean_geno <- as.data.frame(clean_geno)
  colnames(clean_geno) <- gsub(x = colnames(clean_geno), pattern = "\\/", replacement = ".")  
  all_mark <-clean_geno[marker,]
  info_marker <- snps[marker,]
  a <- info_marker$A1 # a1 is the reference allele BB 1
  b <- info_marker$A2 # a2 is the minor allele AA -1 
  genotypes <- t(all_mark)
  
  consensus <-read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv", sep=",")
  allele_freq <-read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/allele_freq_consensus_mus_dom.csv", sep=",")
  rownames(consensus)<- consensus$X 
  consensus$X <- NULL
  mus_dom <- consensus[marker,]
  
  
  # abundances
  #DNAorRNA <- "DNA"
  # otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Miseq_hybrid_mice/DNA_and_RNA/tables_core/",trait, "_table_f2_core.csv"))
  otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/",trait, "_table_f2_core.csv"))
  rownames(otu_table)<- otu_table$X
  otu_table <- otu_table[which(grepl(paste0("_",DNAorRNA),rownames(otu_table))),] # only do the DNA
  rownames(otu_table) <- substr(rownames(otu_table), 1,15)
  
  SV <- as.data.frame(otu_table[,taxa])
  rownames(SV) <- rownames(otu_table)
  tot_with_SV <- merge(genotypes, SV, by="row.names")
  colnames(tot_with_SV)<- c("individuals", "marker", "abundance")
  tot_with_SV$counts <- tot_with_SV$abundance*10000
  tot_with_SV$log_counts <- log10(tot_with_SV$counts+1)
  tot_with_SV$marker<- as.factor(tot_with_SV$marker)
  
  my_comp <- list(c(a, b),c(a, "H"), c(b, "H"))
  
  fig<-ggplot(tot_with_SV, aes(x=marker, y=abundance)) + geom_boxplot(aes(fill=marker))+ 
    stat_compare_means(comparisons = my_comp, label = "p.signif") + theme(legend.position = "none") + 
    labs(x="Genotypes",y="Relative abundance")+ scale_fill_manual(values = cbbPalette)+
    scale_x_discrete(limits=c(a, "H", b),labels=c(paste0(a,a), paste0(a,b), paste0(b,b)))+
    ggtitle(paste0(taxa," abundance for marker ", marker))
  fig <-fig + labs(caption = paste0("Mmm: ", mus_dom[,1], ", Mmd: ", mus_dom[,2]))
  stat_box_data <- function(x, upper_limit = max(tot_with_SV$abundance) * 1.75) {
    return( 
      data.frame(
        y = 0.95 * upper_limit,
        label = paste('count =', 
                      format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE), 
                      '\n',
                      'mean =', 
                      format(round(mean(x), 4), big.mark = ",", decimal.mark = ".", scientific = FALSE))
      )
    )
  }
  fig <-fig + 
    stat_summary(
      fun.data = stat_box_data, 
      geom = "text", 
      hjust = 0.5,
      vjust = 0.9, size=3
    ) 
  fig
  #annotate_figure(fig3, fig.lab = "GG: domesticus, AA: musculus", fig.lab.pos = "bottom.right", fig.lab.size = 12)
  ggsave(file=paste0(taxa,"_marker_", marker,"_",DNAorRNA, ".pdf"), plot=fig)
  
  fig2<-ggplot(tot_with_SV, aes(x=marker, y=log_counts)) + geom_boxplot(aes(fill=marker))+ 
    stat_compare_means(comparisons = my_comp, label = "p.signif")  +theme(legend.position = "none") +
    labs(x="Genotypes",y="Log10 transformed absolute abundance")+ scale_fill_manual(values = cbbPalette)+
    scale_x_discrete(limits=c(a, "H", b),labels=c(paste0(a,a), paste0(a,b), paste0(b,b)))+
    ggtitle(paste0(taxa," abundance for marker ", marker)) 
  #annotate_figure(fig3, fig.lab = "GG: domesticus, AA: musculus", fig.lab.pos = "bottom.right", fig.lab.size = 12)
  ggsave(file=paste0(taxa,"_marker_", marker,"_",DNAorRNA, "log.pdf"), plot=fig2)
  return(fig)
  
}
effect.boxplot <- function(taxa, trait, DNAorRNA, marker){
  require(readxl)
  require(ggplot2)
  require(reshape2)
  require(cowplot)
  require(ggpubr)
  
  cbbPalette <- c("#1f78b4",  "#a6cee3","#b2df8a","#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", '#6a3d9a', "#ffff99", "#b15928", "#40e0d0", "#ffff00", "#ee82ee", "#fb8072", "#191970", "#fdb462", "#CBD4D4", "#fccde5", "#d9d9d9", "#bc80bd", "#32cd32", "#ffed6f", "lightcoral" )
  # load("~/Documents/PhD/Experiments/QTL_mapping_results/cleaning/clean_geno.Rdata")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_geno.Rdata")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata")
  theme_set(theme_test())
  
  clean_geno <- as.data.frame(clean_geno)
  colnames(clean_geno) <- gsub(x = colnames(clean_geno), pattern = "\\/", replacement = ".")  
  all_mark <-clean_geno[marker,]
  info_marker <- genotypes_ref[marker,]
  a <- info_marker$AA_min # a1 is the reference allele BB 1
  b <- info_marker$BB_maj # a2 is the minor allele AA -1 
  genotypes <- t(all_mark)
  
  consensus <-read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv", sep=",")
  allele_freq <-read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/allele_freq_consensus_mus_dom.csv", sep=",")
  rownames(consensus)<- consensus$X 
  consensus$X <- NULL
  mus_dom <- consensus[marker,]
  
  
  # abundances
  #DNAorRNA <- "DNA"
  # otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Miseq_hybrid_mice/DNA_and_RNA/tables_core/",trait, "_table_f2_core.csv"))
  otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/",trait, "_table_f2_core.csv"))
  rownames(otu_table)<- otu_table$X
  otu_table <- otu_table[which(grepl(paste0("_",DNAorRNA),rownames(otu_table))),] # only do the DNA
  rownames(otu_table) <- substr(rownames(otu_table), 1,15)
  
  SV <- as.data.frame(otu_table[,taxa])
  rownames(SV) <- rownames(otu_table)
  tot_with_SV <- merge(genotypes, SV, by="row.names")
  colnames(tot_with_SV)<- c("individuals", "marker", "abundance")
  tot_with_SV$counts <- tot_with_SV$abundance*10000
  tot_with_SV$log_counts <- log10(tot_with_SV$counts+1)
  tot_with_SV$marker<- as.factor(tot_with_SV$marker)
  # 
  # my_comp <- list(c(a, b),c(a, "H"), c(b, "H"))
  # 
  fig<-ggplot(tot_with_SV, aes(x=marker, y=abundance)) + geom_boxplot()+ 
    # stat_compare_means(comparisons = my_comp, label = "p.signif") + 
    theme(legend.position = "none") + geom_jitter(aes(color=marker))+
    labs(x="Genotypes",y="Relative abundance")+ scale_color_manual(values = cbbPalette)+
    scale_x_discrete(limits=c(a, "H", b),labels=c(paste0(a,a), paste0(a,b), paste0(b,b)))+
    ggtitle(paste0(taxa," abundance for marker ", marker))+ scale_y_continuous(labels = scales::percent)
  fig <-fig + labs(caption = paste0("Mmm: ", mus_dom[,1], ", Mmd: ", mus_dom[,2]))
  # stat_box_data <- function(x, upper_limit = max(tot_with_SV$abundance) * 1.75) {
  #   return( 
  #     data.frame(
  #       y = 0.95 * upper_limit,
  #       label = paste('count =', 
  #                     format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE), 
  #                     '\n')
  #     )
  #   )
  # }
  # fig <-fig + 
  #   stat_summary(
  #     fun.data = stat_box_data, 
  #     geom = "text", 
  #     hjust = 0.5,
  #     vjust = 0.9, size=3
  #   ) 
  fig
  #annotate_figure(fig3, fig.lab = "GG: domesticus, AA: musculus", fig.lab.pos = "bottom.right", fig.lab.size = 12)
  ggsave(file=paste0(taxa,"_marker_", marker,"_",DNAorRNA, "_box.pdf"), plot=fig)
  
  fig2<-ggplot(tot_with_SV, aes(x=marker, y=log_counts)) + geom_boxplot()+ 
   # stat_compare_means(comparisons = my_comp, label = "p.signif")  
  theme(legend.position = "none") +  geom_jitter( aes(color=marker))+
    labs(x="Genotypes",y="Log10 transformed absolute abundance")+ scale_color_manual(values = cbbPalette)+
    scale_x_discrete(limits=c(a, "H", b),labels=c(paste0(a,a), paste0(a,b), paste0(b,b)))+
    ggtitle(paste0(taxa," abundance for marker ", marker)) + labs(caption = paste0("Mmm: ", mus_dom[,1], ", Mmd: ", mus_dom[,2]))
  ggsave(file=paste0(taxa,"_marker_", marker,"_",DNAorRNA, "log_box.pdf"), plot=fig2)
  return(fig)
  
}
effect.violin.plot <- function(taxa, trait, DNAorRNA, marker){
  require(readxl)
  require(ggplot2)
  require(reshape2)
  require(cowplot)
  require(ggpubr)
  require(gtools)
  
  cbbPalette <- c("#1f78b4",  "#a6cee3","#b2df8a","#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", '#6a3d9a', "#ffff99", "#b15928", "#40e0d0", "#ffff00", "#ee82ee", "#fb8072", "#191970", "#fdb462", "#CBD4D4", "#fccde5", "#d9d9d9", "#bc80bd", "#32cd32", "#ffed6f", "lightcoral" )
  # load("~/Documents/PhD/Experiments/QTL_mapping_results/cleaning/clean_geno.Rdata")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_geno.Rdata")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Plots/Effect_plots/ref_alt_min_major.Rdata")
  theme_set(theme_test())
  
  clean_geno <- as.data.frame(clean_geno)
  colnames(clean_geno) <- gsub(x = colnames(clean_geno), pattern = "\\/", replacement = ".")  
  all_mark <-clean_geno[marker,]
  info_marker <- genotypes_ref[marker,]
  a <- info_marker$AA_min # a1 is the reference allele BB 1
  b <- info_marker$BB_maj # a2 is the minor allele AA -1 
  genotypes <- t(all_mark)
  
  consensus <-read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/consensus_mus_dom.csv", sep=",")
  allele_freq <-read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/allele_freq_consensus_mus_dom.csv", sep=",")
  rownames(consensus)<- consensus$X 
  consensus$X <- NULL
  mus_dom <- consensus[marker,]
  
  
  # abundances
  #DNAorRNA <- "DNA"
  # otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Miseq_hybrid_mice/DNA_and_RNA/tables_core/",trait, "_table_f2_core.csv"))
  otu_table <- read.csv(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/",trait, "_table_f2_core.csv"))
  rownames(otu_table)<- otu_table$X
  otu_table <- otu_table[which(grepl(paste0("_",DNAorRNA),rownames(otu_table))),] # only do the DNA
  rownames(otu_table) <- substr(rownames(otu_table), 1,15)
  
  SV <- as.data.frame(otu_table[,taxa])
  rownames(SV) <- rownames(otu_table)
  tot_with_SV <- merge(genotypes, SV, by="row.names")
  colnames(tot_with_SV)<- c("individuals", "marker", "abundance")
  tot_with_SV$counts <- tot_with_SV$abundance*10000
  tot_with_SV$log_counts <- log10(tot_with_SV$counts+1)
  tot_with_SV$marker<- as.factor(tot_with_SV$marker)
  tot_with_SV$transformed <- inv.logit(tot_with_SV$abundance)
  
  my_comp <- list(c(a, b),c(a, "H"), c(b, "H"))
  
  viol_fig <-  ggplot(tot_with_SV, aes(x=marker, y=abundance)) + geom_violin(aes(fill=marker))+
    # geom_boxplot(aes(fill=marker), width=0.1)+ 
    theme(legend.position = "none") +
    labs(x="Genotypes",y="Relative abundance")+ scale_fill_manual(values = cbbPalette)+
    scale_x_discrete(limits=c(a, "H", b),labels=c(paste0(a,a), paste0(a,b), paste0(b,b)))+
    ggtitle(paste0(taxa," abundance for marker ", marker)) + scale_y_continuous(labels = scales::percent)+ 
    labs(caption = paste0("Mmm: ", mus_dom[,1], ", Mmd: ", mus_dom[,2]))
  
  # Function to produce summary statistics (mean and +/- sd)
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    min<- min(x)
    if (min>ymin){ymin=min}
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  viol_fig <-viol_fig + stat_summary(fun.data=data_summary) 

  stat_box_data <- function(x, upper_limit = max(tot_with_SV$abundance) * 1.25) {
    return(
      data.frame(
        y = upper_limit,
        label = paste('mean =',
                      format(mean(x), big.mark = ",", decimal.mark = ".", scientific = FALSE),
                      '\n', "count=", length(x))
      )
    )
  }
  viol_fig2<-viol_fig + 
    stat_summary(
      fun.data = stat_box_data,
      geom = "text",
      hjust = 0.5,
      vjust = 0.9, size=3
    )
  ggsave(file=paste0(taxa,"_marker_", marker,"_",DNAorRNA, "_violin.pdf"), plot=viol_fig)
  
  
  viol_fig_log <-  ggplot(tot_with_SV, aes(x=marker, y=log_counts)) + geom_violin(aes(fill=marker))+
    geom_boxplot(aes(fill=marker), width=0.1)+ 
    theme(legend.position = "none") +
    labs(x="Genotypes",y="Log10 transformed absolute abundance")+ scale_fill_manual(values = cbbPalette)+
    scale_x_discrete(limits=c(a, "H", b),labels=c(paste0(a,a), paste0(a,b), paste0(b,b)))+
    ggtitle(paste0(taxa," abundance for marker ", marker)) 
  ggsave(file=paste0(taxa,"_marker_", marker,"_",DNAorRNA, "_log_violin.pdf"), plot=viol_fig_log)
  
  return(viol_fig2)
}
