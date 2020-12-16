setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/genes_overlap/")

library(tidyverse)
mm10 <- read.table("mouse.mm10.genome")
colnames(mm10)<- c("chr", "chr_len")

overlaps <- read_excel("intersect_out.xlsx", sheet="overlaps")
colnames(overlaps)<- c("chr", "start", "end", "width")
chr1 <- overlaps[overlaps$chr=="chr2",]


#goodChrOrder <- paste0("chr",c(1:19,"X"))
goodChrOrder <- c(1:19)

output <- data.frame()
for (chr in goodChrOrder){
  df.tmp <- overlaps[overlaps$chr==paste0("chr",chr),]
  dd<- c()
  for (x in 1:nrow(df.tmp)){
    dd<- c(dd,seq(df.tmp$start[x],df.tmp$end[x], 1e6))
  }
  res<-as.data.frame(table(cut_interval(dd,length=1e6)))
  
  #res <- as.data.frame(table(cut_width(df.tmp$start/1e6,1, boundary = 0)))

  res$start <- as.numeric(as.character(substr(res$Var1,2,regexpr(",",res$Var1)-1)))
  res$end <- res$start +1e6
  res$chr <- chr
  output <- rbind(output,res)
}

# output <- data.frame()
# for (chr in goodChrOrder){
#   df.tmp <- overlaps[overlaps$chr==paste0("chr",chr),]
#   res <- as.data.frame(table(cut_width(df.tmp$start/1e6,1, boundary = 0)))
#   res2<-as.data.frame(table(cut_width(df.tmp$end/1e6,1, boundary = 0)))
#   res$start <- as.numeric(as.character(substr(res$Var1,2,regexpr(",",res$Var1)-1)))
#   res$end <- res$start +1
#   res$chr <- chr
#   output <- rbind(output,res)
# }
mean_freq <- mean(output$Freq, na.rm=T)

output <- output %>% filter(Freq != 0)

chrSizes <- mm10[1:19,]
chrSizes$chr<- 1:19

output$start<- as.numeric(as.character(output$start))/1e6
output$end<- as.numeric(as.character(output$end))/1e6
library(viridis)
snp_plot <- ggplot()+
  geom_segment(data = chrSizes,
               aes(x = chr, xend = chr, y = 0, yend = chr_len/1e6),
               lineend = "round", color = "#DEE1E3", size = 5) + 
  
  geom_segment(data=output, 
               aes(x=chr, xend=chr, y=start, yend=end, color=Freq), size=5)+
  coord_flip() + scale_color_viridis_c(option="magma", direction = -1, limits=c(0,30)) + 
  labs(y="Position (Mb)", x= "Chromosome", color="Overlap count") +theme_test()
snp_plot + scale_x_discrete(limits=factor(goodChrOrder))
ggsave("overlapping_studies_intervals.pdf")
