setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/otu/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_csv2(all_files, file="../otu_allresults.csv")
save(all_files, file="../otu_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/genus/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_csv2(all_files, file="../genus_allresults.csv")
save(all_files, file="../genus_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/")
list_of_files <- list.files(".",pattern="_allresults.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv2(file = fi)
  infile$tax_level <- gsub("_allresults.csv", "", fi)
  all_files <- rbind(all_files, infile)
}
all_files <- all_files %>% 
  drop_na(trait)
save(all_files,file="all_DNA.Rdata")
write_csv2(all_files, file="all_sig_DNA.csv")

med<- median(all_files$length[all_files$length>0], na.rm=T)


library(intervals)
reduced_intervals<- data.frame(matrix(ncol=3))
colnames(reduced_intervals)<- c("chr", "start", "end")
for (chr in 1:19){
  x<- Intervals(as.matrix(all_files[all_files$chr==chr,8:9]))
  y<- reduce(x)
  z<- data.frame(y)
  colnames(z)<- c("start", "end")
  z$chr<- paste0("chr",chr)
  reduced_intervals<- rbind(reduced_intervals,z)
}


setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/otu/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv2(file = fi)
  all_files <- rbind(all_files, infile)
}

write_csv2(all_files, file="../otu_allresults.csv")
save(all_files, file="../otu_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/phylum/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv2(file = fi)
  all_files <- rbind(all_files, infile)
}

write_csv2(all_files, file="../genus_allresults.csv")
save(all_files, file="../genus_allresults.Rdata")
write_csv2(all_files, file="../family_allresults.csv")
save(all_files, file="../family_allresults.Rdata")


write_csv2(all_files, file="../order_allresults.csv")
save(all_files, file="../order_allresults.Rdata")
write_csv2(all_files, file="../class_allresults.csv")
save(all_files, file="../class_allresults.Rdata")
write_csv2(all_files, file="../phylum_allresults.csv")
save(all_files, file="../phylum_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/")
list_of_files <- list.files(".",pattern="_allresults.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv2(file = fi)
  infile$tax_level <- gsub("_allresults.csv", "", fi)
  all_files <- rbind(all_files, infile)
}
all_files <- all_files %>% 
  drop_na(trait)
save(all_files,file="all_RNA.Rdata")
write_csv2(all_files, file="all_sig_RNA.csv")

med<- median(all_files$length[all_files$length>0], na.rm=T)


library(intervals)
reduced_intervals<- data.frame(matrix(ncol=3))
colnames(reduced_intervals)<- c("chr", "start", "end")
for (chr in 1:19){
  x<- Intervals(as.matrix(all_files[all_files$chr==chr,8:9]))
  y<- reduce(x)
  z<- data.frame(y)
  colnames(z)<- c("start", "end")
  z$chr<- paste0("chr",chr)
  reduced_intervals<- rbind(reduced_intervals,z)
}
save(reduced_intervals,file="reduced_intervals_RNA.rdata")
