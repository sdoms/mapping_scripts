setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/otu/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../otu_allresults.csv", delim=";")
save(all_files, file="../otu_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/genus/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../genus_allresults.csv", delim=";")
save(all_files, file="../genus_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/family/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../family_allresults.csv", delim=";")
save(all_files, file="../family_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/order/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../order_allresults.csv", delim=";")
save(all_files, file="../order_allresults.Rdata")


setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/class/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../class_allresults.csv", delim=";")
save(all_files, file="../class_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/phylum/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../phylum_allresults.csv", delim=";")
save(all_files, file="../phylum_allresults.Rdata")
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/")
list_of_files <- list.files(".",pattern="_allresults.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_delim(file = fi, delim=";")
  infile$tax_level <- gsub("_allresults.csv", "", fi)
  all_files <- rbind(all_files, infile)
}
all_files <- all_files %>% 
  drop_na(trait)
save(all_files,file="all_DNA.Rdata")
write_delim(all_files, file="all_sig_DNA.csv", delim=";")

med<- median(all_files$length[all_files$length>0], na.rm=T)


library(intervals)
reduced_intervals<- data.frame(matrix(ncol=3))
colnames(reduced_intervals)<- c("chr", "start", "end")

for (chr in 1:19){
  x<- Intervals(as.matrix(all_files[all_files$chr==chr,9:10]))
  y<- reduce(x)
  z<- data.frame(y)
  colnames(z)<- c("start", "end")
  z$chr<- paste0("chr",chr)
  reduced_intervals<- rbind(reduced_intervals,z)
}

save(reduced_intervals,file="reduced_intervals_DNA.rdata")


############################################################################################3
###### RNA #####

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/otu/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../otu_allresults.csv", delim=";")
save(all_files, file="../otu_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/phylum/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../phylum_allresults.csv", delim=";")
save(all_files, file="../phylum_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/class/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../class_allresults.csv", delim=";")
save(all_files, file="../class_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/genus/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../genus_allresults.csv", delim=";")
save(all_files, file="../genus_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/family/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}



write_delim(all_files, file="../family_allresults.csv", delim=";")
save(all_files, file="../family_allresults.Rdata")

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/order/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}


write_delim(all_files, file="../order_allresults.csv", delim=";")
save(all_files, file="../order_allresults.Rdata")


setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries/")
list_of_files <- list.files(".",pattern="_allresults.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_delim(file = fi, delim=";")
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
  set<- all_files[all_files$chr==chr,9:10] %>% arrange(stop.LD.pos) %>% arrange(start.LD.pos)
  x<- Intervals(as.matrix(set))
  y<- reduce(x)
  z<- data.frame(y)
  colnames(z)<- c("start", "end")
  z$chr<- paste0("chr",chr)
  reduced_intervals<- rbind(reduced_intervals,z)
}
save(reduced_intervals,file="reduced_intervals_RNA.rdata")
