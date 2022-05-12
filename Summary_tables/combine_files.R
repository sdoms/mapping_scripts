# setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/otu/")
setwd("~/Documents/Research//Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries_pval_corr/otu/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

# write_delim(all_files, file="../otu_allresults.csv", delim=";")
# save(all_files, file="../otu_allresults.Rdata")

write_delim(all_files, file="../otu_allresults_pval_corr.csv", delim=";")
save(all_files, file="../otu_allresults_pval_corr.Rdata")

# setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/genus/")
setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries_pval_corr/genus/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../genus_allresults_pval_corr.csv", delim=";")
save(all_files, file="../genus_allresults_pval_corr.Rdata")

setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries_pval_corr/family/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../family_allresults_pval_corr.csv", delim=";")
save(all_files, file="../family_allresults_pval_corr.Rdata")

setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries_pval_corr/order/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../order_allresults_pval_corr.csv", delim=";")
save(all_files, file="../order_allresults_pval_corr.Rdata")


setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries_pval_corr/class/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../class_allresults_pval_corr.csv", delim=";")
save(all_files, file="../class_allresults_pval_corr.Rdata")

setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries_pval_corr/phylum/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../phylum_allresults.csv", delim=";")
save(all_files, file="../phylum_allresults.Rdata")
setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries_pval_corr/")
list_of_files <- list.files(".",pattern="_allresults_pval_corr.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_delim(file = fi, delim=";")
  infile$tax_level <- gsub("_allresults.csv", "", fi)
  all_files <- rbind(all_files, infile)
}
all_files <- all_files %>% 
  drop_na(trait)
save(all_files,file="all_DNA_pval_corr.Rdata")
write_delim(all_files, file="all_sig_DNA_pval_corr.csv", delim=";")

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

save(reduced_intervals,file="reduced_intervals_DNA_pval_corr.rdata")


############################################################################################3
###### RNA #####
setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries_pval_corr/otu/")

# setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/otu/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

# write_delim(all_files, file="../otu_allresults.csv", delim=";")
# save(all_files, file="../otu_allresults.Rdata")

write_delim(all_files, file="../otu_allresults_pval_corr.csv", delim=";")
save(all_files, file="../otu_allresults_pval_corr.Rdata")


setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries_pval_corr/phylum/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../phylum_allresults_pval_corr.csv", delim=";")
save(all_files, file="../phylum_allresults_pval_corr.Rdata")

setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries_pval_corr/class/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../class_allresults_pval_corr.csv", delim=";")
save(all_files, file="../class_allresults_pval_corr.Rdata")

setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries_pval_corr/genus/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}

write_delim(all_files, file="../genus_allresults_pval_corr.csv", delim=";")
save(all_files, file="../genus_allresults_pval_corr.Rdata")

setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries_pval_corr/family/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}



write_delim(all_files, file="../family_allresults_pval_corr.csv", delim=";")
save(all_files, file="../family_allresults_pval_corr.Rdata")

setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries_pval_corr/order/")
require(tidyverse)

list_of_files <- list.files(".",pattern="_gwscan_sigIntervals.annot_bon.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_csv(file = fi)
  all_files <- rbind(all_files, infile)
}


write_delim(all_files, file="../order_allresults_pval_corr.csv", delim=";")
save(all_files, file="../order_allresults_pval_corr.Rdata")


setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/RNA/sig.summaries_pval_corr/")
list_of_files <- list.files(".",pattern="_allresults_pval_corr.csv")

all_files <- tibble()
for (fi in list_of_files){
  infile<- read_delim(file = fi, delim=";")
  infile$tax_level <- gsub("_allresults_pval_corr.csv", "", fi)
  all_files <- rbind(all_files, infile)
}
all_files <- all_files %>% 
  drop_na(trait)
save(all_files,file="all_RNA.Rdata")
write_csv2(all_files, file="all_sig_RNA_pval_corr.csv")

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
save(reduced_intervals,file="reduced_intervals_RNA_pval_corr.rdata")
