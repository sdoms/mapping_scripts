setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/")
library(argyle)
library(qtl2convert)



load("../../QTL_mapping_results/Genotype_data/GM_snps.Rdata")

musdom <- read.csv("~/Documents/PhD/Mice/PayseurLab_GigaMUGA_data/BretPayseurGigaMUGA_May2015.csv", sep=",")
sample_key <- read.csv("~/Documents/PhD/Mice/PayseurLab_GigaMUGA_data/SampleKey.csv")


# mismatch between marker names in geno and snps files
# replace "." with "-" in rownames and marker names
# these are markers like "B6_01-074963079-S" which got corrupted to "B6_01.074963079.S"
rownames(GM_snps) <- gsub("\\.", "-", rownames(GM_snps))
GM_snps[,1] <- gsub("\\.", "-", GM_snps[,1])
musdom$Marker <-gsub("\\.", "-", musdom$Marker)
rownames(musdom) <- musdom$Marker

# getting only mus and dom separately
mus <- sample_key[which(sample_key$Subspp=="mus"), "ColnameOrigGenoFile" ]
mus <- gsub("[\\(\\)\\-]","\\.", mus)
dom <- sample_key[which(sample_key$Subspp=="dom"), "ColnameOrigGenoFile" ]
dom <- gsub("[\\(\\)\\-]","\\.", dom)
geno_mus <- musdom[,mus]
geno_dom <- musdom[,dom]


#finding the consensus sequence along mus and dom
consensus <-matrix(ncol=2, nrow=nrow(musdom))
colnames(consensus) <- c("mus", "dom")
rownames(consensus)<- rownames(musdom)
consensus[,1]<- qtl2convert::find_consensus_geno(as.matrix(geno_mus), na.strings="N")
consensus[,2]<- qtl2convert::find_consensus_geno(as.matrix(geno_dom), na.strings="N")

# replace the weird numeric values with NA
#library(SOfun)
#cons <-makemeNA(as.data.frame(consensus), "[ACGTH]", FALSE)
write.csv(consensus, file="consensus_mus_dom.csv")


allele_freq <- data.frame(consensus)
allele_freq$mm_freq <- "NA"
allele_freq$dd_frew <- "NA"
for (marker in rownames(consensus)){
  mm <- consensus[marker,1]
  dd <- consensus[marker,2]
  if (is.na(mm)){
    mm <- "N"
  }
  if (is.na(dd)){
    dd <- "N"
  }
  mm_tab <- as.data.frame(table(unlist(geno_mus[marker,])))
  dd_tab <- as.data.frame(table(unlist(geno_dom[marker,])))
  
  mm_freq <- mm_tab[mm_tab$Var1==mm,][[2]]/20
  dd_freq <- dd_tab[mm_tab$Var1==dd,][[2]]/20
  allele_freq[marker,"mm_freq"] <- mm_freq
  allele_freq[marker, "dd_freq"] <- dd_freq
}
write.csv(allele_freq, file="allele_freq_consensus_mus_dom.csv")

library(matrixStats)
counts_dom <- data.frame(consensus)
counts_dom$A <-rowCounts(as.matrix(geno_dom), value = "A")
counts_dom$C <-rowCounts(as.matrix(geno_dom), value = "C")
counts_dom$H <-rowCounts(as.matrix(geno_dom), value = "H")
counts_dom$T <-rowCounts(as.matrix(geno_dom), value = "T")
counts_dom$G <-rowCounts(as.matrix(geno_dom), value = "G")
counts_dom$N <-rowCounts(as.matrix(geno_dom), value = "N")
write.csv(counts_dom, file="counts_dom.csv")


counts_mus <- data.frame(consensus)
counts_mus$A <-rowCounts(as.matrix(geno_mus), value = "A")
counts_mus$C <-rowCounts(as.matrix(geno_mus), value = "C")
counts_mus$H <-rowCounts(as.matrix(geno_mus), value = "H")
counts_mus$T <-rowCounts(as.matrix(geno_mus), value = "T")
counts_mus$G <-rowCounts(as.matrix(geno_mus), value = "G")
counts_mus$N <-rowCounts(as.matrix(geno_mus), value = "N")

write.csv(counts_mus, file="counts_mus.csv")
