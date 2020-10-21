
###############################################################################################
##                               Calculate consensus genotypes                               ##
##  This script will calculate the consensus genotype for each of the ancestral mouse lines  ##
###############################################################################################


library(argyle)
library(qtl2convert)
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/")



##----------------------------------------------------------------
##                       For cleaned SNPs                       --
##----------------------------------------------------------------


load("clean_snps.Rdata")
load("parents_all_snps_geno.RData")
parents <- as.matrix(as.data.frame(parents))

# mismatch between marker names in geno and snps files
# replace "." with "-" in rownames and marker names
# these are markers like "B6_01-074963079-S" which got corrupted to "B6_01.074963079.S"
rownames(snps) <- gsub("\\.", "-", rownames(snps))
snps[,1] <- gsub("\\.", "-", snps[,1])
rownames(parents) <- gsub("\\.", "-", rownames(parents))

markers <- snps[,1]

clean_parents <- merge(snps, parents, by="marker")
clean_parents <- clean_parents[,c(1,17:47)]
colnames(clean_parents) <- c("marker", "chr", "pos", "cM", "A1", "A2", colnames(clean_parents)[7:32])
names <-colnames(clean_parents)[7:32]
codes <- substr(names, 1,4)

geno <- clean_parents 
rm(clean_parents)
rownames(geno)<- geno$marker


# consensus genotypes
reduced_geno <- matrix(ncol=8, nrow=nrow(geno))
let <- LETTERS[1:8]
dimnames(reduced_geno) <- list(rownames(geno), unique(codes))
codes_count <- table(codes)
for(i in unique(codes))
  reduced_geno[,i] <- qtl2convert::find_consensus_geno(geno[,substr(colnames(geno),1,4)==i])


save(reduced_geno, file="consensus_F0_clean.Rdata")
write.csv(reduced_geno, file="consensus_F0_clean.csv")


##----------------------------------------------------------------
##                         For all SNPs                         --
##----------------------------------------------------------------

geno <- as.data.frame(parents)
#rm(parents)
rownames(geno)<- geno$marker


# consensus genotypes
reduced_geno <- matrix(ncol=8, nrow=nrow(geno))
let <- LETTERS[1:8]
dimnames(reduced_geno) <- list(rownames(geno), unique(codes))
codes_count <- table(codes)
for(i in unique(codes))
  reduced_geno[,i] <- qtl2convert::find_consensus_geno(geno[,substr(colnames(geno),1,4)==i])


save(reduced_geno, file="consensus_F0_all.Rdata")
write.csv(reduced_geno, file="consensus_F0_all.csv")
