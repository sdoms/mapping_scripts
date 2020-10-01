writeKinship <- function (plink, pheno,outputdir, individuals){
  require(argyle)
  gemma <- "~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/bin/gemma"
  #geno <- read.plink("../../QTL_mapping_results/cleaning/clean_f2")
  geno <- read.plink(plink)
  # Map file ####
  cat("Reading map file.\n")
  load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata")
  snps$cM <- NULL
  snps$tier <- NULL
  snps$pos <- snps$pos * 10e5
  
  #Change the names from / to . for the genotype file 
  colnames(geno) <- gsub(x = colnames(geno), pattern = "\\/", replacement = ".")  

  # pheno table 
  write.table(pheno[individuals,2],"pheno.txt", quote=F, row.names=F, col.names=F)
  # relatedness matrix ####
  for (chr in c(1:19, "X"){
    cat("Calculating kinship matrix for chromosome ", chr, "\n")
    kinship_markers <- snps[which(snps$chr!=chr&snps$chr!="X" &snps$chr != "Y" & snps$chr != "M"),1]
    # write geno
    gts <- as.data.frame(geno, check.names=F)
    gts <- gts[kinship_markers,c("marker", "A1", "A2", individuals)]

    write.table(gts,"kinship_geno.txt",sep = " ",quote = FALSE,row.names = FALSE,
                col.names = FALSE)
    # write snps
    write.table(snps[kinship_markers,c("marker","pos","chr")],"kinship_map.txt",sep = " ",quote = FALSE,
                row.names = FALSE,col.names = FALSE)
    
    
    # Compute the kinship matrices LOCO
    system(paste0(gemma, " -g kinship_geno.txt -a kinship_map.txt -gk 1 -p pheno.txt -debug -outdir ", outputdir," -o kinship_chr",chr ))
    
  }
  
}
DNAorRNA <- "RNA"
taxa <- read.csv("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotyping_27.02.2020/tables_core/otu_table_f2_core.csv", header=T, row.names = 1)
taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),]
names_taxa<- substr(rownames(taxa), 1,nchar(rownames(taxa))-4)
# Logit Transform 
library(gtools)
taxa <- as.data.frame(sapply(taxa, function(x) inv.logit(x), USE.NAMES = T))
rownames(taxa)<- names_taxa
pheno<- taxa
individuals <- names_taxa
plink <- "~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_f2"
outputdir <- paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Kinship_", DNAorRNA, "/")
writeKinship("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_f2", taxa,paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Kinship_", DNAorRNA, "/"), names_taxa)









