library(argyle)
library(ggplot2)
library(cowplot)
library("optparse")

set.seed(314)
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/gemma_test/")



# ONLY use this for testing -----
analysis <-"body_weight"
#DNAorRNA <-"DNA"
#covar <- "time"
source("~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/plotting_functions_gemma.r")
source("~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/function_for_gemma.r")
F2 <- read.plink("../../Cleaning_snps/clean_f2")
gemma <- "~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/bin/gemma"
load("../../Cleaning_snps/clean_snps.Rdata")
map <- snps
covar_file <- read.covar("../../Phenotypes_histology_new.csv")
taxa <- covar_file[1:327,]

rownames(taxa) <- taxa$mouse_name
# -------------------------------




cat("Initiating Analysis ",analysis,".\n",sep="")

# GENOTYPE DATA ####
cat("Reading in genotyping data.\n")
# READING INPUT FILE


F2 <- recode(F2, "relative")

# Map file ####
cat("Reading map file.\n")
#map <- read.table(paste0(inputdir,"genotypes/SNP_info_all.csv"), sep=",", header=TRUE) 
map$cM <- NULL
map$tier <- NULL
rownames(map)<-map$marker
map$pos <- map$pos * 10e5


# PHENOTYPE and COVARIATE data ####
#read in phenotypes and covariate file 
cat("Loading in phenotype data.\n")

# TAXA files
# taxa <- read.pheno(paste0(inputdir,"phenotypes/",analysis,"_transformed_f2_core.csv"))
# 
# taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),] # only do the DNA
# rownames(taxa) <- substr(rownames(taxa), 1,15) # remove the _DNA
# Change the names from / to . for the genotype file 
colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")  

# only where we have genotype info for
indi <- colnames(F2)
taxa <- taxa[indi,]
# select the individuals where we have phenotype data for
ind <- rownames(taxa)
covar_file <- covar_file[covar_file$id %in% ind,]

#choose what phenotype to analyse and what to use as covariate(s)
all_covars <- colnames(covar_file)
all_taxa <- colnames(taxa)




#write.gemma.covariates("covar.txt", covar, covar_file)

#Select the individuals from the genotype file 
F2 <- as.matrix(as.data.frame(F2))
individuals <- ncol(F2)
F2_part <- apply(F2[,7:individuals], 2, as.numeric)# make sure that the genotype information(0,1,2) is numeric

F2_part <- F2_part[,colnames(F2_part) %in% ind]

F2 <- cbind(F2[,1:6],F2_part)



# Running GEMMA ####

for (trait in all_taxa[c(15,16, 20:27)]) {
  trait <- "body_weight"
  covar <- "body_length"
  
  write.gemma.covariates("covar.txt", covar, covar_file)
  
  cat("Running GEMMA for ", trait,".\n", sep="")
  write.gemma.pheno("pheno.txt", trait, taxa)
  
  chromosomes <- c(1:19)
  scans        <- vector("list",length(chromosomes))
  names(scans) <- chromosomes
  pve          <- rep(NA,length(chromosomes))
  
  for (chr in chromosomes) {
    #compute the kinship matrix using all the markers that are not on the chromosome
    cat("Mapping QTLs on chromosome ", chr, ".\n",sep="")
    
    
    kinship_markers <- which(map$chr!=chr&map$chr!="X" &map$chr != "Y" & map$chr != "M")
    #write geno file
    write.gemma.geno("kinship_geno.txt",F2[kinship_markers,])
    write.gemma.map("kinship_map.txt", map[kinship_markers,])
    
    # Compute the kinship matrices LOCO
    system(paste(gemma, "-g kinship_geno.txt -a kinship_map.txt -gk 1 -p pheno.txt -outdir . -o kinship" ))
    
    #write geno file
    markers <- which(map$chr==chr)
    x <- paste0("chr",chr)
    geno <- subset(F2, F2[,"chr"]==x)
    write.gemma.geno("geno.txt",geno)
    write.gemma.map("map.txt", map[markers,])
    
    # run gemma
    system(paste(gemma, "-g geno.txt -a map.txt -p pheno.txt -k ../../Kinship_RNA/kinship_chr1.cXX.txt -c covar.txt -lmm 4 -outdir ./output" ))
    
    all_results <- merge(gemma_results, results, by.x="rs", by.y="marker")
    All_pvals <- all_results[,c(13,14,15,32,33, 34)]

    rcorr(as.matrix(All_pvals), type = "spearman")
# Load the results of the GEMMA association analysis.
    #scans[[chr]] <-read.gemma.assoc("output/result.assoc.txt")
    
    gemma_results <- read.table("output/result.assoc.txt", header = T)
    # Get the estimate of the proportion of variance explained by the
    # polygenic effects on all chromosomes except the current
    # chromosome.
    out      <- scan("output/result.log.txt",
                     what = "character",sep = " ",quiet = TRUE)
    pve[chr] <- out[which(out == "pve") + 7]
    
    
  }
  
  gwscan           <- do.call(rbind,scans)
  rownames(gwscan) <- do.call(c,lapply(scans,rownames))
  
  save(gwscan,file=paste0("output/gwscan_", trait, ".rdata" ))
  
  trait_table <- data.frame()
  for (chr in 1:19){
    tx <-  readRDS(paste0("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Sterility/", trait, "/", trait, "_chr_", chr, "with_add_dom.rds"))
    trait_table <- rbind(trait_table, tx)
  }
  tot <- merge(trait_table, gwscan, by.x="marker", by.y="row.names")
  pvalues <- tot[,c(1,18, 19,20,24)]
  rownames(pvalues)<-pvalues$marker
  pvalues$marker <- NULL
  colnames(pvalues)<-c("P_glm", "add_P", "dom_P", "P_gemma")
  
  M <-cor(pvalues,use = "pairwise.complete.obs", method = "spearman")
  pdf(paste0(trait, "_correlation_P.pdf"))
  corrplot(M, method = "number", type="upper") +title(trait)
  dev.off()
}
