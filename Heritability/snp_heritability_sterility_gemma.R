library(argyle)
library(ggplot2)
library(cowplot)
library(MASS)

set.seed(314)
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/snp_heritability/")
source("~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/plotting_functions_gemma.r")
source("~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/function_for_gemma.r")



gemma <- "~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/bin/gemma"

    # parameters -----
    #analysis <-"otu"
    load("../../Cleaning_snps/clean_snps.Rdata")
    map <- snps
    covar_file <- read.covar("../../Phenotypes_histology_new.csv")
    cat("Initiating Analysis.\n")
    
    # GENOTYPE DATA ####
    cat("Reading in genotyping data.\n")
    # READING INPUT FILE
    
    
    F2 <- read.plink("../../Cleaning_snps/clean_f2")
    F2 <- recode(F2, "relative")
    
    # Map file ####
    cat("Reading map file.\n")
    #map <- read.table(paste0(inputdir,"genotypes/SNP_info_all.csv"), sep=",", header=TRUE) 
    map$cM <- NULL
    map$tier <- NULL
    #rownames(map)<-map$marker
    map$pos <- map$pos * 10e5
    
    
    # PHENOTYPE and COVARIATE data ####
    #read in phenotypes and covariate file 
    cat("Loading in phenotype data.\n")
    
    # # TAXA files
    # taxa <- read.pheno(paste0("../../Phenotyping_27.02.2020/tables_core/",analysis,"_table_f2_core.csv"))
    # # 
    # taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),] # only do the DNA
    # rownames(taxa) <- substr(rownames(taxa), 1,15) # remove the _DNA
    # # Change the names from / to . for the genotype file 
    colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")  
    
    # only where we have genotype info for
    indi <- colnames(F2)
    # taxa <- taxa[rownames(taxa) %in% indi,]
    # # select the individuals where we have phenotype data for
    # ind <- rownames(taxa)
    covar_file <- covar_file[covar_file$id %in% indi,]
    
    #choose what phenotype to analyse and what to use as covariate(s)
    all_covars <- colnames(covar_file)
    # all_taxa <- colnames(taxa)
    
    
    
    
    #write.gemma.covariates("covar.txt", covar, covar_file)
    
    #Select the individuals from the genotype file 
    F2 <- as.matrix(as.data.frame(F2))
    individuals <- ncol(F2)
    F2_part <- apply(F2[,7:individuals], 2, as.numeric)# make sure that the genotype information(0,1,2) is numeric
    
    F2_part <- F2_part[,colnames(F2_part) %in% indi]
    
    F2 <- cbind(F2[,1:6],F2_part)
    
    # # compute kinship matrix
    # write.gemma.geno("kinship_geno.txt",F2)
    # write.gemma.map("kinship_map.txt", map)
    # system(paste(gemma, "-g kinship_geno.txt -a kinship_map.txt -gk 1 -p pheno.txt -outdir . -o kinship" ))
    
    
    all_pve <- data.frame()

    trait <-"combined_testis_weight"
    covar <- "body_weight"
    
    # running GEMMA
    #for (trait in all_covars) {
      pheno_file <- cbind(covar_file, covar_file[trait], covar_file[covar] )
      colnames(pheno_file)<- c(colnames(covar_file),"trait", "covar")
      # regress matingpair on phenotypes and use residuals as new phenotypes
      dd <- glm.nb(trait ~ covar, data=pheno_file)
      y <- round(residuals(dd),digits = 6)
      write.table(y,"pheno.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      # run gemma
      system(paste(gemma, "-g geno.txt -a map.txt -p pheno.txt -k ../../Kinship_DNA/kinship.cXX.txt -lmm 4 -outdir ./output" ))
      out      <- scan("output/result.log.txt",
                       what = "character",sep = " ",quiet = TRUE)
      pve<- out[which(out == "pve") + 7]
      
      all_pve[trait,1]<- pve
   # }
      
    colnames(all_pve)<- analysis
    write.table(all_pve, file=paste0("pve_", DNAorRNA, "_",analysis, ".txt" ))
    
