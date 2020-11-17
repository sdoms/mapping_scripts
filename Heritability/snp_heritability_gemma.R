library(argyle)
library(ggplot2)
library(cowplot)
library(MASS)

set.seed(314)
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/snp_heritability/")
source("~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/plotting_functions_gemma.r")
source("~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/function_for_gemma.r")

aa <- c("DNA", "RNA")
all_analyses <- c("phylum", "class", "order", "family",  "genus", "otu")

#gemma <- "~/Documents/PhD/Experiments/QTL_mapping/Old_code/GEMMA/bin/gemma"
gemma <-"/usr/local/bin/gemma"

for (DNAorRNA in aa){
  DNAorRNA <-"DNA"
  #for (analysis in all_analyses){
# parameters -----
  analysis <-"otu"
    load("../../Cleaning_snps/clean_snps.Rdata")
    map <- snps
    covar_file <-  read_delim("~/Documents/PhD/Experiments/Final_QTL_mapping/Phenotypes_histology_new.csv", 
                              ";", escape_double = FALSE, trim_ws = TRUE)
    #rownames(covar_file) <- covar_file$mouse_name
cat("Initiating Analysis ",analysis,".\n",sep="")

# GENOTYPE DATA ####
cat("Reading in genotyping data.\n")
# READING INPUT FILE

    
F2 <- read.plink("../../Cleaning_snps/clean_f2")
F2 <- recode(F2, "relative")

# Map file ####
cat("Reading map file.\n")
#map <- read.table(paste0(inputdir,"genotypes/SNP_info_all.csv"), sep=",", header=TRUE) 
#map$cM <- NULL
map$tier <- NULL
#rownames(map)<-map$marker
map$pos <- map$pos * 10e5


# PHENOTYPE and COVARIATE data ####
#read in phenotypes and covariate file 
cat("Loading in phenotype data.\n")

# TAXA files
taxa <- read.pheno(paste0("../../Phenotyping_27.02.2020/tables_core/",analysis,"_table_f2_core.csv"))
# 
taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),] # only do the DNA
rownames(taxa) <- substr(rownames(taxa), 1,15) # remove the _DNA
# Change the names from / to . for the genotype file 
colnames(F2) <- gsub(x = colnames(F2), pattern = "\\/", replacement = ".")  

# only where we have genotype info for
indi <- colnames(F2)
taxa <- taxa[rownames(taxa) %in% indi,]
# select the individuals where we have phenotype data for
ind <- rownames(taxa)
covar_file <- covar_file[covar_file$Mouse_Name %in% ind,]

#choose what phenotype to analyse and what to use as covariate(s)
all_covars <- colnames(covar_file)
all_taxa <- colnames(taxa)




#write.gemma.covariates("covar.txt", covar, covar_file)

#Select the individuals from the genotype file 
F2 <- as.matrix(as.data.frame(F2))
individuals <- ncol(F2)
F2_part <- apply(F2[,7:individuals], 2, as.numeric)# make sure that the genotype information(0,1,2) is numeric

F2_part <- F2_part[,colnames(F2_part) %in% ind]

F2 <- cbind(map[,1:6],F2_part)

write.gemma.geno("geno.txt",F2)
write.gemma.map("map.txt", map)

# # compute kinship matrix

# system(paste(gemma, "-g geno.txt -a map.txt -gk 1 -p pheno.txt -outdir . -o kinship_RNA" ))

all_pve <- data.frame()

if (analysis=="otu"&DNAorRNA=="RNA"){
  all_taxa<- all_taxa[-8]
}

scans        <- vector("list",length(all_taxa))
names(scans) <- all_taxa

# running GEMMA
for (trait in all_taxa) {

  pheno_file <- cbind(covar_file, taxa[trait] )
  colnames(pheno_file)<- c(colnames(covar_file),"tax")
  # regress matingpair on phenotypes and use residuals as new phenotypes
  dd <- glm.nb(tax ~ mating.pair, data=pheno_file)
  y <- round(residuals(dd),digits = 6)
  write.table(y,"pheno.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
  # run gemma
  system(paste0(gemma, " -g geno.txt -a map.txt -p pheno.txt -k kinship_",DNAorRNA, ".cXX.txt -lmm 4 -outdir ./output " ))
  out      <- scan("output/result.log.txt",
                  what = "character",sep = " ",quiet = TRUE)
  pve<- out[which(out == "pve") + 7]
  se <- out[which(out == "se(pve)") + 6]

  all_pve[trait,1]<- pve
  all_pve[trait,2]<- se
  
  scans[[trait]] <-read.gemma.assoc("output/result.assoc.txt")
}
colnames(all_pve)<- analysis
write.table(all_pve, file=paste0("pve_", DNAorRNA, "_",analysis, ".txt" ))


save(scans,file=paste0("output/gwscan_", analysis, ".rdata" ))


#}
}
