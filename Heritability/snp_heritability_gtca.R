setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/")
library(gtools)
library(filesstrings)
library(MASS)

# does not converge !!!
DNAorRNA <-"RNA"
trait <- "otu"

dir.create(paste0("../snp_heritability/",trait, "_", DNAorRNA,"/"), showWarnings = F)

taxa <- read.csv(paste0("~/Documents/PhD/Experiments/Miseq_hybrid_mice/DNA_and_RNA/tables_core/", trait, "_table_f2_core.csv"), row.names = "X")
taxa <- taxa[which(grepl(paste0("_",DNAorRNA),rownames(taxa))),]
rownames(taxa)<- substr(rownames(taxa), 1,nchar(rownames(taxa))-4)


taxa <- as.data.frame(sapply(taxa, function(x) inv.logit(x), USE.NAMES = T),row.names=rownames(taxa))
ids <- rownames(taxa)
ids2 <- gsub(pattern = "\\.H", replacement = "/H", x = ids)
famId <- substr(ids2, 1,2)

# grm file
#system(paste0("/Users/doms/Software/gcta_1.93.0beta_mac/bin/gcta64 --bfile clean_f2 --make-grm --out gtcagrm"))
#system((paste0("/Users/doms/Software/gcta_1.93.0beta_mac/bin/gcta64 --grm gtcagrm --make-bK 0.05 --out grm2_bK")))
# write the pheno files
for (tax in colnames(taxa)){
  pheno <- data.frame(famId, ids2, taxa[,tax])
  colnames(pheno)<- c("mating.pair", "mouseID", "tax")
  dd <- glm.nb(tax ~ mating.pair, data=pheno)
  phe_out<- pheno[,1:2]
  phe_out$resids <- residuals(dd)
  write.table(phe_out, paste0("../snp_heritability/", trait,"_",DNAorRNA,"/",tax,"_", DNAorRNA, ".pheno"), col.names = FALSE, row.names = FALSE, quote=F)
  try(system(paste0("/Users/doms/Software/gcta_1.93.0beta_mac/bin/gcta64 --reml --mgrm multi_grm.txt --reml-maxit 500 --pheno ../snp_heritability/",trait,"_",DNAorRNA,"/",tax,"_", DNAorRNA, ".pheno --out ",tax,"_bKsK")))
  try(file.move(paste0(tax,"_bKsK.hsq"), paste0("../../QTL_mapping/snp_heritability/",trait, "_", DNAorRNA,"/" )))
  try(file.move(paste0(tax,"_bKsK.log"), paste0("../../QTL_mapping/snp_heritability/",trait, "_", DNAorRNA,"/" )))
  try(system(paste0("/Users/doms/Software/gcta_1.93.0beta_mac/bin/gcta64 --reml --mgrm add_domi_grm.txt --pheno ../../QTL_mapping/snp_heritability/",trait,"_",DNAorRNA,"/",tax,"_", DNAorRNA, ".pheno --thread-num 10 --reml-maxit 250 --out ",tax,"_add_domi")))
  try(file.move(paste0(tax,"_add_domi.hsq"), paste0("../../QTL_mapping/snp_heritability/",trait, "_", DNAorRNA,"/" )))
  try(file.move(paste0(tax,"_add_domi.log"), paste0("../../QTL_mapping/snp_heritability/",trait, "_", DNAorRNA,"/" )))
  
}

