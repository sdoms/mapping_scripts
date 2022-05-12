## Code snippets from https://doi.org/10.1371/journal.pgen.1005607 
## Mapping of Craniofacial Traits in Outbred Mice Identifies Major Developmental Genes Involved in Shape Determination
## Luisa F. Pallares,Peter Carbonetto,Shyam Gopalakrishnan,Clarissa C. Parker,Cheryl L. Ackert-Bicknell,Abraham A. Palmer,Diethard Tautz

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Write the phenotype data to a file in the format used by GEMMA. Each
# line of the file contains one phenotype observation.
write.gemma.pheno <- function (file, phenotype, pheno) {
  y <- pheno[,phenotype]
  if (is.numeric(y))
    y <- round(y,digits = 6)
  write.table(y,file,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

# ----------------------------------------------------------------------
# Write the covariate data to a file in the format used by GEMMA. Each
# line corresponds to a sample. We must include an additional
# covariate for the intercept.
write.gemma.covariates <- function (file, covariates, pheno) {
  if (is.null(covariates)) {
    write.table(data.frame(rep(1,nrow(pheno))),file,sep = " ",
                quote = FALSE,row.names = FALSE,col.names = FALSE)
  } else {
    write.table(cbind(1,data.frame(lapply(subset(pheno, select=covariates),function (x) {
      if (is.numeric(x))
        round(x,digits=6)
      else
        x
    }))),file,sep = " ",quote = FALSE,row.names = FALSE,col.names = FALSE)
  }
}

# ----------------------------------------------------------------------
# Write the SNP information to a space-delimited text file in the
# format used by GEMMA. This file contains one line per SNP, with
# three columns: (1) SNP label, (2) base-pair position, (3)
# chromosome.
write.gemma.map <- function (file, map)
  write.table(map[c("marker","pos","chr")],file,sep = " ",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

# ----------------------------------------------------------------------
# Store the mean genotypes as a space-delimited text file in the
# format used by GEMMA, in which we have one row per SNP, and one
# column per sample. The first three column give the SNP label, and
# the two alleles.
write.gemma.geno <- function (file, geno) {
  geno <- as.data.frame(geno,check.names = FALSE)
  geno$pos <- NULL
  geno$cM <- NULL
  geno$chr <-NULL
  write.table(geno,file,sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
}

# ----------------------------------------------------------------------
# Reads in the covariates stored in a CSV file and returns a dataframe
read.covar <- function(file) {
  covar <- read.table(file, header=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep = ";")

  covar <- transform(covar,
                     id = as.character(id),
                     family= as.factor(family),
                     mating.pair = factor(mating.pair),
                     litter.id = factor(litter.id)
                     )
  return(covar)
}
# ----------------------------------------------------------------------
# Reads in the phenotypes stored in a CSV file and returns a dataframe
read.pheno <- function(file) {
  pheno <- read.table(file, header=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep = ",", row.names = "")
  return(pheno)
}

# ----------------------------------------------------------------------
# Reads in the association results from GEMMA, and returns a data
# frame containing three columns: chromosome number ("chr"); base-pair
# position ("pos"); and the base-10 logarithm of the p-value ("log10p").
read.gemma.assoc <- function (file) {
  gwscan <- read.table(file,sep = "\t",header = TRUE,check.names = FALSE,
                       quote = "",stringsAsFactors = FALSE)
  rownames(gwscan) <- gwscan$rs
  gwscan           <- gwscan[c("chr","ps","p_lrt")]
  gwscan           <- transform(gwscan,p_lrt = -log10(p_lrt))
  colnames(gwscan) <- c("chr","pos","log10p")
  return(gwscan)
}

# ----------------------------------------------------------------------
# Remove outliers from the phenotype data for the given phenotype,
# optionally conditioning on covariates. If covariates are specified,
# outliers are determined according to the residual of the phenotype
# after regressing on the covariates.
remove.outliers <- function (pheno, phenotype, covariates, outliers,
                            verbose = TRUE) {

 # Specify the function for removing the outliers.
 is.outlier <- function (x) {
   y           <- outliers(x)
   y[is.na(y)] <- FALSE
   return(y)
 }

 # If we are conditioning on one or more covariates, get the
 # residuals of the phenotype conditioned on these covariates.
 if (length(covariates) > 0) {
   f <- formula(paste(phenotype,"~",paste(covariates,collapse="+")))
   r <- resid(lm(f,pheno,na.action = na.exclude))
 } else
   r <- pheno[[phenotype]]

 # If requested, report how many outlying data points are removed.
 if (verbose) {
   n <- sum(is.outlier(r))
   if (n == 0)
     cat("No outliers for",phenotype)
   else
     cat(n,"outliers for",phenotype)
   if (length(covariates) == 0)
     cat(" are removed.\n")
   else
     cat(" conditioned on",paste(covariates,collapse=" + "),"are removed.\n")
 }

 # Remove the outlying data points.
 pheno[is.outlier(r),phenotype] <- NA

 # Return the updated phenotype data table.
 return(pheno)
}
# Centers the columns of matrix X so that the entries in each column
# of X add up to zero.
center.columns <- function (X) {
  mu <- matrix(colMeans(X),1,ncol(X))
  X  <- X - repmat(mu,nrow(X),1)
  return(X)
}

# ----------------------------------------------------------------------
# Does the same thing as repmat(A,m,n) in MATLAB.
repmat <- function (A,m,n)
  return(kronecker(matrix(1,m,n),A))
# ----------------------------------------------------------------------
# Centers the rows of matrix X so that the entries in each row
# of X add up to zero.
center.rows <- function (X) {
  mu <- matrix(rowMeans(X),1,nrow(X))
  X  <- X - t(repmat(mu,ncol(X),1))
  return(X)
}

