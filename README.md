Key features of the genetic architecture and evolution of host-microbe interactions revealed by high-resolution genetic mapping of the mucosa-associated gut microbiome in hybrid mice
================

This github repository contains the original code from the article
“Key features of the genetic architecture and evolution of host-microbe interactions revealed by high-resolution genetic mapping of the mucosa-associated gut microbiome in hybrid mice” (https://www.biorxiv.org/content/10.1101/2021.09.28.462095v3). This code repository aims to provide inspiration to conduct a similar GWAS analysis. Feel free to contact me when you are having questions by email (doms@evolbio.mpg.de) or through github.

# Data

The data folder contains the genotypes and the phenotypes (bacterial
abundance tables). The genotypes folder contains the raw genotypes
(all\_genotypes\_all\_samples.Rdata), the pre-LD-filtered SNPs
“hybrid\_f2\_5”, and the cleaned genotypes “clean\_f2”.

# 1. Cleaning
First, we start with a quality control of the genotype data.
Script *1\_cleaning\_snps.R* uses Argyle R package and PLINK to clean the
genotype data.
We remove:
* Non-biallelic SNPs
* Individuals and SNPs with high levels of missingness (>0.1)
* SNPs with low MAF frequency (<0.05)
* SNPs out of Hardy-Weinberg equilibrium (P<.00001)
* SNPs pairwise LD > 0.9 in a window size of 5 SNPs by SNP

Script *2\_consensus\_F0.R* calculates the consensus in the
G0 generation.


# 2. Heritability
* The script snp_heritability_lme4qtl.R uses GEMMA to calculate a kinship matrix and will then use lme4QTL to calculate the SNP-based heritability estimates and make the plot (Fig. 1 A-B).
* *function_for_gemma.r* contains code snippets to load/write data frames correctly for GEMMA. Code originates from Pallares et al (2015): https://doi.org/10.1371/journal.pgen.1005607
* *cospec_vs_herit.R* combines the SNP-based heritability estimates with the cospeciation rates calculates by Groussin et al (2017) to make Fig 1C-D.


# 3. Association mapping

Folder *"Run\_script"* contains separate scripts for the autosomes and the
X chromosome. This script takes the phenotypes and genotypes as input. The taxa abundances are first inverse.logit transformed.
GEMMA is used for calculating kinship matrices using the leave-one-chromosome-out approach.
Using the lme4QTl R package, we calculate a null model with the mating pair and the kinship matrix as random effects and a snp model that has the additive and the dominance effect of the SNP. We use an anova test to determine the significance.

# 4. Genomic inflation
Calculate the genomic inflation and 'correct' the p values: script *revisions/genomic_inflation.R*

# 5. Summary Tables
* *gwascanSummaries.R* and *gwascanSummaries_p_values_corrected.R* combines SNPs within 10Mb from each other into single
regions and expands those regions with SNPs in LD (>0.9). It then
uses biomart and mm10 mouse genome to determine the genes within those
intervals.
* *combine_files.R* is used to combine all the genome-wide significant regions from the individual taxa by taxonomic level.
* *combine_files_genomic_control.R* combines the results of p-value corrected when lambdaGC>1.05 with not corrected for lambdaGC<1.05
* *smallest_intervals.R* filters the smallest intervals out and annotates
* *summary_tables_all_snps_SW_22.R* makes overview tables by significant SNPs
* *summary_numbers.R* calculates the summary statistics (ie median interval size, number of significant loci, ...)


# 5. Enrichment analysis

Enrichment analysis of the genes closest to the significant SNPs using
the clusterProfiler R package.

# 6. Plotting

## Manhattan plots

* *Script run\_manhattan.R* will make the Manhattan plots using the
pretty\_manhattan.R script for each trait. Script
* *run\_overlay\_manhattan.r* takes the minimum P values of all the traits
within one taxonomic level.

## Region Plots

To make region plots run the run\_pretty\_region\_plot.R script. It
calls the New\_pretty\_plot\_script.R. New\_pretty\_plot.R first
calculates the confidence interval calculating the pairwise LD to the
peak SNP and takes the positions from the two outer SNPs, that are in
LD>0.9 with the peak SNP, as the borders. It then uses biomart to
find the genes in the interval and plot this in a locuszoom-like manner.
At the end, it will write a result file for the trait with the
information about the significant regions.

## Effect plots

*effect\_plot.R* takes a taxa, a taxonomic level, DNA or RNA and a marker
name as input and will plot a boxplot of the abundance for each genotype
and also prints the consensus genotype for musculus and domesticus mice.

*effect\_plot\_error\_bar.R* plots it with an error bar instead of a box
plot.

## Genotyping Plots

This contains scripts for a PCA, MAF distribution and to visualize the
location of the SNPs on the chromomes.

## Overview Plots

This folder contains scripts to plot the results on a karyotype like
plot to visualize the density of genomic regions associated with
bacterial traits.
* *Intervals_density_plot.R* will produce panel A of Figure 2.
* *density_results_snps.R* will produce panel B of Figure 2.

# 7. Other

## Percentage of variance explained (PVE)

* *PVE\_SNP\_model.R* shows how to calculate to PVE of each SNP using a
model and PVE\_snp.R using a formula.
* *PVE\_all\_snp\_RNA.R* and *PVE\_all\_snp\_DNA.R* are used to calculate a polygenic score equivalent in order to determine the percentage of variance explained by all the significant SNPs for one trait.

## matSpDLite

This script uses matSpDLite
(<http://gump.qimr.edu.au/general/daleN/matSpDlite/>) from Nyholt DR
(2004) to calculate the number of independent taxa tested in order to
calculate a study-wide significance threshold.

## Genes

Calculate the overlap of genes/regions with other studies using
poverlap.
