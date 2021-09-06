Genetic mapping of host loci influencing gut microbiome composition in
hybrid mice
================

This github repository contains the original code from the article
“Genetic mapping of host loci influencing gut microbiome composition in
hybrid mice” published in …

# Data

The data folder contains the genotypes and the phenotypes (bacterial
abundance tables). The genotypes folder contains the raw genotypes
(all\_genotypes\_all\_samples.Rdata), the pre-LD-filtered SNPs
“hybrid\_f2\_5”, and the cleaned genotypes “clean\_f2”.

# Cleaning

Script 1\_cleaning\_snps.R uses Argyle R package and PLINK to clean the
genotype data. Script 2\_consensus\_F0.R calculates the consensus in the
G0 generation.

# Kinship

Using GEMMA to calculate kinship matrices.

# Heritability

Calculation of the heritability estimates and the correlation with the
cospeciation rates by Groussin et al. 2017.

# Association mapping

Folder ’Run\_script" contains separate scripts for the autosomes and the
X chromosome.

# Plotting

## Manhattan plots

Script run\_manhattan.R will make the Manhattan plots using the
pretty\_manhattan.R script for each trait. Script
run\_overlay\_manhattan.r takes the minimum P values of all the traits
within one taxonomic level.

## Region Plots

To make region plots run the run\_pretty\_region\_plot.R script. It
calls the New\_pretty\_plot\_script.R. New\_pretty\_plot.R first
calculates the confidence interval calculating the pairwise LD to the
peak SNP and takes the positions from the two outer SNPs, that are in
LD&gt;0.9 with the peak SNP, as the borders. It then uses biomart to
find the genes in the interval and plot this in a locuszoom-like manner.
At the end, it will write a result file for the trait with the
information about the significant regions.

## Effect plots

effect\_plot.R takes a taxa, a taxonomic level, DNA or RNA and a marker
name as input and will plot a boxplot of the abundance for each genotype
and also prints the consensus genotype for musculus and domesticus mice.

effect\_plot\_error\_bar.R plots it with an error bar instead of a box
plot.

## Genotyping Plots

This contains scripts for a PCA, MAF distribution and to visualize the
location of the SNPs on the chromomes.

## Overview Plots

This folder contains scripts to plot the results on a karyotype like
plot to visualize the density of genomic regions associated with
bacterial traits.

# Percentage of variance explained (PVE)

PVE\_SNP\_model.R shows how to calculate to PVE of each SNP using a
model and PVE\_snp.R using a formula.

# matSpDLite

This script uses matSpDLite
(<http://gump.qimr.edu.au/general/daleN/matSpDlite/>) from Nyholt DR
(2004) to calculate the number of independent taxa tested in order to
calculate a study-wide significance threshold.

# Summary Tables

gwascanSummaries.R combines SNPs within 10Mb from each other into single
regions and expands those regions with SNPs in LD (&gt;0.9). It then
uses biomart and mm10 mouse genome to determine the genes within those
intervals.

# Enrichment

Enrichment analysis of the genes closest to the significant SNPs using
the clusterProfiler R package.

# Genes

Calculate the overlap of genes/regions with other studies using
poverlap.
