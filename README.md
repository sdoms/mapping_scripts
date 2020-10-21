Association mapping of the mouse gut microbiome
================

# Plotting

## Manhattan plots

Script run\_manhattan(\_histo).R will make the Manhattan plots using the
pretty\_manhattan.R script for each trait. Script
run\_overlay\_manhattan.r takes the minimum P values of all the traits
within one taxonomic level.

## Region Plots

To make region plots run the run\_pretty\_region\_plot.R script. It
calls the New\_pretty\_plot\_script.R. New\_pretty\_plot.R first
calculates the confidence interval calculating the pairwise LD to the
peak SNP and takes the positions from the two outer SNPs, that are in
LD\>0.9 with the peak SNP, as the borders. It then uses biomart to find
the genes in the interval and plot this in a locuszoom-like manner. At
the end, it will write a result file for the trait with the information
about the significant regions.

## Effect plots

abundance\_effect\_plot.R takes a taxa, a taxonomic level, DNA or RNA
and a marker name as input and will plot a boxplot of the abundance for
each genotype and also prints the consensus genotype for musculus and
domesticus mice.

## Including Code

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
