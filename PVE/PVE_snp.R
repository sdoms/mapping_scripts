setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/PVE_snps")

library(readxl)
library(tidyverse)
library(ggpmisc)
load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/sig.summaries/all_DNA.Rdata")
load("../DNA/sig.summaries/otu/SV35_gwscan_allchr.Rdata")

# to get the counts 
input<- all_files[,1:23] %>% 
  left_join(SV35_gwscan[,1:12], by=c("peak.snp"="marker"))

input$MAF<- ((2*input$AA/input$n) +input$AB/input$n)/2

# additive effect ####
input$add.num <- 2*(input$add.Beta_peak_snp^2)*input$MAF * (1-input$MAF)
input$add.denom <- 2*(input$add.Beta_peak_snp^2)*input$MAF * (1-input$MAF) + (input$add.StdErr_peak_snp^2)* 2 * input$n *input$MAF * (1-input$MAF)

input$add.pve.lead.snp <- input$add.num/input$add.denom

input$dom.num <- 2*(input$dom.Beta_peak_snp^2)*input$MAF * (1-input$MAF)
input$dom.denom <- 2*(input$dom.Beta_peak_snp^2)*input$MAF * (1-input$MAF) + (input$dom.StdErr_peak_snp^2)* 2 * input$n *input$MAF * (1-input$MAF)

input$dom.pve.lead.snp <- input$dom.num/input$dom.denom

ggplot(input, aes(x=add.pve.lead.snp, y=dom.pve.lead.snp))+geom_point() +geom_smooth(method = "lm")+
  stat_poly_eq( formula=y ~ x,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)

# calculate per tax 
input.per.taxa <- input %>% 
  group_by(trait) %>% 
  distinct(chr.num, start.LD.pos, stop.LD.pos, .keep_all=T) %>% 
  distinct(peak.snp, .keep_all=T) %>% 
  select(-tax, -index, -pos, -chr.x, -chr.y, -cM,-X1) %>% 
  ungroup()
sum.input.per.taxa <-  input.per.taxa %>% 
  group_by(trait) %>% 
  summarise(sum(add.num))
  

###### heritabilit ######
herit <- read_delim("../../snp_heritability/lme4qtl/sign_heritability_dna.csv", delim=";")

herit.pve <- herit %>% 
  inner_join(sum.input.per.taxa, by=c("taxa"="trait"))
ggplot(herit.pve, aes(x=h2_id, y=`sum(add.num)`)) + geom_point() +geom_smooth(method="lm" ) +geom_text(aes(label=taxa))
