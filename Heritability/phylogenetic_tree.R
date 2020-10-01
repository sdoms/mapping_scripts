setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/snp_heritability/cospeciation/")
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(ggtree)
library(readxl)
library(tidytree)
library(viridis)
set.seed(21-09-20)
sequences <- read.table("SV_for_fasta.tsv")
colnames(sequences) <- c("SV", "sequence")

tax_table <- read.csv("../../../Phenotyping_27.02.2020/tables_core/otu_tax_table_f2_core.csv")
tax_table$name <- paste0(tax_table$X,"_", tax_table$Genus)
tax_table <- tax_table[,c(1,8)]
colnames(tax_table)<- c("SV", "Name")

seq2 <- merge(tax_table, sequences, by="SV")
rownames(seq2)<- seq2$Name
n <- as.character(seq2$sequence)
names(n) <- seq2$Name

alignment <- AlignSeqs(DNAStringSet(n), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
tree <- fitGTR$tree

herit <- read_excel("../heritability_cospeciation.xlsx")
herit <- herit[,c(1:6)]


herit$tip <- paste0(herit$SV, "_", herit$Genus)

x <- as_tibble(tree)

d <- tibble(label=herit$tip, trait= herit$DNA)
y <- full_join(x,d, by='label') 

y
tree2 <- as.treedata(y)
ggtree(tree2,branch.length = "none", aes(color=trait)) +geom_tiplab(size=2) + scale_color_viridis() + theme(legend.position="top") + ggplot2::xlim(0, 30)
ggsave("SNP_heritability_DNA_tree.pdf")

d <- tibble(label=herit$tip, trait= herit$RNA)
y <- full_join(x,d, by='label') 

y
tree2 <- as.treedata(y)
ggtree(tree2,branch.length = "none", aes(color=trait)) +geom_tiplab(size=2) + scale_color_viridis() + theme(legend.position="top") + ggplot2::xlim(0, 30)
ggsave("SNP_heritability_RNA_tree.pdf")

d <- tibble(label=herit$tip, trait= herit$`TM score`)
y <- full_join(x,d, by='label') 

y
tree2 <- as.treedata(y)
ggtree(tree2,branch.length = "none", aes(color=trait)) +geom_tiplab(size=2) + scale_color_viridis() + theme(legend.position="top") + ggplot2::xlim(0, 30)
ggsave("TM_score_tree.pdf")

d <- tibble(label=herit$tip, trait= herit$CospeciationScoreGenus)
y <- full_join(x,d, by='label') 

y
tree2 <- as.treedata(y)
ggtree(tree2,branch.length = "none", aes(color=trait)) +geom_tiplab(size=2) + scale_color_viridis() + theme(legend.position="top") + ggplot2::xlim(0, 30)
ggsave("Cospeciation_tree.pdf")
