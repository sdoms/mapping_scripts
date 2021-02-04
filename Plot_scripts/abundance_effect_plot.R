source('~/Documents/PhD/Experiments/Final_QTL_mapping/Scripts/Plot_scripts/Effect_plots.R')


genotypes_ref <- data.frame()
for (marker in rownames(clean_geno)){
  all_mark<- clean_geno[marker,]
  a1 <- snps[marker, "A1"]
  a2<- snps[marker, "A2"]
  if (is.na(a1)|is.na(a2)){
    tt<- table(as.factor(all_mark))
    genotypes_ref[marker, "A1_name"]<- names(which.min(tt))
    genotypes_ref[marker, "A1_num"]<- min(tt)
    genotypes_ref[marker, "A2_name"]<- names(which.max(tt))
    genotypes_ref[marker, "A2_num"]<- max(tt)
    genotypes_ref[marker, "AA_min"]<- names(which.min(tt))
    genotypes_ref[marker, "BB_maj"]<- names(which.max(tt))
  } else{
    
    
    tt<- table(factor(all_mark, levels=c(a1,"H", a2)))
    
    a1_num <- tt[a1]
    a2_num<- tt[a2]
    genotypes_ref[marker, "A1_name"]<- a1
    genotypes_ref[marker, "A1_num"]<- a1_num
    genotypes_ref[marker, "A2_name"]<- a2
    genotypes_ref[marker, "A2_num"]<- a2_num
    if (a1_num>a2_num){
      genotypes_ref[marker, "AA_min"]<- names(a2_num)
      genotypes_ref[marker, "BB_maj"]<- names(a1_num)
    } else if (a2_num>a1_num){
      genotypes_ref[marker, "AA_min"]<- names(a1_num)
      genotypes_ref[marker, "BB_maj"]<- names(a2_num)
      
    } else if (a1_num==a2_num){
      genotypes_ref[marker, "AA_min"]<- a1
      genotypes_ref[marker, "BB_maj"]<-a2
    }
  }
  
}

genotypes_ref$marker <- rownames(genotypes_ref)
save(genotypes_ref, file="ref_alt_min_major.Rdata")
allP_all_minmaj<- merge(allP_all, genotypes_ref, by="marker")
save(allP_all_minmaj, file="allP_all_minmaj.Rdata")

write.csv(allP_all_minmaj, "allP_all_minmaj.csv",sep = ";" )

setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Plots/Effect_plots/")
abundance.plot("SV1", "otu", "RNA", "UNC11806630")
abundance.plot("SV1", "otu", "RNA", "UNC8462988")
abundance.plot("SV2", "otu", "RNA", "UNC24866496")
abundance.plot("SV2", "otu", "RNA", "UNCHS032891")
abundance.plot("SV19", "otu", "RNA", "UNC3694204")
abundance.plot("SV29", "otu", "RNA", "UNC26246433")
abundance.plot("SV49", "otu", "RNA", "JAX00150794")
abundance.plot("SV52", "otu", "RNA", "UNCHS001411")
abundance.plot("SV91", "otu", "RNA", "UNCHS011150")
abundance.plot("SV36", "otu", "RNA", "JAX00633860")
abundance.plot("SV36", "otu", "DNA", "JAX00633860")
abundance.plot("G_Oscillibacter", "genus", "RNA", "JAX00633860")
abundance.plot("G_Acetatifactor", "genus", "DNA", "UNCHS006872")
abundance.plot("SV142", "otu", "RNA", "UNCHS006872")
abundance.plot("SV29", "otu", "DNA", "UNC24166386")
abundance.plot("SV17", "otu", "RNA", "UNC1045279")


# UNCHS039968
abundance.plot("SV197", "otu", "DNA", "UNCHS039968")
abundance.plot("SV243", "otu", "DNA", "UNCHS039968")
abundance.plot("SV139", "otu", "DNA", "UNCHS039968")
abundance.plot("SV61", "otu", "DNA", "UNCHS039968")
abundance.plot("SV173", "otu", "DNA", "UNCHS039968")
abundance.plot("SV79", "otu", "DNA", "UNCHS039968")
abundance.plot("SV246", "otu", "DNA", "UNCHS039968")
abundance.plot("G_Anaerostipes", "genus", "DNA", "UNCHS039968")
abundance.plot("unclassified_F_Prevotellaceae", "genus", "DNA", "UNCHS039968")
abundance.plot("G_Roseburia", "genus", "DNA", "UNCHS039968")

abundance.plot("SV197", "otu", "RNA", "UNCHS039968")
abundance.plot("SV243", "otu", "RNA", "UNCHS039968")
abundance.plot("SV139", "otu", "RNA", "UNCHS039968")
abundance.plot("SV61", "otu", "RNA", "UNCHS039968")
abundance.plot("SV173", "otu", "RNA", "UNCHS039968")
abundance.plot("SV79", "otu", "RNA", "UNCHS039968")
abundance.plot("SV246", "otu", "RNA", "UNCHS039968")
abundance.plot("G_Anaerostipes", "genus", "RNA", "UNCHS039968")
abundance.plot("unclassified_F_Prevotellaceae", "genus", "RNA", "UNCHS039968")
abundance.plot("G_Roseburia", "genus", "RNA", "UNCHS039968")


# UNC11687867
abundance.plot("G_Lactobacillus", "genus", "RNA", "UNC11687867")

# UNC11691749
abundance.plot("G_Lactobacillus", "genus", "RNA", "UNC11691749")

# slc26a3 UNC20836314
abundance.plot("unclassified_C_Deltaproteobacteria", "genus", "DNA", "UNC20836314")
abundance.plot("unclassified_C_Deltaproteobacteria", "genus", "DNA", "UNCJPD004947")

abundance.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNC20836314")
abundance.plot("F_Staphylococcaceae", "family", "DNA", "UNC20836314")
abundance.plot("F_Staphylococcaceae", "family", "RNA", "UNC20836314")
abundance.plot("F_Staphylococcaceae", "family", "RNA", "UNCHS033197")

abundance.plot("SV26", "otu", "DNA", "UNC20836314")
abundance.plot("SV26", "otu", "DNA", "UNC20865798")

abundance.plot("SV29", "otu", "DNA", "UNC20836314")
abundance.plot("SV29", "otu", "DNA", "UNC20863257")

abundance.plot("SV88", "otu", "DNA", "UNC20836314")


abundance.plot("SV396", "otu", "RNA", "UNC20836314")
abundance.plot("SV396", "otu", "RNA", "UNC20861803")
abundance.plot("G_Acetatifactor", "genus", "RNA", "UNC20836314")
abundance.plot("G_Acetatifactor", "genus", "RNA", "UNCJPD004947")


#### adcy2
abundance.plot("SV36", "otu", "RNA", "UNCHS036333")

# Foxp1
#DNA: F_Rikenellaceae | SV264 | SV20 | G_Alistipes | SV31 | SV63 | SV248 | SV21 
# RNA: G_Lactobacillus | O_Lactobacillales | C_Bacilli | SV31 | SV181 | SV20 | F_Lactobacillaceae | O_Bacillales
abundance.plot("SV20", "otu", "DNA", "UNC11687867")
abundance.plot("SV31", "otu", "DNA", "UNC11687867")
abundance.plot("SV63", "otu", "DNA", "UNC11687867")
abundance.plot("SV248", "otu", "DNA", "UNC11687867")
abundance.plot("SV21", "otu", "DNA", "UNC11687867")
abundance.plot("F_Rikenellaceae", "family", "DNA", "UNC11687867")
abundance.plot("G_Alistipes", "genus", "DNA", "UNC11687867")

abundance.plot("SV31", "otu", "RNA", "UNC11687867")
abundance.plot("SV181", "otu", "RNA", "UNC11687867")
abundance.plot("SV20", "otu", "RNA", "UNC11687867")
abundance.plot("G_Lactobacillus", "genus", "RNA", "UNC11687867")
abundance.plot("O_Lactobacillales", "order", "RNA", "UNC11687867")
abundance.plot("C_Bacilli", "class", "RNA", "UNC11687867")
abundance.plot("F_Lactobacillaceae", "family", "RNA", "UNC11687867")
abundance.plot("O_Bacillales", "order", "RNA", "UNC11687867")


# slc3a1 UNCHS045125
# DNA: SV7 | G_Paraprevotella | F_Prevotellaceae | RNA: F_Prevotellaceae | SV7 | G_Paraprevotella

abundance.plot("SV7", "otu", "DNA", "UNCHS045125")
abundance.plot("G_Paraprevotella", "genus", "DNA", "UNCHS045125")
abundance.plot("F_Prevotellaceae", "family", "DNA", "UNCHS045125")
abundance.plot("SV7", "otu", "RNA", "UNCHS045125")
abundance.plot("G_Paraprevotella", "genus", "RNA", "UNCHS045125")
abundance.plot("F_Prevotellaceae", "family", "RNA", "UNCHS045125")


#### effect plots for z scores ####
abundance.plot("SV184","otu", "RNA","UNC7414459")
abundance.plot("SV184","otu", "RNA","UNC26145702")
# only one dom ~ 0 dom.t=0.19
effect.violin.plot("SV91","otu", "RNA","UNCHS043099")
effect.boxplot("SV91","otu", "RNA","UNCHS043099")
# dom.T =-0.81
effect.violin.plot("unclassified_C_Deltaproteobacteria", "genus", "RNA", "UNC26246433")
abundance.plot("unclassified_C_Deltaproteobacteria", "genus", "RNA", "UNC26246433")

# add.T = 0.36
effect.violin.plot("SV264", "otu", "RNA", "UNCHS047355")

# degree of dominance -0.41, add.T=6.5, dom.T=-2.5
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNC6440308")
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNCJPD001537")


# 
effect.violin.plot("G_Odoribacter", "genus", "RNA", "JAX00260563")
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNC6440308")
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "JAX00499850")
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "RNA", "UNC6454815")
effect.violin.plot("SV52", "otu", "RNA", "JAX00260563")
effect.violin.plot("G_Odoribacter", "genus", "RNA", "UNC6535447")
effect.violin.plot("G_Odoribacter", "genus", "RNA", "UNC6535901")
effect.violin.plot("SV142", "otu", "DNA", "UNC7121581")
effect.violin.plot("unclassified_C_Deltaproteobacteria", "genus", "RNA", "UNC26246433")
effect.violin.plot("G_Eisenbergiella", "genus", "RNA", "ICR1886")
effect.violin.plot("SV264", "otu", "RNA", "UNCHS047355")
effect.violin.plot("SV91", "otu", "RNA", "UNCHS043099")
effect.violin.plot("SV91", "otu", "RNA", "UNCHS011150")
effect.violin.plot("SV264", "otu", "RNA", "UNC20743110")
effect.violin.plot("SV145", "otu", "DNA", "JAX00214148")
effect.violin.plot("SV145", "otu", "DNA", "UNCJPD008271")
effect.violin.plot("SV17", "otu", "DNA", "UNCJPD001371")
effect.violin.plot("SV142", "otu", "DNA", "UNCJPD001666")
effect.violin.plot("SV142", "otu", "RNA", "UNCJPD001666")
effect.violin.plot("SV49", "otu", "RNA", "UNC2393187")

#### Maj > H > min
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNC6440308")
effect.boxplot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNC6440308")
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNCJPD001537")
effect.boxplot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "UNCJPD001537")
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "JAX00499850")
effect.boxplot("unclassified_F_Porphyromonadaceae", "genus", "DNA", "JAX00499850")
effect.violin.plot("unclassified_F_Porphyromonadaceae", "genus", "RNA", "UNC26246433")
effect.boxplot("unclassified_F_Porphyromonadaceae", "genus", "RNA", "UNC26246433")
# Maj > H > min deg~0
effect.violin.plot("SV91", "otu", "RNA", "UNCHS043099")
effect.boxplot("SV91", "otu", "RNA", "UNCHS043099")

#min >H>maj
effect.violin.plot("SV264", "otu", "RNA", "UNC20743110")
effect.boxplot("SV264", "otu", "RNA", "UNC20743110")
effect.violin.plot("SV21", "otu", "DNA", "JAX00277190")
effect.boxplot("SV21", "otu", "DNA", "JAX00277190")
effect.violin.plot("SV7", "otu", "DNA", "UNCHS014190")
effect.boxplot("SV7", "otu", "DNA", "UNCHS014190")
effect.violin.plot("C_Bacilli", "class", "DNA", "UNC29132165")
effect.boxplot("C_Bacilli", "class", "DNA", "UNC29132165")
effect.violin.plot("C_Deferribacteres", "class", "DNA", "UNC10048917")
effect.boxplot("C_Deferribacteres", "class", "DNA", "UNC10048917")
effect.violin.plot("O_Bacillales", "order", "DNA", "UNC20972762")
effect.boxplot("O_Bacillales", "order", "DNA", "UNC20972762")

#min > maj > H
effect.violin.plot("G_Lactobacillus", "genus", "RNA", "UNC11691749")
effect.boxplot("G_Lactobacillus", "genus", "RNA", "UNC11691749")
effect.violin.plot("SV7", "otu", "DNA", "JAX00312557")
effect.boxplot("SV7", "otu", "DNA", "JAX00312557")
#JAX00028034
effect.violin.plot("G_Dorea", "genus", "RNA", "JAX00028034")
effect.boxplot("G_Dorea", "genus", "RNA", "JAX00028034")
effect.violin.plot("G_Dorea", "genus", "DNA", "JAX00028034")
effect.boxplot("G_Dorea", "genus", "DNA", "JAX00028034")
effect.violin.plot("SV184", "otu", "RNA", "JAX00028034")
effect.boxplot("SV184", "otu", "RNA", "JAX00028034")
effect.violin.plot("SV293", "otu", "RNA", "JAX00028034")
effect.boxplot("SV293", "otu", "RNA", "JAX00028034")
