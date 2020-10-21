setwd("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps")

library(data.table)
library(argyle)
library(tidyverse)
#6.10.20 added payseur data to PCA plots
#color_palettes
cbbPalette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", '#6a3d9a', "#ffff99", "#b15928", "#40e0d0", "#ffff00", "#ee82ee", "#fb8072", "#191970", "#fdb462", "#CBD4D4", "#fccde5", "#d9d9d9", "#bc80bd", "#32cd32", "#ffed6f", "lightcoral" )
colors_palette = c("#a6cee3", "#1f78b4", "#f0f8ff", "#b2df8a", "#33a02c", "#fb9a99", "#ff69b4", "#e31a1c", "#ffe4e1", "#cd853f", "#f0bf6f", "#ff7f00", "#cab2d6", "#ba3d9a","#6a3d94","#ffff99", "#B15928", "#ffb90f", "#8b4513")


load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/all_genotypes_all_samples.Rdata")
load("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/clean_snps.Rdata")
load("parents_all_snps_geno.RData")
load("~/Documents/PhD/Experiments/Final_QTL_mapping/Genotype_data/GM_snps.Rdata")



payseur <- fread("../Genotype_data/BretPayseurGigaMUGA_May2015_copy.csv")
payseur <- payseur[,-c(2:4)]
all_markers <- rownames(geno)
payseur <- payseur[Marker %in% all_markers,]
rownames(payseur)<- payseur$Marker

all <- cbind(geno, payseur)
all <- as.matrix(all,rownames="Marker")
GM_snps <- GM_snps[all_markers,]
# all snps #####




#make fam file
fam_hybrids <- fread("clean_f2.fam")
colnames(fam_hybrids) <- c("fid", "iid", "father", "mother", "sex", "pheno")
fam_hybrids$fid <- "HZ_F2"

# payseur 
sample_key <- fread("../Genotype_data/SampleKey.csv")
fam_payseur <- sample_key %>% 
  mutate(gender=ifelse(Sex=="M", 1,ifelse(Sex=="F", 2,0))) %>% 
  select(fid=Pop, iid=ColnameOrigGenoFile, sex=gender) %>% 
  mutate(father=0, mother=0, pheno=0)


# parents
fam_parents <- as.data.frame(colnames(parents))%>% 
  select(iid=`colnames(parents)`)
fam_parents$iid <- as.character(fam_parents$iid)
fam_parents$sex <- substr(fam_parents$iid, nchar(fam_parents$iid), nchar(fam_parents$iid))
fam_parents <- fam_parents %>% 
  mutate(sex=ifelse(sex=="M", 1,ifelse(sex=="F", 2,0)),fid="HZ_F0",father=0, mother=0, pheno=0)

#representative males
males_rep <- colnames(geno)[which(!(colnames(geno)%in% (c(fam_hybrids$iid,fam_parents$iid))))]

fam_males_rep <- as.data.frame(males_rep) %>% 
  select(iid=males_rep)

fam_males_rep$iid <- as.character(fam_males_rep$iid)
fam_males_rep$sex <- substr(fam_males_rep$iid, nchar(fam_males_rep$iid), nchar(fam_males_rep$iid))
fam_males_rep <- fam_males_rep %>% 
  mutate(sex=ifelse(sex=="M", 1,ifelse(sex=="F", 2,0)),fid="HZ_rep",father=0, mother=0, pheno=0)

#together

fam1 <- rbind(fam_parents, fam_hybrids)
fam2 <- rbind(fam1, fam_payseur)
fam <- rbind(fam2, fam_males_rep)
colnames(fam) <- c("iid", "sex", "fid", "dad", "mom", "pheno")
rownames(fam)<- fam$iid
#map file 
snps <- GM_snps %>% 
  select(chr,marker, cM, pos) %>% 
  drop_na(marker)
rownames(snps)<- snps$marker
#as genotypes object
all <- all[rownames(all)%in%snps$marker,]
all_geno <- genotypes(all, snps,alleles = "native" )
#write.plink(all_geno, "payseur_and_hybrids", map=snps, fam=fam)



pc <- pca(all_geno, K=5)
plot(pc, screeplot = T)
pc_cor <- merge(fam, pc[,c(1,7:11)], by="iid")
pc_cor$fid <- as.factor(pc_cor$fid)



library(ggplot2)
library(ggsci)

tron_p2 <- c("#FF410DFF", "36EE2FFFF", "#F7c530ff", "#95cc5eff", "#d0dfe6ff", "#f79d1eff", "#748aa6ff", "#ff99cc", "#4dd2ff")
colpal1 <- c("#26547C","#EF476F","#FFD166","#06D6A0","#90be6d","#43aa8b","#577590")
colpal2 <- c("e2e2df","d2d2cf","e2cfc4","f7d9c4","faedcb","c9e4de","c6def1","dbcdf0","f2c6de","f9c6c9")
colpal3 <- c("f94144","f3722c","f8961e","f9844a","f9c74f","90be6d","43aa8b","4d908e","577590","277da1")
colpal <- c("#e2e2df","#d2d2cf","#e2cfc4","#f7d9c4","#faedcb","#c9e4de","#c6def1","#dbcdf0","#f2c6de",
            "#f9c6c9","#f94144","#f3722c","#f8961e","#f9844a","#f9c74f","#90be6d","#43aa8b","#4d908e",
            "#577590","#277da1","#11151c","#212d40","#364156","#7d4e57","#d66853", "black")
colpal26 <- c("#df8730","#606cd4","#9fba3a","#975ed5","#519f2f","#df6dd0","#57c771","#9c40a2","#c8a83e",
              "#5e71b5","#d04827","#4cbfb2","#d4438c","#3d8756","#cd3955", "#5fa3da","#e77060","#99b36b",
              "#bd8dd5","#67752a","#9a4b7a","#8f6f2f","#e486a5","#a2552b","#dc986d","#ac585a")
colpal22 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', 
              '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', 
              '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')
shapes <- c(0,2,15,1,10,21,13,24)

#PC1 en PC2
ggplot(pc_cor)+ geom_point(aes(x=PC1, y=PC2, color=fid)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC1_PC2_hybrids_vs_payseur_all_snps.pdf")


info <- fread("../Genotype_data/HZsampleInfo.csv")
#add pop 
pc_cor <- pc_cor %>% 
  separate(iid,c("head", "tail"),"\\(") %>% 
  separate(head, c("head2, tail2","UW_"),remove = F)
pc_cor$head[370:439]<- pc_cor$UW_[370:439]
pc_cor$tail <- NULL
colnames(pc_cor)<- c("iid", "site", "HZ2", colnames(pc_cor)[4:13])
sites <- info[match(pc_cor$iid[370:439], info$SampleID),"Site"]
pc_cor$site[370:439]<- sites$Site
sites <- as.character(pc_cor[which(is.na(pc_cor$site)),"fid"])
pc_cor[which(is.na(pc_cor$site)),"site"] <-sites

pc_cor <- pc_cor %>% 
  mutate(subsp=ifelse(fid=="HZ","hyb",
                      ifelse(fid=="HZ_F2","hyb",
                                             ifelse(fid=="AL","mus", 
                                                                       ifelse(fid=="CR", "mus", 
                                                                                               ifelse(fid=="CB", "dom", 
                                                                                                      ifelse(fid=="MC","dom", ifelse(fid=="HZ_rep", "hyb",ifelse(fid=="HZ_F0", "hyb", NA)))))))))
write.csv(pc_cor, "PCA_val_all_snps.csv")

ggplot(pc_cor)+ geom_point(aes(x=PC1, y=PC2, color=site,shape=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_colour_manual(values = colpal26)
ggsave("PCA_all_site_all_snps.pdf")
ggplot(pc_cor)+ geom_point(aes(x=PC1, y=PC2, color=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw() +scale_colour_manual(values = colpal1)
ggsave("PCA_all_subsp_all_snps.pdf")
#PC2 en PC3
ggplot(pc_cor)+ geom_point(aes(x=PC2, y=PC3, color=site, fill=site,shape=fid)) + labs(x="PC2", y="PC3", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC2_PC3_site_all_snps.pdf")

#PC3 en PC4
ggplot(pc_cor)+ geom_point(aes(x=PC3, y=PC4, color=site, fill=site,shape=fid)) + labs(x="PC3", y="PC4", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC3_PC4_site_all_snps.pdf")

#PC4 en PC5
ggplot(pc_cor)+ geom_point(aes(x=PC4, y=PC5, color=site, fill=site,shape=fid)) + labs(x="PC4", y="PC5", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC4_PC5_site_all_snps.pdf")

# only autosomes

geno_auto <- autosomes(all_geno)

pc_auto <- pca(geno_auto, K=5)

pc_cor_auto <- pc_cor %>% 
  select(iid, site, sex, fid, subsp) %>% 
  mutate(pc_auto[7:11])
write.csv(pc_cor_auto, "PCA_val_all_snps_autosomes.csv")

#PC1 en PC2
ggplot(pc_cor_auto)+ geom_point(aes(x=PC1, y=PC2, color=fid)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC1_PC2_hybrids_vs_payseur_all_snps_autosome.pdf")

ggplot(pc_cor_auto)+ geom_point(aes(x=PC1, y=PC2, color=site,shape=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_colour_manual(values = colpal26)
ggsave("PCA_all_site_all_snps_autosome.pdf")
ggplot(pc_cor_auto)+ geom_point(aes(x=PC1, y=PC2, color=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw() +scale_colour_manual(values = colpal1)
ggsave("PCA_all_subsp_all_snps_autosome.pdf")
#PC2 en PC3
ggplot(pc_cor_auto)+ geom_point(aes(x=PC2, y=PC3, color=site, fill=site,shape=fid)) + labs(x="PC2", y="PC3", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC2_PC3_site_all_snps_autosome.pdf")

#PC3 en PC4
ggplot(pc_cor_auto)+ geom_point(aes(x=PC3, y=PC4, color=site, fill=site,shape=fid)) + labs(x="PC3", y="PC4", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC3_PC4_site_all_snps_autosome.pdf")

#PC4 en PC5
ggplot(pc_cor_auto)+ geom_point(aes(x=PC4, y=PC5, color=site, fill=site,shape=fid)) + labs(x="PC4", y="PC5", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC4_PC5_site_all_snps_autosome.pdf")

# only autosomes

geno_X <-xchrom(all_geno)

pc_X <- pca(geno_X, K=5)

pc_cor_X <- pc_cor %>% 
  select(iid, site, sex, fid, subsp) %>% 
  mutate(pc_X[7:11])
write.csv(pc_cor_X, "PCA_val_all_snps_Xchrom.csv")

#PC1 en PC2
ggplot(pc_cor_X)+ geom_point(aes(x=PC1, y=PC2, color=fid)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC1_PC2_hybrids_vs_payseur_all_snps_Xchrom.pdf")

ggplot(pc_cor_X)+ geom_point(aes(x=PC1, y=PC2, color=site,shape=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_colour_manual(values = colpal26)
ggsave("PCA_all_site_all_snps_Xchrom.pdf")
ggplot(pc_cor_X)+ geom_point(aes(x=PC1, y=PC2, color=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw() +scale_colour_manual(values = colpal1)
ggsave("PCA_all_subsp_all_snps_Xchrom.pdf")
#PC2 en PC3
ggplot(pc_cor_X)+ geom_point(aes(x=PC2, y=PC3, color=site, fill=site,shape=fid)) + labs(x="PC2", y="PC3", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC2_PC3_site_all_snps_Xchrom.pdf")

#PC3 en PC4
ggplot(pc_cor_X)+ geom_point(aes(x=PC3, y=PC4, color=site, fill=site,shape=fid)) + labs(x="PC3", y="PC4", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC3_PC4_site_all_snps_Xchrom.pdf")

#PC4 en PC5
ggplot(pc_cor_X)+ geom_point(aes(x=PC4, y=PC5, color=site, fill=site,shape=fid)) + labs(x="PC4", y="PC5", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC4_PC5_site_all_snps_Xchrom.pdf")


###### cleaned snps #####
load("./clean_snps.Rdata")

clean_geno <- subset(all_geno, markers(all_geno)$marker%in%snps$marker)

clean_geno_auto <- autosomes(clean_geno)
clean_geno_Xchrom <- xchrom(clean_geno)
pc_cor <- read.csv("PCA_val_all_snps.csv")

pc_auto_clean <- pca(clean_geno_auto, K=5)

pc_cor_auto_clean <- pc_cor %>% 
  select(iid, site, sex, fid, subsp) %>% 
  mutate(pc_auto_clean[7:11])
write.csv(pc_cor_auto_clean, "PCA_val_clean_snps_autosomes.csv")

#PC1 en PC2
ggplot(pc_cor_auto_clean)+ geom_point(aes(x=PC1, y=PC2, color=fid)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC1_PC2_hybrids_vs_payseur_clean_snps_autosome.pdf")

ggplot(pc_cor_auto_clean)+ geom_point(aes(x=PC1, y=PC2, color=site,shape=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_colour_manual(values = colpal26)
ggsave("PCA_all_site_clean_snps_autosome.pdf")
ggplot(pc_cor_auto_clean)+ geom_point(aes(x=PC1, y=PC2, color=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw() +scale_colour_manual(values = colpal1)
ggsave("PCA_all_subsp_clean_snps_autosome.pdf")
#PC2 en PC3
ggplot(pc_cor_auto_clean)+ geom_point(aes(x=PC2, y=PC3, color=site, fill=site,shape=fid)) + labs(x="PC2", y="PC3", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC2_PC3_site_clean_snps_autosome.pdf")

#PC3 en PC4
ggplot(pc_cor_auto_clean)+ geom_point(aes(x=PC3, y=PC4, color=site, fill=site,shape=fid)) + labs(x="PC3", y="PC4", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC3_PC4_site_clean_snps_autosome.pdf")

#PC4 en PC5
ggplot(pc_cor_auto_clean)+ geom_point(aes(x=PC4, y=PC5, color=site, fill=site,shape=fid)) + labs(x="PC4", y="PC5", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC4_PC5_site_clean_snps_autosome.pdf")

# only autosomes


pc_X_clean <- pca(clean_geno_Xchrom, K=5)

pc_cor_X_clean <- pc_cor %>% 
  select(iid, site, sex, fid, subsp) %>% 
  mutate(pc_X_clean[7:11])
write.csv(pc_cor_X_clean, "PCA_val_clean_snps_Xchrom.csv")

#PC1 en PC2
ggplot(pc_cor_X_clean)+ geom_point(aes(x=PC1, y=PC2, color=fid)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC1_PC2_hybrids_vs_payseur_clean_snps_Xchrom.pdf")

ggplot(pc_cor_X_clean)+ geom_point(aes(x=PC1, y=PC2, color=site,shape=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_colour_manual(values = colpal26)
ggsave("PCA_all_site_clean_snps_Xchrom.pdf")
ggplot(pc_cor_X_clean)+ geom_point(aes(x=PC1, y=PC2, color=subsp)) + labs(x="PC1", y="PC2", color="") + theme_bw() +scale_colour_manual(values = colpal1)
ggsave("PCA_all_subsp_clean_snps_Xchrom.pdf")
#PC2 en PC3
ggplot(pc_cor_X_clean)+ geom_point(aes(x=PC2, y=PC3, color=site, fill=site,shape=fid)) + labs(x="PC2", y="PC3", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC2_PC3_site_clean_snps_Xchrom.pdf")

#PC3 en PC4
ggplot(pc_cor_X_clean)+ geom_point(aes(x=PC3, y=PC4, color=site, fill=site,shape=fid)) + labs(x="PC3", y="PC4", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC3_PC4_site_clean_snps_Xchrom.pdf")

#PC4 en PC5
ggplot(pc_cor_X_clean)+ geom_point(aes(x=PC4, y=PC5, color=site, fill=site,shape=fid)) + labs(x="PC4", y="PC5", color="", fill="", shape="") + 
  theme_bw()+scale_color_manual(values=colpal26)+scale_shape_manual(values=shapes) +scale_fill_manual(values=colpal26)
ggsave("PCA_PC4_PC5_site_clean_snps_Xchrom.pdf")
