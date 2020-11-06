setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps/")
library(flashpcaR)
library(reshape2)
library(ggplot2)
library(ggsci)
library(patchwork)
fl<-flashpca("../Cleaning_snps/clean_f2", ndim=10, stand="binom2", do_loadings = T  )

meta <- read.table("../Cleaning_snps/clean_f2.fam")
meta$V2 <- gsub("/", ".", meta$V2)

rownames(fl$vectors)<- meta$V2
rownames(fl$projection)<- meta$V2
colnames(fl$projection)<- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
colnames(fl$vectors)<- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

dat <- fl$projection
dat <- cbind(dat,substr(rownames(dat), 1,7))
dat <- cbind(dat,substr(rownames(dat), 1,10))
dat <- cbind(dat,substr(rownames(dat), 1,12))
colnames(dat)<- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "family", "group", "litter")
dat <- as.data.frame(dat)
dat$PC1<- as.numeric(dat$PC1)
dat$PC2<- as.numeric(dat$PC2)
dat$PC3<- as.numeric(dat$PC3)
dat$PC4<- as.numeric(dat$PC4)
dat$PC5<- as.numeric(dat$PC5)
dat$PC6<- as.numeric(dat$PC6)
dat$PC7<- as.numeric(dat$PC7)
dat$PC8<- as.numeric(dat$PC8)
dat$PC9<- as.numeric(dat$PC9)
dat$PC10<- as.numeric(dat$PC10)
save(fl, file="fl_pca_family.Rdata")
write.table(fl$vectors, "PC1-10.txt")
p <- ggplot(data=dat, aes(x=PC1, y=PC2)) + geom_point(aes(color=family)) + theme_test() +  xlab(paste0("PC1 (", sprintf("%.1f", 100*fl$pve[1]), "%)")) + 
  ylab(paste0("PC2 (", sprintf("%.1f", 100*fl$pve[2]), "%)")) +scale_color_d3()+labs(col="")
p 
ggsave("PC1_PC2_hybrids_family.pdf")

q <- ggplot(data=dat, aes(x=PC2, y=PC3)) + geom_point(aes(color=family)) + theme_test() +  xlab(paste0("PC2 (", sprintf("%.1f", 100*fl$pve[2]), "%)")) + 
  ylab(paste0("PC3 (", sprintf("%.1f", 100*fl$pve[3]), "%)")) +scale_color_d3()+labs(col="")

q
ggsave("PC2_PC3_hybrids_family.pdf")


r <- ggplot(data=dat, aes(x=PC3, y=PC4)) + geom_point(aes(color=family)) + theme_test() +  xlab(paste0("PC3 (", sprintf("%.1f", 100*fl$pve[3]), "%)")) + 
  ylab(paste0("PC4 (", sprintf("%.1f", 100*fl$pve[4]), "%)")) +scale_color_d3()+labs(col="")
 
r
ggsave("PC3_PC4_hyrbids_family.pdf")

s <- ggplot(data=dat, aes(x=PC4, y=PC5)) + geom_point(aes(color=family)) + theme_test() +  xlab(paste0("PC4 (", sprintf("%.1f", 100*fl$pve[4]), "%)")) + 
  ylab(paste0("PC5 (", sprintf("%.1f", 100*fl$pve[5]), "%)")) +scale_color_d3()+labs(col="")

s
ggsave("PC4_PC5_hybrids_family.pdf")


(p+q)/(r+s)+plot_layout(guides="collect")+plot_annotation(tag_levels = "A")
ggsave("PCA_hybrids_family_PC1-PC5.pdf")


t <- ggplot(data=dat, aes(x=PC5, y=PC6)) + geom_point(aes(color=family)) + theme_test() +  xlab(paste0("PC1 (", sprintf("%.1f", 100*fl$pve[1]), "%)")) + 
  ylab(paste0("PC2 (", sprintf("%.1f", 100*fl$pve[2]), "%)")) +scale_color_d3()+labs(col="")
p 
t
ggsave("PC5_PC6_family.pdf")

u <- ggplot(data=dat, aes(x=PC6, y=PC7)) + geom_point(aes(color=family)) 
u
ggsave("PC6_PC7_family.pdf")

v <- ggplot(data=dat, aes(x=PC7, y=PC8)) + geom_point(aes(color=family)) 
v
ggsave("PC7_PC8_family.pdf")

w <- ggplot(data=dat, aes(x=PC8, y=PC9)) + geom_point(aes(color=family)) 
w
ggsave("PC8_PC9_family.pdf")

x <- ggplot(data=dat, aes(x=PC9, y=PC10)) + geom_point(aes(color=family)) 
x
ggsave("PC9_PC10_family.pdf")

library(corrplot)
