setwd("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Cleaning_snps")

eigenvec <- read.table("plink.eigenvec")
eigenvec$family <- substr(eigenvec$V2, 1,7)
eigenvec$matingpair <- substr(eigenvec$V2, 1, 10)
eigenvec$id <- substr(eigenvec$V2, 9, 10)


library(ggplot2)
library(ggsci)

tron_p2 <- c("#FF410DFF", "36EE2FFFF", "#F7c530ff", "#95cc5eff", "#d0dfe6ff", "#f79d1eff", "#748aa6ff", "#ff99cc", "#4dd2ff")

#PC1 en PC2
ggplot(eigenvec)+ geom_point(aes(x=V3, y=V4, color=family)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC1_PC2_family.pdf")

ggplot(eigenvec)+ geom_point(aes(x=V3, y=V4, color=id, alpha=0.3)) + labs(x="PC1", y="PC2", color="") + theme_bw()+scale_color_manual(values=tron_p2) +facet_wrap(~family) + theme(legend.position = "none")
ggsave("PCA_PC1_PC2_facet_family.pdf")

#PC2 en PC3
ggplot(eigenvec)+ geom_point(aes(x=V4, y=V5, color=family)) + labs(x="PC2", y="PC3", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC2_PC3_family.pdf")

ggplot(eigenvec)+ geom_point(aes(x=V4, y=V5, color=id, alpha=0.3)) + labs(x="PC2", y="PC3", color="") + theme_bw()+scale_color_manual(values=tron_p2) +facet_wrap(~family) + theme(legend.position = "none")
ggsave("PCA_PC2_PC3_facet_family.pdf")

#PC3 en PC4
ggplot(eigenvec)+ geom_point(aes(x=V5, y=V6, color=family)) + labs(x="PC3", y="PC4", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC3_PC4_family.pdf")

ggplot(eigenvec)+ geom_point(aes(x=V5, y=V6, color=id, alpha=0.3)) + labs(x="PC3", y="PC4", color="") + theme_bw()+scale_color_manual(values=tron_p2) +facet_wrap(~family) + theme(legend.position = "none")
ggsave("PCA_PC3_PC4_facet_family.pdf")


#PC4 en PC5
ggplot(eigenvec)+ geom_point(aes(x=V6, y=V7, color=family)) + labs(x="PC4", y="PC5", color="") + theme_bw()+scale_color_manual(values=tron_p2)
ggsave("PCA_PC4_PC5_family.pdf")

ggplot(eigenvec)+ geom_point(aes(x=V6, y=V7, color=id, alpha=0.3)) + labs(x="PC4", y="PC5", color="") + theme_bw()+scale_color_manual(values=tron_p2) +facet_wrap(~family) + theme(legend.position = "none")
ggsave("PCA_PC4_PC5_facet_family.pdf")
