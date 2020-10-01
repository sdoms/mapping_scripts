setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/snp_heritability/cospeciation/")

library(readxl)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpmisc)
library(patchwork)
#color_palettes
cbbPalette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", '#6a3d9a', "#ffff99", "#b15928", "#40e0d0", "#ffff00", "#ee82ee", "#fb8072", "#191970", "#fdb462", "#CBD4D4", "#fccde5", "#d9d9d9", "#bc80bd", "#32cd32", "#ffed6f", "lightcoral" )
colors_palette = c("#a6cee3", "#1f78b4", "#f0f8ff", "#b2df8a", "#33a02c", "#fb9a99", "#ff69b4", "#e31a1c", "#ffe4e1", "#cd853f", "#f0bf6f", "#ff7f00", "#cab2d6", "#ba3d9a","#6a3d94","#ffff99", "#B15928", "#ffb90f", "#8b4513")

# genus ####
genus <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 2,col_names =T )


# Plot
genus$RNA <- genus$RNA*-1
gen_melt <-melt(genus, id.vars = "genus", measure.vars = c("DNA", "RNA"))
# X Axis Breaks and Labels 
brks <- seq(-1, 1, 0.2)
lbls = as.character(c(seq(1, 0, -0.2), seq(0.2, 1, 0.2)))
tot_gen <-ggplot(gen_melt, aes(x = reorder(genus, abs(value)), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(y="SNP heritability", x="", fill="") +
  theme_bw() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks.y  = element_blank()) +   # Centre plot title
  scale_fill_brewer(palette = "Dark2")  # Color palette
tot_gen
ggsave("genus_SNP_heritability.pdf", plot=tot_gen)


# cospeciation
genus <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 2,col_names =T )
cospec_genus <-genus[complete.cases(genus),]

#Total DNA
gen_DNA <-ggplot(cospec_genus, aes(x=CospeciationScore, y=DNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.2, size = 3)+ geom_text_repel(aes(label = genus), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.85, parse = TRUE, size = 3)+theme_bw() + ggtitle("DNA")
gen_DNA
gen_RNA <-ggplot(cospec_genus, aes(x=CospeciationScore, y=RNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.4, size = 3)+ geom_text_repel(aes(label = genus), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.85, parse = TRUE, size = 3)+theme_bw() + ggtitle("RNA")
gen_RNA

library(cowplot)
p1<-plot_grid(gen_DNA,gen_RNA, nrow = 1, align = "v")
plot_grid(tot_gen, p1, ncol = 1,  axis="l")
ggsave("Cospeciation_DNA_Genus.pdf", height = 8, width=10)


# TM score

# Plot
genus$RNA <- genus$RNA*-1
gen_melt <-melt(genus, id.vars = "genus", measure.vars = c("DNA", "RNA"))
# X Axis Breaks and Labels 
brks <- seq(-1, 1, 0.2)
lbls = as.character(c(seq(1, 0, -0.2), seq(0.2, 1, 0.2)))
tot_gen <-ggplot(gen_melt, aes(x = reorder(genus, abs(value)), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(y="SNP heritability", x="", fill="") +
  theme_bw() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks.y  = element_blank()) +   # Centre plot title
  scale_fill_brewer(palette = "Dark2")  # Color palette
tot_gen
ggsave("genus_SNP_heritability.pdf", plot=tot_gen)


# cospeciation
genus <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 2,col_names =T )
genus$CospeciationScore <- NULL
TM_genus <-genus[complete.cases(genus),]

#Total DNA
gen_DNA <-ggplot(TM_genus, aes(x=`TM score`, y=DNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="TM score", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.2, size = 3)+ geom_text_repel(aes(label = genus), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.85, parse = TRUE, size = 3)+theme_bw() + ggtitle("DNA")
gen_DNA
gen_RNA <-ggplot(TM_genus, aes(x=`TM score`, y=RNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="TM score", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.4, size = 3)+ geom_text_repel(aes(label = genus), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.85, parse = TRUE, size = 3)+theme_bw() + ggtitle("RNA")
gen_RNA

library(cowplot)
p1<-plot_grid(gen_DNA,gen_RNA, nrow = 1, align = "v")
plot_grid(tot_gen, p1, ncol = 1,  axis="l")
ggsave("TM_score_Genus.pdf", height = 8, width=10)

gen_DNA
ggsave("TM_score_Genus_DNA.pdf")

gen_RNA
ggsave("TM_score_Genus_RNA.pdf")

# FAMILY ####
family <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 3,col_names =T )


# Plot
family$RNA <- family$RNA*-1
fam_melt <-melt(family, id.vars = "family", measure.vars = c("DNA", "RNA"))
# X Axis Breaks and Labels 
brks <- seq(-1, 1, 0.2)
lbls = as.character(c(seq(1, 0, -0.2), seq(0.2, 1, 0.2)))
tot_fam <-ggplot(fam_melt, aes(x = reorder(family, abs(value)), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(y="SNP heritability", x="", fill="") +
  theme_bw() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks.y  = element_blank()) +   # Centre plot title
  scale_fill_brewer(palette = "Dark2")  # Color palette
ggsave("family_SNP_heritability.pdf", plot = tot_fam)


# cospeciation
family <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 3,col_names =T )
cospec_family <-family[complete.cases(family),]

#Total DNA
fam_DNA <-ggplot(cospec_family, aes(x=cospeciation_score, y=DNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.3, size = 3)+ geom_text_repel(aes(label = family), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.9, parse = TRUE, size = 3)+theme_bw() + ggtitle("DNA")
fam_DNA
fam_RNA <-ggplot(cospec_family, aes(x=cospeciation_score, y=RNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.5, size = 3)+ geom_text_repel(aes(label = family), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.9, parse = TRUE, size = 3)+theme_bw() + ggtitle("RNA")
fam_RNA

p1<-plot_grid(fam_DNA,fam_RNA, nrow = 1, align = "v")
plot_grid(tot_fam, p1, ncol = 1,  axis="l")
ggsave("Cospeciation_DNA_family.pdf", height = 8, width=10)

# order ####
order <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 4,col_names =T )


# Plot
order$RNA <- order$RNA*-1
ord_melt <-melt(order, id.vars = "order", measure.vars = c("DNA", "RNA"))
# X Axis Breaks and Labels 
brks <- seq(-1, 1, 0.2)
lbls = as.character(c(seq(1, 0, -0.2), seq(0.2, 1, 0.2)))
tot_ord <-ggplot(ord_melt, aes(x = reorder(order, abs(value)), y = value, fill = variable)) +   # 
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(y="SNP heritability", x="", fill="") +
  theme_bw() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks.y  = element_blank()) +   # Centre plot title
  scale_fill_brewer(palette = "Dark2")  # Color palette
ggsave("order_SNP_heritability.pdf", plot = tot_ord)


# cospeciation
order <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 4,col_names =T )
cospec_order <-order[complete.cases(order),]

#Total DNA
ord_DNA <-ggplot(cospec_order, aes(x=cospeciation, y=DNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.15, size = 3)+ geom_text_repel(aes(label = order), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.9, parse = TRUE, size = 3)+theme_bw() + ggtitle("DNA")
ord_DNA
ord_RNA <-ggplot(cospec_order, aes(x=cospeciation, y=RNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="Cospeciation rate", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 1, size = 3)+ geom_text_repel(aes(label = order), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.85, parse = TRUE, size = 3)+theme_bw() + ggtitle("RNA")
ord_RNA

p1<-plot_grid(ord_DNA,ord_RNA, nrow = 1, align = "v")
plot_grid(tot_ord, p1, ncol = 1,  axis="l")
ggsave("Cospeciation_DNA_order.pdf", height = 8, width=10)

#######

SV <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 1,col_names =T )
SV$name <- paste(SV$SV, SV$Genus, sep="_")
# Plot
SV$RNA <- SV$RNA*-1
SV_melt <-melt(SV, id.vars = c("name","SV", "Genus"), measure.vars = c("DNA", "RNA"))
# X Axis Breaks and Labels 
brks <- seq(-1, 1, 0.2)
lbls = as.character(c(seq(1, 0, -0.2), seq(0.2, 1, 0.2)))
tot_SV <-ggplot(SV_melt, aes(x = reorder(name, abs(value)), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(y="SNP heritability", x="", fill="") +
  theme_bw() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks.y  = element_blank(), axis.text.y = element_text(size=6)) +   # Centre plot title
  scale_fill_brewer(palette = "Dark2")  # Color palette
tot_SV
ggsave("SV_SNP_heritability.pdf", plot = tot_SV, height=10)


tot_SV2 <-ggplot(SV) +   # Fill column
  geom_bar(aes(x = reorder(SV, DNA), y = DNA, fill = CospeciationScoreGenus),stat = "identity", position="dodge", color="grey") + scale_fill_viridis_c(na.value="#F8F8FF")
tot_SV2 + coord_flip()+
  labs(y="SNP heritability", x="", fill="Cospeciation score genus") +
  theme_bw() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks.y  = element_blank())  # Color palette

ggsave("SV_SNP_heritability_cospec.pdf", height=8, width = 8)
# cospeciation
SV <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 1,col_names =T )
SV <- SV[,c(1:4, 6)]
tm_score_SV <-SV[complete.cases(SV),]


#Total DNA
SV_DNA <-ggplot(tm_score_SV, aes(x=`TM score`, y=DNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="TM score", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.1, size = 3)+ geom_text_repel(aes(label = SV), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.8, parse = TRUE, size = 3)+theme_bw() + ggtitle("DNA")
SV_DNA
SV_RNA <-ggplot(tm_score_SV, aes(x=`TM score`, y=RNA)) + geom_point()+
  geom_smooth(method = lm) + 
  labs(x="TM score", y="SNP heritability") +
  stat_fit_glance(method = 'lm',geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "" ), col="pink"),label.x  = 'middle', label.y = 0.3, size = 3)+ geom_text_repel(aes(label = SV), size = 3)+
  stat_poly_eq(aes(label = paste(..rr.label..), col="pink"), formula=y~x, label.x = "middle", label.y = 0.9, parse = TRUE, size = 3)+theme_bw() + ggtitle("RNA")
SV_RNA

p1<-plot_grid(SV_DNA,SV_RNA, nrow = 1, align = "v")
plot_grid(tot_SV, p1, ncol = 1,  axis="l", rel_heights = c(2, 1))
ggsave("TM_score_SV.pdf", height = 8, width=10)
SV_DNA
ggsave("TM_score_SV_DNA.pdf")

SV_RNA
ggsave("TM_score_SV_RNA.pdf")

# boxplot 
SV <- read_xlsx("../heritability_cospeciation.xlsx", sheet = 1,col_names =T )
SV$name <- paste(SV$SV, SV$Genus, sep="_")
melted <- melt(SV, measure.vars = c("DNA", "RNA"))
#
ggplot(melted,aes(x=reorder(Genus,CospeciationScoreGenus), y=value, fill=variable) ) + geom_bar(aes(y=CospeciationScoreGenus/4), stat="identity", position = "dodge")+geom_boxplot() + scale_fill_brewer(palette="Dark2")  + labs(x="", y="SNP heritability") + scale_y_continuous(sec.axis = sec_axis(~.*4, name = "Cospeciation score Genus")) + coord_flip()
ggsave("SV_boxplot_SNP_heritability_per_genus.pdf")
