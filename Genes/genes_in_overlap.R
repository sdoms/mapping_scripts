setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/overlap_with_other_studies/")
load("../../../all_genes_intervals/SW/genes_DNA_RNA_intervals_SW.Rdata")
library(ggplot2)
library(patchwork)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db) 
library(tidyverse)
library(readxl)
library(intervals)


# 
# reduced_intervals<- data.frame(matrix(ncol=3))
# colnames(reduced_intervals)<- c("chr", "start", "end")
# for (chr in 1:19){
#   x<- Intervals(as.matrix(input_SW[input_SW$chr==chr,c("start.LD.pos", "stop.LD.pos")]))
#   y<- reduce(x)
#   z<- data.frame(y)
#   colnames(z)<- c("start", "end")
#   z$chr<- paste0("chr",chr)
#   reduced_intervals<- rbind(reduced_intervals,z)
# }
# 
# # reduced_intervals<- reduced_intervals[-1,]
# 
# write.table(reduced_intervals, file="reduced_intervals_our_study_all.txt", row.names = F)
input_SW$length<- (input_SW$stop.LD.pos-input_SW$start.LD.pos)/1e6
input_SW<- transform(input_SW, length_num=as.numeric(length))
input_SW_unique_int <- input_SW %>% 
  group_by(chr) %>% 
  distinct(start.LD.pos, stop.LD.pos, length_num)

med_uniq_int<- median(input_SW_unique_int$length_num[input_SW_unique_int$length_num>0], na.rm=T)
mean_uniq_int<- mean(input_SW_unique_int$length_num[input_SW_unique_int$length_num>0], na.rm=T)

input_SW_10mb <- input_SW_unique_int[input_SW_unique_int$length_num<10,c("chr", "start.LD.pos", "stop.LD.pos")] 
input_SW_10mb<- input_SW_10mb %>% 
  arrange(chr, start.LD.pos, stop.LD.pos)
input_SW_10mb$chr<- paste0("chr", input_SW_10mb$chr)
write.table(input_SW_10mb, file="all_intervals_our_study_10mb.txt", row.names = F, col.names = F)


# reduced_intervals_10mb<- data.frame(matrix(ncol=3))
# colnames(reduced_intervals_10mb)<- c("chr", "start", "end")
# for (chr in 1:19){
#   x<- Intervals(as.matrix(input_SW[input_SW$chr==chr&input_SW$length_num<10,c("start.LD.pos", "stop.LD.pos")]))
#   y<- reduce(x)
#   z<- data.frame(y)
#   colnames(z)<- c("start", "end")
#   z$chr<- paste0("chr",chr)
#   reduced_intervals_10mb<- rbind(reduced_intervals_10mb,z)
# }
# 
# reduced_intervals_10mb<- reduced_intervals_10mb[-1,]
# 
# write.table(reduced_intervals_10mb, file="reduced_intervals_our_study_10mb.txt", row.names = F)

system("~/poverlap/poverlap2.py poverlap --a ../../other_studies/all_intervals_mouse_studies_10mb.txt --b all_intervals_our_study_10mb.txt -g ../../mouse.mm10.genome --n 9999")
system("bedtools intersect -wo -a ../../other_studies/all_intervals_mouse_studies_10mb.txt -b all_intervals_our_study_10mb.txt > bed_intersect_out.txt")

overlaps <- read_delim("bed_intersect_out.txt",delim = "\t", col_names = F)
colnames(overlaps)<- c("chrA", "startA", "endA","chrB", "startB", "endB", "width")
overlap_intervals <- overlaps %>% 
  rowwise() %>% 
  mutate(start=max(startA, startB), end=min(endA, endB)) %>% 
  mutate(length=end-start)

overlap_intervals_uniq <- overlap_intervals %>% 
  distinct(chrA, start, end)

target <- with(overlap_intervals_uniq,
               GRanges( seqnames = Rle(chrA),
                        ranges   = IRanges(start, end=end),
                        strand   = NULL ) )
loc <- VariantAnnotation::locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants())
loc
names(loc) <- NULL


out <- as.data.frame(loc)
# out<- out[, - c(10:11, 13:14)]
out<- out %>% 
  distinct(GENEID, QUERYID, .keep_all=T) %>%  drop_na(GENEID)

Symbol2id <- as.list( org.Mm.egSYMBOL2EG )
id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
names(id2Symbol) <- unlist(Symbol2id)

x <- unique( with(out, c(levels(as.factor(GENEID)))))
table( x %in% names(id2Symbol)) # good, all found

out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]


final_out <- out %>% 
  group_by(QUERYID) %>% 
  mutate(genes= paste0(unique(GENESYMBOL), collapse= ", ")) %>% 
  ungroup() %>% 
  distinct(seqnames, start, end, width, genes)
write.csv(final_out, "genes_in_overlaps_SW.csv")

genes_list_SW_in_overlaps<- out %>% 
  distinct(GENESYMBOL)

# loc_coding <- VariantAnnotation::locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, CodingVariants())
# 
# names(loc_coding) <- NULL
# 
# 
# out_coding <- as.data.frame(loc_coding)
# # out<- out[, - c(10:11, 13:14)]
# out_coding<- out_coding %>% 
#   distinct(GENEID, QUERYID, .keep_all=T) %>%  drop_na(GENEID)
# 
# x_coding <- unique( with(out_coding, c(levels(as.factor(GENEID)))))
# table( x_coding %in% names(id2Symbol)) # good, all found
# 
# out_coding$GENESYMBOL <- id2Symbol[ as.character(out_coding$GENEID) ]
# 
# 
# final_out <- out %>% 
#   group_by(QUERYID) %>% 
#   mutate(genes= paste0(unique(GENESYMBOL), collapse= ", ")) %>% 
#   ungroup() %>% 
#   distinct(seqnames, start, end, width, genes)
# write.csv(final_out, "genes_in_overlaps_SW.csv")
# 
# genes_list_SW_in_overlaps<- out %>% 
#   distinct(GENESYMBOL)

reduced_intervals<- data.frame(matrix(ncol=3))
colnames(reduced_intervals)<- c("chr", "start", "end")
for (chr in 1:19){
  x<- Intervals(as.matrix(overlap_intervals[overlap_intervals$chrA==paste0("chr",chr),c("start", "end")]))
  y<- reduce(x)
  z<- data.frame(y)
  colnames(z)<- c("start", "end")
  z$chr<- paste0("chr",chr)
  reduced_intervals<- rbind(reduced_intervals,z)
}

reduced_intervals<- reduced_intervals[-1,]
target <- with(reduced_intervals,
               GRanges( seqnames = Rle(chr),
                        ranges   = IRanges(start=start, end=end),
                        strand   = Rle(strand("*")) ) )
loc <- locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants())
loc
names(loc) <- NULL


out2 <- as.data.frame(loc)
out2<- out2[, - c(10:11, 13:14)]
out2<- unique(out2) %>%  drop_na(GENEID)

x <- unique( with(out2, c(levels(as.factor(GENEID)))))
table( x %in% names(id2Symbol)) # good, all found

out2$GENESYMBOL <- id2Symbol[ as.character(out2$GENEID) ]


final_out2 <- out2 %>% 
  group_by(QUERYID) %>% 
  mutate(genes= paste0(unique(GENESYMBOL), collapse= ", ")) %>% 
  ungroup() %>% 
  distinct(seqnames, start, end, width, genes)
write.csv(final_out2, "genes_in_overlap_reduced_SW.csv")
genes_reduced <- unique(out2$GENESYMBOL)
write.table(genes_reduced, file="genes_list_reduced-SW.txt", row.names = F, quote = F, col.names = F)

marker_genes<- read.csv("../../../../Shared/all_markers_RNA_DNA-with-genes_SW.csv", sep=",")
all_marker_genes<- stringr::str_to_upper(na.omit(unique(unlist(str_split(marker_genes$all_genes, " | ")))))
all_marker_genes<- c(all_marker_genes, stringr::str_to_upper(na.omit(unique(unlist(str_split(marker_genes$all_preceding_genes, " | "))))))
all_marker_genes<- c(all_marker_genes, stringr::str_to_upper(na.omit(unique(unlist(str_split(marker_genes$all_following_genes, " | "))))))
all_marker_genes<- unique(all_marker_genes)

library(data.table)
markers_in_overlap<-data.frame(matrix(ncol=ncol(marker_genes)))
colnames(markers_in_overlap)<- colnames(marker_genes)
markers_in_overlap<- markers_in_overlap[-1,]

for (chr in 1:19){
  spans<- reduced_intervals[reduced_intervals$chr==paste0("chr",chr),]
  markers_chr<- marker_genes[marker_genes$chr==chr,]
  markers_chr$pos<- markers_chr$pos * 1e6
  for (i in 1:nrow(spans)){
    
    aa<- markers_chr%>% filter(between(pos, spans[i,2], spans[i,3]))
    markers_in_overlap<- rbind(markers_in_overlap, aa)
  }
  
}
write_csv(markers_in_overlap, "markers_in_overlap_SW.csv")


# snp in genes#####
markers_in_overlap_genes<- stringr::str_to_upper(na.omit(unique(unlist(str_split(markers_in_overlap$all_genes, " | ")))))
# over/ under expression GF
OE_GF_COL <- read_excel("../mills/../../../overexpressed_protein_GF_con.xlsx", "GF", col_names = F)
colnames(OE_GF_COL)<- "genes"
# all closest overlapping genes in GF snp in gene
OE_GF_overlap_COL <- OE_GF_COL %>% 
  filter(genes %in% markers_in_overlap_genes)

OE_CON_COL <- read_excel("../mills/../../../overexpressed_protein_GF_con.xlsx", "con", col_names = F)
colnames(OE_CON_COL)<- "genes"
# all closest overlapping genes in CON snp in gene
OE_CON_overlap_COL <- OE_CON_COL %>% 
  filter(genes %in% markers_in_overlap_genes)

# preceding genes 
markers_in_overlap_preced_genes<- stringr::str_to_upper(na.omit(unique(unlist(str_split(markers_in_overlap$all_preceding_genes, " | ")))))
# over/ under expression GF
OE_GF_overlap_preced_COL <- OE_GF_COL %>% 
  filter(genes %in% markers_in_overlap_preced_genes)


OE_CON_overlap_preced_COL <- OE_CON_COL %>% 
  filter(genes %in% markers_in_overlap_preced_genes)


# following genes 
markers_in_overlap_follow_genes<- stringr::str_to_upper(na.omit(unique(unlist(str_split(markers_in_overlap$all_following_genes, " | ")))))
# over/ under expression GF
OE_GF_overlap_follow_COL <- OE_GF_COL %>% 
  filter(genes %in% markers_in_overlap_follow_genes)


OE_CON_overlap_follow_COL <- OE_CON_COL %>% 
  filter(genes %in% markers_in_overlap_follow_genes)

all_unique_marker_genes <- unique(c(markers_in_overlap_follow_genes,markers_in_overlap_genes,markers_in_overlap_preced_genes))

OE_CON_overlap_COL <- OE_CON_COL %>% 
  filter(genes %in% all_unique_marker_genes)
OE_GF_overlap_COL <- OE_GF_COL %>% 
  filter(genes %in% all_unique_marker_genes)

all_overlap_COL<- unique(c(OE_CON_overlap_COL$genes,OE_GF_overlap_COL$genes))

OE_CON_COL_snps <- OE_CON_COL %>% 
  filter(genes %in% all_marker_genes)
OE_GF_COL_snps <- OE_GF_COL %>% 
  filter(genes %in% all_marker_genes)
OE_COL_snps <- unique(c(OE_CON_COL_snps$genes, OE_GF_COL_snps$genes))

# all genes in intervals
# over/ under expression GF
OE_GF_overlap_intervals_COL <- OE_GF_COL %>% 
  filter(genes %in% stringr::str_to_upper(genes_reduced))

OE_CON_overlap_intervals_COL <- OE_CON_COL %>% 
  filter(genes %in% stringr::str_to_upper(genes_reduced))
OE_overlap_intervals_COL <- unique(c(OE_GF_overlap_intervals_COL$genes, OE_CON_overlap_intervals_COL$genes))
# write to files 
write_tsv(OE_GF_overlap_COL, "OE_GF_overlap_COL_SW.tsv")
write_tsv(OE_CON_overlap_COL, "OE_CON_overlap_COL_SW.tsv")
write_tsv(OE_GF_overlap_preced_COL, "OE_GF_overlap_preceding_genes_COL_SW.tsv")
write_tsv(OE_GF_overlap_follow_COL, "OE_GF_overlap_following_genes_COL_SW.tsv")
write_tsv(OE_CON_overlap_preced_COL, "OE_CON_overlap_preceding_genes_COL_SW.tsv")
write_tsv(OE_CON_overlap_follow_COL, "OE_CON_overlap_following_genes_COL_SW.tsv")


write_tsv(OE_GF_overlap_intervals_COL, "OE_GF_overlap_interval_genes_COL_SW.tsv")
write_tsv(OE_CON_overlap_intervals_COL, "OE_CON_overlap_interval_genes_COL_SW.tsv")

# 
# # TOTAL (not just overlap) #####
# load("/Users/doms/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/genes_DNA_RNA.Rdata")
# all_genes<- unique(stringr::str_to_upper(na.omit(unique(unlist(str_split(genes_DNA_RNA$list_total_genes, ", "))))))
# all_genes_COLON_GF<- OE_GF_COL %>% 
#   filter(genes %in% stringr::str_to_upper(all_genes))
# all_genes_COLON_CON<- OE_CON_COL %>% 
#   filter(genes %in% stringr::str_to_upper(all_genes))
# write_tsv(all_genes_COLON_CON, "../Results/Bacterial traits/protein_overexpression/all_genes_CON_COLON.txt")
# write_tsv(all_genes_COLON_GF, "../Results/Bacterial traits/protein_overexpression/all_genes_GF_COLON.txt")
# 
# all_genes_COLON <- unique(c(all_genes_COLON_CON$genes, all_genes_COLON_GF$genes))
# 
# all_genes_markers_COLON_GF<- OE_GF_COL%>% 
#   filter(genes %in% stringr::str_to_upper(marker_genes$all_genes))
# all_genes_markers_COLON_CON<- OE_CON_COL%>% 
#   filter(genes %in% stringr::str_to_upper(marker_genes$all_genes))
# 
# all_genes_markers_preced_COLON_GF<- OE_GF_COL%>% 
#   filter(genes %in% stringr::str_to_upper(marker_genes$all_preceding_genes))
# all_genes_markers_preced_COLON_CON<- OE_CON_COL %>% 
#   filter(genes %in% stringr::str_to_upper(marker_genes$all_preceding_genes))
# 
# all_genes_markers_follow_COLON_GF<- OE_GF_COL%>% 
#   filter(genes %in% stringr::str_to_upper(marker_genes$all_following_genes))
# all_genes_markers_follow_COLON_CON<- OE_CON_COL%>% 
#   filter(genes %in% stringr::str_to_upper(marker_genes$all_following_genes))
# 
# write_tsv(all_genes_markers_COLON_CON, "../Results/Bacterial traits/protein_overexpression/all_genes_markers_CON_COLON.txt")
# write_tsv(all_genes_markers_COLON_GF, "../Results/Bacterial traits/protein_overexpression/all_genes_markers_GF_COLON.txt")
# write_tsv(all_genes_markers_COLON_CON, "../Results/Bacterial traits/protein_overexpression/all_genes_markers_CON_COLON.txt")
# write_tsv(all_genes_markers_COLON_GF, "../Results/Bacterial traits/protein_overexpression/all_genes_markers_GF_COLON.txt")
# write_tsv(all_genes_markers_COLON_CON, "../Results/Bacterial traits/protein_overexpression/all_genes_markers_CON_COLON.txt")
# write_tsv(all_genes_markers_COLON_GF, "../Results/Bacterial traits/protein_overexpression/all_genes_markers_GF_COLON.txt")

########### small intestine ###########

OE_GF_SI <- read_excel("../mills/../../../overexpressed_protein_GF_con.xlsx", "GF_SI", col_names = F)
colnames(OE_GF_SI)<- c("genes", "score")
OE_GF_overlap_SI <- OE_GF_SI %>% 
  filter(genes %in% markers_in_overlap_genes)
OE_GF_overlap_follow_SI <- OE_GF_SI %>% 
  filter(genes %in% markers_in_overlap_follow_genes)
OE_GF_overlap_preced_SI <- OE_GF_SI %>% 
  filter(genes %in% markers_in_overlap_preced_genes)

OE_CON_SI <- read_excel("../mills/../../../overexpressed_protein_GF_con.xlsx", "CON_SI", col_names = F)
colnames(OE_CON_SI)<- c("genes", "score")
OE_CON_overlap_SI <- OE_CON_SI %>% 
  filter(genes %in% markers_in_overlap_genes)
OE_CON_overlap_follow_SI <- OE_CON_SI %>% 
  filter(genes %in% markers_in_overlap_follow_genes)
OE_CON_overlap_preced_SI <- OE_CON_SI %>% 
  filter(genes %in% markers_in_overlap_preced_genes)

OE_CON_overlap_SI_snps <- OE_CON_SI %>% 
  filter(genes %in% all_unique_marker_genes)
OE_GF_overlap_SI_snps <- OE_GF_SI %>% 
  filter(genes %in% all_unique_marker_genes)
OE_overlap_SI_snps <- unique(c(OE_CON_overlap_SI_snps$genes, OE_GF_overlap_SI_snps$genes))

OE_GF_overlap_intervals_SI <- OE_GF_SI %>% 
  filter(genes %in% stringr::str_to_upper(genes_reduced))

OE_CON_overlap_intervals_SI <- OE_CON_SI %>% 
  filter(genes %in% stringr::str_to_upper(genes_reduced))

OE_overlap_intervals_SI <- unique(c(OE_GF_overlap_intervals_SI$genes, OE_CON_overlap_intervals_SI$genes))

write_tsv(OE_GF_overlap_SI, "OE_GF_overlap_SI.tsv")
write_tsv(OE_CON_overlap_SI, "OE_CON_overlap_SI.tsv")
write_tsv(OE_GF_overlap_preced_SI, "OE_GF_overlap_preceding_genes_SI.tsv")
write_tsv(OE_GF_overlap_follow_SI, "OE_GF_overlap_following_genes_SI.tsv")
write_tsv(OE_CON_overlap_preced_SI, "OE_CON_overlap_preceding_genes_SI.tsv")
write_tsv(OE_CON_overlap_follow_SI, "OE_CON_overlap_following_genes_SI.tsv")


write_tsv(OE_GF_overlap_intervals_SI, "OE_GF_overlap_interval_genes_SI.tsv")
write_tsv(OE_CON_overlap_intervals_SI, "OE_CON_overlap_interval_genes_SI.tsv")

write(OE_overlap_intervals_SI, "OE_overlap_interval_genes_SI.txt")

write(OE_overlap_intervals_COL, "OE_overlap_interval_genes_COL.txt")


# snp in genes
# all_genes<- stringr::str_to_upper(na.omit(unique(unlist(str_split(genes_DNA_RNA$list_total_genes, ", ")))))
all_genes_SI_GF<- OE_GF_SI %>% 
  filter(genes %in% stringr::str_to_upper(all_genes))
all_genes_SI_CON<- OE_CON_SI %>% 
  filter(genes %in% stringr::str_to_upper(all_genes))
write_tsv(all_genes_SI_CON, "./snps/all_genes_CON_SI.txt")
write_tsv(all_genes_SI_GF, "./snps/all_genes_GF_SI.txt")

# all genes shared: mine + OE
all_genes_SI <- unique(c(all_genes_SI_CON$genes, all_genes_SI_GF$genes))

all_genes_markers_SI_GF<- OE_GF_SI %>% 
  filter(genes %in% stringr::str_to_upper(marker_genes$all_genes))
all_genes_markers_SI_CON<- OE_CON_SI %>% 
  filter(genes %in% stringr::str_to_upper(marker_genes$all_genes))


all_genes_markers_preced_SI_GF<- OE_GF_SI %>% 
  filter(genes %in% stringr::str_to_upper(marker_genes$all_preceding_genes))
all_genes_markers_preced_SI_CON<- OE_CON_SI %>% 
  filter(genes %in% stringr::str_to_upper(marker_genes$all_preceding_genes))

all_genes_markers_follow_SI_GF<- OE_GF_SI %>% 
  filter(genes %in% stringr::str_to_upper(marker_genes$all_following_genes))
all_genes_markers_follow_SI_CON<- OE_CON_SI %>% 
  filter(genes %in% stringr::str_to_upper(marker_genes$all_following_genes))


OE_CON_SI_snps <- OE_CON_SI %>% 
  filter(genes %in% all_marker_genes)
OE_GF_SI_snps <- OE_GF_SI %>% 
  filter(genes %in% all_marker_genes)
OE_SI_snps <- unique(c(OE_CON_SI_snps$genes, OE_GF_SI_snps$genes))


write_tsv(all_genes_markers_SI_CON, "snps/all_genes_markers_CON_SI.txt")
write_tsv(all_genes_markers_SI_GF, "snps/all_genes_markers_GF_SI.txt")
write_tsv(all_genes_markers_SI_CON, "snps/all_genes_markers_CON_SI.txt")
write_tsv(all_genes_markers_SI_GF, "snps/all_genes_markers_GF_SI.txt")
write_tsv(all_genes_markers_SI_CON, "snps/all_genes_markers_CON_SI.txt")
write_tsv(all_genes_markers_SI_GF, "snps/all_genes_markers_GF_SI.txt")

######## TOTAL #####
all_OE <- unique(c(OE_GF_SI$genes, OE_CON_SI$genes, OE_CON_COL$genes, OE_GF_COL$genes))

# intervals 
OE_overlap_intervals_GF <- unique(c(OE_GF_overlap_intervals_COL$genes, OE_GF_overlap_intervals_SI$genes))
OE_overlap_intervals_CON <- unique(c(OE_CON_overlap_intervals_COL$genes, OE_CON_overlap_intervals_SI$genes))

OE_overlap_intervals_all <- unique(c(OE_overlap_intervals_SI, OE_overlap_intervals_COL))
# all_genes_OE <- unique(c(all_genes_COL, all_genes_SI))

# snps 
OE_snps <- unique(c(OE_COL_snps, OE_SI_snps))
OE_overlap_snps <- unique(c(OE_overlap_SI_snps, all_overlap_COL))

write.table(OE_overlap_intervals_all, "Oe_overlap_intervals_all.txt", row.names = F, col.names = F, quote = F)
write.table(OE_overlap_snps, "Oe_overlap_snps.txt", row.names = F, col.names = F, quote = F)

# result<- data.frame(matrix(ncol=ncol(final_out2)))
# colnames(result)<- colnames(final_out2)
# result$gene <- NA
# for (gene in OE){
#   res <- final_out2[grepl(gene, list_all_genes),]
#   res$gene <- gene
#   result <- rbind(result, res)
# }
# 
# result$start <- result$start/ 1e6
# result$end <-result$end/ 1e6
# result$width<- result$width/1e6
write(OE_overlap_intervals_GF, "OE_overlap_intervasls_GF.txt")
write(OE_overlap_intervals_CON, "OE_overlap_intervasls_CON.txt")
write(all_unique_marker_genes, "all_genes_from_markers_in_overlap.txt")



###### Venn diagram #####
library(VennDiagram)
library(stringr)


# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

all_other_mouse_studies <- read.table("../../other_studies/all_intervals_mouse_studies_10mb.txt")
colnames(all_other_mouse_studies)<- c("chr", "start", "end")
target <- with(all_other_mouse_studies,
               GRanges( seqnames = Rle(chr),
                        ranges   = IRanges(start, end=end),
                        strand   = NULL ) )
loc <- VariantAnnotation::locateVariants(target, TxDb.Mmusculus.UCSC.mm10.knownGene, AllVariants())
loc
names(loc) <- NULL


out3 <- as.data.frame(loc)
# out<- out[, - c(10:11, 13:14)]
out3<- out3 %>% 
  distinct(GENEID, QUERYID, .keep_all=T) %>%  drop_na(GENEID)

Symbol2id <- as.list( org.Mm.egSYMBOL2EG )
id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
names(id2Symbol) <- unlist(Symbol2id)

x <- unique( with(out3, c(levels(as.factor(GENEID)))))
table( x %in% names(id2Symbol)) # good, all found

out3$GENESYMBOL <- id2Symbol[ as.character(out3$GENEID) ]


final_out3 <- out3 %>% 
  group_by(QUERYID) %>% 
  mutate(genes= paste0(unique(GENESYMBOL), collapse= ", ")) %>% 
  ungroup() %>% 
  distinct(seqnames, start, end, width, genes)
write.csv(final_out3, "genes_in_other_studies.csv")

genes_other_studies<- out3 %>% 
  distinct(GENESYMBOL)

# Chart
venn.diagram(
  x = list(na.omit(all_genes),na.omit(all_marker_genes), na.omit(all_OE), 
           na.omit(str_to_upper(genes_other_studies$GENESYMBOL))),
  category.names = c("Present study" , "Closest genes to SNP" , "Mills et al., 2020", "Previous mouse studies"),
  filename = 'vendiagram_overlap.png',
  output=TRUE,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff', "green"),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3), alpha("green", 0.3)),
  cex = 0.6,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.prompts=T,
  cat.pos = c(45, 90, 135, 180),
  cat.dist = c(0.055, 0.055, 0.085, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#f7cb44ff', "dark green"),
  rotation = 1

)

venn.plot <- venn.diagram(
  x = list(
    "Genes intervals" = na.omit(all_genes),
    "Genes SNPs" = na.omit(all_marker_genes),
    "Mills, et al., 2020" =  na.omit(all_OE),
    "Previous studies" = na.omit(str_to_upper(genes_list_SW_in_overlaps$GENESYMBOL))
  ),
  filename = "Venn_4set_pretty.tiff",
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", 
                "white", "white", "white", "white", "darkblue", "white", 
                "white", "white", "white", "darkgreen", "white"),
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.07,
  cat.fontfamily = "serif",
  rotation.degree = 270,
  margin = 0.2
)

ss<- calculate.overlap(x = list(
  "Genes intervals" = na.omit(all_genes),
  "Genes SNPs" = na.omit(all_marker_genes),
  "Mills, et al., 2020" =  na.omit(all_OE),
  "Previous studies" = na.omit(str_to_upper(genes_other_studies$GENESYMBOL)), 
  "Overlap"= str_to_upper(genes_list_SW_in_overlaps$GENESYMBOL)))

  