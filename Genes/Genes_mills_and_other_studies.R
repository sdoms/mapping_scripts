setwd("~/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/genes_overlap/SW/")
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_18032021_SW.Rdata")
load("/Users/doms/Documents/Research/Experiments/Final_QTL_mapping/Results/Bacterial traits/all_DNA_RNA_22022022_SW_pval_corr.Rdata")

library(ggplot2)
library(patchwork)
library(VariantAnnotation)
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # for annotation
library(org.Mm.eg.db) 
library(tidyverse)
library(readxl)
library(intervals)
input_SW<- all_dna_rna_SW
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
write.table(input_SW_10mb, file="all_intervals_our_study_10mb_pval_corr.txt", row.names = F, col.names = F)

reduced_intervals_10mb<- data.frame(matrix(ncol=3))
colnames(reduced_intervals_10mb)<- c("chr", "start", "end")
for (chr in 1:19){
  x<- Intervals(as.matrix(input_SW[input_SW$chr==chr&input_SW$length_num<10,c("start.LD.pos", "stop.LD.pos")]))
  y<- reduce(x)
  z<- data.frame(y)
  colnames(z)<- c("start", "end")
  z$chr<- paste0("chr",chr)
  reduced_intervals_10mb<- rbind(reduced_intervals_10mb,z)
}

reduced_intervals_10mb<- reduced_intervals_10mb[-1,]

write.table(reduced_intervals_10mb, file="reduced_intervals_our_study_10mb_pval_corr.txt", row.names = F)



system("~/poverlap/poverlap2.py poverlap --a ../other_studies/all_intervals_mouse_studies_10mb.txt --b all_intervals_our_study_10mb_pval_corr.txt -g ../mouse.mm10.genome --n 9999")
system("bedtools intersect -wo -a ../other_studies/all_intervals_mouse_studies_10mb.txt -b all_intervals_our_study_10mb_pval_corr.txt > bed_intersect_out_pval_corr.txt")
system("~/poverlap/poverlap2.py poverlap --a ../other_studies/all_intervals_mouse_studies_10mb.txt --b reduced_intervals_our_study_10mb_pval_corr.txt -g ../mouse.mm10.genome --n 9999")
system("bedtools intersect -wo -a ../other_studies/all_intervals_mouse_studies_10mb.txt -b reduced_intervals_our_study_10mb_pval_corr.txt > bed_reduced_intersect_out_pval_corr.txt")

overlaps <- read_delim("bed_intersect_out_pval_corr.txt",delim = "\t", col_names = F)
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
write_delim(final_out, "intervals_in_other_studies_with_genes_SW_pval_corr.csv", delim=";")

genes_list_SW_in_overlaps<- out %>% 
  distinct(GENESYMBOL)

write.table(genes_list_SW_in_overlaps, file="genes_intervals_in_other_studies-SW_pval_corr.txt", row.names = F, quote = F, col.names = F)

#### All closest genes to sig SNPs ####
marker_genes<- read.csv("../../all_genes_snps/all_markers_RNA_DNA-with-genes_SW_pval_corr.csv", sep=",")
all_marker_genes<- stringr::str_to_upper(na.omit(unique(unlist(str_split(marker_genes$all_genes, " | ")))))
all_marker_genes<- c(all_marker_genes, stringr::str_to_upper(na.omit(unique(unlist(str_split(marker_genes$all_preceding_genes, " | "))))))
all_marker_genes<- c(all_marker_genes, stringr::str_to_upper(na.omit(unique(unlist(str_split(marker_genes$all_following_genes, " | "))))))
all_marker_genes<- unique(all_marker_genes)
all_marker_genes<- all_marker_genes[-2] # remove NA
all_marker_genes<- all_marker_genes[-2] # remove |
write(all_marker_genes, "../../all_marker_genes_pval_corr.txt")

### which of our significant SNPS are in other studies? ####
library(data.table)
markers_in_overlap<-data.frame(matrix(ncol=ncol(marker_genes)))
colnames(markers_in_overlap)<- colnames(marker_genes)
markers_in_overlap<- markers_in_overlap[-1,] # remove NA row

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
for (chr in 1:19){
  spans<- reduced_intervals[reduced_intervals$chr==paste0("chr",chr),]
  markers_chr<- marker_genes[marker_genes$chr==chr,]
  markers_chr$pos<- markers_chr$pos
  for (i in 1:nrow(spans)){
    
    aa<- markers_chr%>% filter(between(pos, spans[i,2], spans[i,3]))
    markers_in_overlap<- rbind(markers_in_overlap, aa)
  }
  
}
write_csv(markers_in_overlap, "markers_in_other_studies_SW_pval_corr.csv")
# which closest genes to SNPs are in other studies ####
marker_genes_in_other_studies <- stringr::str_to_upper(na.omit(unique(unlist(str_split(markers_in_overlap$all_genes, " | ")))))
marker_genes_in_other_studies<- c(marker_genes_in_other_studies, stringr::str_to_upper(na.omit(unique(unlist(str_split(markers_in_overlap$all_preceding_genes, " | "))))))
marker_genes_in_other_studies<- c(marker_genes_in_other_studies, stringr::str_to_upper(na.omit(unique(unlist(str_split(markers_in_overlap$all_following_genes, " | "))))))
marker_genes_in_other_studies<- unique(marker_genes_in_other_studies)
marker_genes_in_other_studies <- marker_genes_in_other_studies[-3] # remove |
marker_genes_in_other_studies <- marker_genes_in_other_studies[-3] # remove NA
write(marker_genes_in_other_studies, "marker_genes_in_other_studies.txt")
# which genes from intervals are in Mills ####
all_genes<- unique(stringr::str_to_upper(na.omit(unique(unlist(str_split(input_SW$list_total_genes, ", "))))))

all_mills_genes <- read_excel("../../overexpressed_protein_GF_con.xlsx", "ALL_DF", col_names = F)
colnames(all_mills_genes)<- "mills"
genes_intervals_in_mills <- all_mills_genes %>% 
  filter(mills %in% stringr::str_to_upper(all_genes))
write(genes_intervals_in_mills$mills, "genes_intervals_in_mills_pval_corr.txt")
# which genes from snps are in Mills ####
closest_genes_snps_in_mills <- all_mills_genes %>% 
 filter(mills %in% stringr::str_to_upper(all_marker_genes))
write(closest_genes_snps_in_mills$mills, "closest_genes_snps_in_mills_pval_corr.txt")
# which genes from intervals are in other studies and Mills ####
genes_intervals_in_mills_and_other_studies <- all_mills_genes %>% 
  filter(mills %in% stringr::str_to_upper(genes_list_SW_in_overlaps$GENESYMBOL)) %>% drop_na() 
write(genes_intervals_in_mills_and_other_studies$mills, "genes_intervals_in_mills_and_other_studies_pval_corr.txt")
# which closest genes snps are in other studies and Mills ####
closest_genes_snps_mills_and_other_studies <- all_mills_genes %>% 
  filter(mills %in% stringr::str_to_upper(marker_genes_in_other_studies))
write(closest_genes_snps_mills_and_other_studies$mills, "closest_genes_snps_in_mills_and_other_studies_pval_corr.txt")
