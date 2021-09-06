if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(tidyverse)
library("AnnotationDbi")
library(stringr)
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/Genes/all_genes_snps/enrichment/")
# results_file<- read_csv2("sig.summaries/results_DNA_SW.csv")
# 
# genes <- results_file %>% 
#   select(all_genes) %>% 
#   distinct() %>% 
#   drop_na()
genes<- read_delim("../../all_marker_genes.txt", col_names =F, delim=" ")
genes<- genes$X1


# WIKIPATHWAYS #####
wp2gene <- read.gmt("~/Documents/PhD/Experiments/Final_QTL_mapping/wikipathways-20201010-gmt-Mus_musculus.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

# convert symbols back to entrez ids
#BiocManager::install("org.Mm.eg.db")
# BiocManager::install("org.Hs.eg.db")
library(org.Mm.eg.db) # remember to install it if you don't have it already
genes_entrez <- mapIds(org.Mm.eg.db, keys = str_to_title(genes), keytype = "SYMBOL", column="ENTREZID")
# genes_entrez<- transform(genes_entrez, gene=as.numeric(genes_entrez))
ewp <- enricher(genes_entrez, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
ewp <- setReadable(ewp, org.Mm.eg.db, keyType = "ENTREZID")
head(ewp)


cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))

cell_markers
y <- enricher(genes_entrez, TERM2GENE=cell_markers, minGSSize=1)
DT::datatable(as.data.frame(y))
y<- setReadable(y, org.Mm.eg.db, keyType = "ENTREZID")
DT::datatable(as.data.frame(y))

# MgSigDb analysis ####

library(msigdbr)
msigdbr_species()
m_df <- msigdbr(species = "Mus musculus")
head(m_df, 2) %>% as.data.frame
m_t2g <- msigdbr(species = "Mus musculus", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)

# HAllmark gene sets: are coherently expressed signatures derived 
# by aggregating many MSigDB gene sets to represent well-defined biological states or processes.
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
# curated gene sets  from online pathway databases, publications in PubMed, and knowledge of domain experts.
m_t2g <- msigdbr(species = "Mus musculus", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c2<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")
write_delim(as.data.frame(em_c2), "hallmark_curated_gene_sets.csv", delim=";")
# regulatory target gene sets  based on gene target predictions for microRNA seed sequences and predicted transcription factor binding sites.
m_t2g <- msigdbr(species = "Mus musculus", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c3<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")
write_delim(as.data.frame(em_c3), "hallmark_regulatory_target_gene_sets.csv", delim=";")
# computational gene sets  defined by mining large collections of cancer-oriented microarray data.
m_t2g <- msigdbr(species = "Mus musculus", category = "C4") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c4<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")

# ontology gene sets  consist of genes annotated by the same ontology term.
m_t2g <- msigdbr(species = "Mus musculus", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c5<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")
write_delim(as.data.frame(em_c5), "hallmark_ontology_gene_sets.csv", delim=";")
# oncogenic signature gene sets  defined directly from microarray gene expression data from cancer gene perturbations.

m_t2g <- msigdbr(species = "Mus musculus", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c6<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")

# immunologic signature gene sets  defined directly from microarray gene expression data from immunologic studies.
m_t2g <- msigdbr(species = "Mus musculus", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c6<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID") # nothing significant
write_delim(as.data.frame(em_c6), "hallmark_immunologic_sign_microarray.csv", delim=";")
# cell type signature gene sets  curated from cluster markers identified in single-cell sequencing studies of human tissue.
m_t2g <- msigdbr(species = "Mus musculus", category = "C8") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c8<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")
write_delim(as.data.frame(em_c8), "hallmark_cell_type_signature_sets_from_SCSseq.csv", delim=";")
# Disease analysis ####
# find human orthologs of mouse genes

# Basic function to convert mouse to human gene names
gene.df <- bitr(genes_entrez, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Mm.eg.db)

convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
humanx<- convertMouseGeneList(gene.df$SYMBOL)
library(org.Hs.eg.db)
genes_entrez_h <- mapIds(org.Hs.eg.db, keys = humanx, keytype = "SYMBOL", column="ENTREZID")


library(DOSE)
x <- enrichDO(gene          = genes_entrez_h,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE, 
              universe=NULL)
x <- setReadable(x, org.Hs.eg.db, keyType = "ENTREZID")
head(x)

# NEtwork of Cancer Gene ####
ncg <- enrichNCG(genes_entrez_h)
head(ncg)
# nothing

# DisGenNET ####
dgn <- enrichDGN(unname(genes_entrez_h))
dgn <- setReadable(dgn, org.Hs.eg.db, keyType = "ENTREZID")
library(enrichplot)
bar_dgn <-barplot(dgn, showCategory=20)
ggsave("disgennet_barplot.pdf")
# GO classifier #####
# back to mouse genes
ggo <- groupGO(gene     = as.character(genes_entrez),
               OrgDb    = org.Mm.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)
write_delim(as.data.frame(ggo),"GO_terms_CC.csv", delim=";")

# GO overrepresentation ####
ego <- enrichGO(gene          = as.character(genes_entrez),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
pdf("go_plot_BP.pdf")
goplot(ego2)
dev.off()
pdf("Cnetplot_GO_terms_BP.pdf")
cnetplot(ego2)
dev.off
p1 <- dotplot(ego2, showCategory=30) + ggtitle("dotplot for BP")
ego_mf <- enrichGO(gene          = genes_entrez,
                   OrgDb         = org.Mm.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
ego2_mf <- clusterProfiler::simplify(ego_mf, cutoff=0.7, by="p.adjust", select_fun=min)
pdf("cnetplot_MF_GO.pdf")
cnetplot(ego2_mf)
dev.off()
p2 <- dotplot(ego2_mf, showCategory=30) + ggtitle("dotplot for MF")

ego_cc <- enrichGO(gene          = genes_entrez,
                   OrgDb         = org.Mm.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
ego2_cc <- clusterProfiler::simplify(ego_cc, cutoff=0.7, by="p.adjust", select_fun=min)
pdf("cnetplot_GO_CC.pdf")
cnetplot(ego2_cc)
dev.off()
p3 <- dotplot(ego2_cc, showCategory=30) + ggtitle("dotplot for CC")
library(patchwork)
p1/p2/p3 +plot_layout(guides="collect")
ggsave("dotplot_GO_terms.pdf")
# KEGG analysis #####
search_kegg_organism('mmu', by='kegg_code')
kk <- enrichKEGG(gene         = as.character(genes_entrez),
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
head(kk)
kk<- setReadable(kk, org.Mm.eg.db, keyType = "ENTREZID")
write_delim(as.data.frame(kk), "kegg_pathways.csv", delim=";")
cnetplot(kk, showCategory = 14, colorEdge=T)
ggsave("netplot_kegg_pathways.pdf", width=15, height=12)
mkk <- enrichMKEGG(gene = as.character(genes_entrez),
                   organism = 'mmu')
head(mkk)
mkk<- setReadable(mkk, org.Mm.eg.db, keyType = "ENTREZID")
write_delim(as.data.frame(mkk), "mkegg_pathways.csv", delim=";")

library(pathview)

for (id in kk$ID){
  pathway <- pathview(gene.data  = as.character(genes_entrez),
                      pathway.id = id,
                      species    = "mmu")
}

library(meshes)
library(meshes)
x <- enrichMeSH(genes_entrez, MeSHDb = "MeSH.Mmu.eg.db", database='gendoo', category = 'C')
head(x)
x<- setReadable(x, org.Mm.eg.db, keyType = "ENTREZID")

write_delim(as.data.frame(x), "mesh.csv", delim=";")



edox <- setReadable(dgn, 'org.Mm.eg.db', 'ENTREZID')
p1 <- cnetplot(edox)
p1
ggsave("disease_genes.pdf")
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue")
p2
p3 <- cnetplot(edox, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

bar_dgn+p1
