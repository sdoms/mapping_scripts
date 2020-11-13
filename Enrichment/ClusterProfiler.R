if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(tidyverse)
setwd("~/Documents/PhD/Experiments/Final_QTL_mapping/Results/Bacterial traits/DNA/")

results_file<- read.csv("summary_tables/genome_wide_significant/all_signif_snps/markers_with_genes_DNA.csv")
# 
# genes <- results_file %>% 
#   select(all_genes) %>% 
#   distinct() %>% 
#   drop_na()

genes <- na.omit(unique(as.vector(do.call('rbind', strsplit(as.character(results_file$all_genes),'|',fixed=TRUE))) ))
# WIKIPATHWAYS #####
wp2gene <- read.gmt("../../../wikipathways-20201010-gmt-Mus_musculus.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

# convert symbols back to entrez ids
library(org.Mm.eg.db) # remember to install it if you don't have it already
genes_entrez <- mapIds(org.Mm.eg.db, keys = genes, keytype = "SYMBOL", column="ENTREZID")

ewp <- enricher(genes_entrez, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
ewp <- setReadable(ewp, org.Mm.eg.db, keyType = "ENTREZID")
head(ewp)

# Cellmarker ####

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

# regulatory target gene sets  based on gene target predictions for microRNA seed sequences and predicted transcription factor binding sites.
m_t2g <- msigdbr(species = "Mus musculus", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c3<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")

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
#em_c6<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID") # nothing significant

# cell type signature gene sets  curated from cluster markers identified in single-cell sequencing studies of human tissue.
m_t2g <- msigdbr(species = "Mus musculus", category = "C8") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(genes_entrez, TERM2GENE=m_t2g)
head(em)
em_c8<- setReadable(em, org.Mm.eg.db, keyType = "ENTREZID")

# Disease analysis ####
library(DOSE)
x <- enrichDO(gene          = genes_entrez,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
x <- setReadable(x, org.Mm.eg.db, keyType = "ENTREZID")
head(x)

 # NEtwork of Cancer Gene ####
ncg <- enrichNCG(genes_entrez)
head(ncg)
# nothing

# DisGenNET ####
dgn <- enrichDGN(unname(genes_entrez))


# GO classifier #####
gene.df <- bitr(genes_entrez, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
ggo <- groupGO(gene     = genes_entrez,
               OrgDb    = org.Mm.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

# GO overrepresentation ####
ego <- enrichGO(gene          = genes_entrez,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
head(ego2)

# KEGG analysis #####
search_kegg_organism('mmu', by='kegg_code')
kk <- enrichKEGG(gene         = genes_entrez,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
head(kk)
kk<- setReadable(kk, org.Mm.eg.db, keyType = "ENTREZID")

mkk <- enrichMKEGG(gene = genes_entrez,
                   organism = 'mmu')
head(mkk)

library(meshes)
x <- enrichMeSH(genes_entrez, MeSHDb = "MeSH.Mmu.eg.db", database='gendoo', category = 'C')
head(x)
