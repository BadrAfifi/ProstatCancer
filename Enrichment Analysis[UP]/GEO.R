library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(cowplot)
library(ggnewscale)
require(pathview)
library(KEGGREST)
library(topGO)

BiocManager::install("KEGGREST",force = TRUE)


uni <- readLines("DEGsuniverse.txt")
gene_names <-  read.table("Overlapped_up_DEGs.txt")
names(gene_names)[1] = "SYMBOL"


hs <- org.Hs.eg.db
my.symbols <- gene_names$SYMBOL
IDs<- AnnotationDbi::select(hs, 
       keys = my.symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

geneID <- IDs$ENTREZID
length(geneID)
##  decreasing order
geneID = sort(geneID, decreasing = TRUE)

write(gene_names$up_names,"genes.txt")

#enrich : over-representation analysis (ORA)
# gse : gene set enrichment analysis (GSEA)


{
# Univarse genes
uni <- readLines("DEGsuniverse.txt")
uniIDs<- AnnotationDbi::select(hs, 
                            keys = uni,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")

uni <- uniIDs$ENTREZID
length(uni)
##  decreasing order
uni = sort(uni, decreasing = TRUE)
}

#        over-representation analysis

#==============[ Bp ]====================================

ego_CC <- enrichGO(gene = geneID,
                   keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db,
                ont = "CC",
                pAdjustMethod= "fdr",
                qvalueCutoff = 0.05,
                minGSSize = 10,
                maxGSSize = 500,
                readable = T)

# [1] Bar Plot
barplot(ego_CC, showCategory = 10)
# [2] Dot Plot
dotcc <- dotplot(ego_CC, showCategory = 10)
# [3] Encrichment map:
ego_CCx2 <- pairwise_termsim(ego_CC)
em_CC <- emapplot(ego_CCx2, showCategory = 10)
# [4] Category Netplot:
circ_CC <- cnetplot(ego_CC,circular = TRUE, colorEdge = TRUE)


# [5] Selecting high count sets
head(summary(ego_CC))
top.CC <- as.data.frame(ego_CC)
top.CC<- subset(top.CC, Count>10)
top.CC[1:10, 1:6]

ggplot(top.CC, aes(y = Description, x = GeneRatio, color= p.adjust , size =Count)) +
  geom_point(aes(color = p.adjust), alpha = 1.0) +
  scale_colour_gradient(low="blue", high="red")+
  theme(panel.grid.major = element_line(color = "steelblue4", size = 0.5, linetype = "dotted",),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect( fill = NA) #linetype = "dashed"
  ) 


#==============[ Bp ]====================================

ego_Bp <- enrichGO(gene = geneID,
                   keyType = "ENTREZID",
                   ont = "BP",
                   OrgDb = org.Hs.eg.db,
                   pAdjustMethod= "fdr",
                   qvalueCutoff = 0.05,
                   minGSSize = 10,
                   maxGSSize = 500,
                   readable = T)

barplot(ego_Bp, showCategory = 10)
dotBp <- dotplot(ego_Bp, showCategory = 10)
ego_Bpx2 <- pairwise_termsim(ego_Bp)
em_BP <- emapplot(ego_Bpx2, showCategory = 10)
# [4] Category Netplot:
circ_Bp <- cnetplot(ego_Bp,circular = TRUE, colorEdge = TRUE)

summary(ego_Bp)[1:10,1:8]

# [5] Selecting high count sets
top.Bp <- as.data.frame(ego_Bp)
top.Bp<- subset(top.Bp, Count>10)
top.Bp[1:5, 1:6]

ggplot(top.Bp, aes(y = Description, x = GeneRatio, color= p.adjust , size =Count)) +
  geom_point(aes(color = p.adjust), alpha = 1.0) +
  scale_colour_gradient(low="blue", high="red")+
  theme(panel.grid.major = element_line(color = "steelblue4", size = 0.5, linetype = "dotted",),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect( fill = NA) #linetype = "dashed"
  ) 

#=================================[ MF ]===================================================

ego_MF <- enrichGO(gene = geneID,
                   keyType = "ENTREZID", #ENSEMBL
                   ont = "MF" ,
                   OrgDb = org.Hs.eg.db,
                   pAdjustMethod= "fdr",
                   qvalueCutoff = 0.05,
                   minGSSize = 10,
                   maxGSSize = 500,
                   readable = T)

barplot(ego_MF, showCategory = 10)
dotMF <- dotplot(ego_MF, showCategory = 10)
ego_MFx2 <- pairwise_termsim(ego_MF)
em_MF <- emapplot(ego_MFx2, showCategory = 10)
# [4] Category Netplot:
circ_MF <- cnetplot(ego_MF,circular = TRUE, colorEdge = TRUE)

summary(ego_MF)


# [5] Selecting high count sets
top.MF <- as.data.frame(ego_MF)
top.MF <- subset(top.MF, Count>5)
top.MF [1:10, 1:6]

ggplot(top.MF, aes(y = Description, x = GeneRatio, color= p.adjust , size =Count)) +
  geom_point(aes(color = p.adjust), alpha = 1.0) +
  scale_colour_gradient(low="blue", high="red")+
  theme(panel.grid.major = element_line(color = "steelblue4", size = 0.5, linetype = "dotted",),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect( fill = NA) #linetype = "dashed"
  ) 
#===================[ GO ALL ]===========================================

cowplot::plot_grid(dotBp, dotcc, dotMF, ncol=2, labels=LETTERS[1:3])
cowplot::plot_grid(em_BP, em_CC, em_MF, ncol=2, labels=LETTERS[1:3])
cowplot::plot_grid(circ_Bp, circ_CC, circ_MF, ncol=2, labels=LETTERS[1:3]) # [4] Category Netplot: cnetplot(ego_MF)


#===================[ GO ALL ]===========================================

# Get the GO annotations for the DEGs
eg2go_hs <- enrichGO(gene          = geneID,
                     keyType = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "ALL",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff = 0.05,
                     readable      = TRUE,
                     minGSSize = 10,
                     maxGSSize = 500,)

barplot(eg2go_hs, showCategory = 10)
dotplot(eg2go_hs, showCategory = 10)
gox2 <- pairwise_termsim(eg2go_hs)
emapplot(gox2, showCategory = 10)
# [4] Category Netplot:
GOALL <- cnetplot(eg2go_hs, categorySize = "pvalue")
#[6]
cnetplot(eg2go_hs, circular = TRUE, colorEdge = TRUE) 

summary(eg2go_hs)[1:10,1:8]

# [5] Selecting high count sets
top.egoAll <- as.data.frame(eg2go_hs)
summary(top.egoAll$Count)
dim(top.egoAll)
top.egoAll <- subset(top.egoAll, Count>10)
summary(top.egoAll$Count)
dim(top.egoAll)
top.egoAll [1:5, 1:6]

ggplot(top.egoAll, aes(y = Description, x = ONTOLOGY, color= p.adjust , size =Count)) +
  geom_point(aes(color = p.adjust), alpha = 1.0) +
  scale_colour_gradient(low="blue", high="red")+
  theme(panel.grid.major = element_line(color = "steelblue4", size = 0.5, linetype = "dotted",),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect( fill = NA) #linetype = "dashed"
  ) 

#================[KEGG Pathways]=========================================


kk <- enrichKEGG(gene = geneID,
                 organism = "hsa",
                 keyType = "kegg",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr",
                 minGSSize = 10,
                 maxGSSize = 500,
                 qvalueCutoff = 0.05,
                 use_internal_data = T)
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(summary(kk))
pairwise_termsim(kk)
#[1]
barplot(kk, showCategory = 10)
#[2]
dotplot(kk, showCategory = 10)
#[3]
x2 <- pairwise_termsim(kk)
emapplot(x2)
emapplot(kk,showCategory = 10)
#[4]
KEGG <- cnetplot(kk, categorySize = "geneNum")
# [5]
cnetplot(kk2, foldChange= geness, circular = TRUE, colorEdge = TRUE) 



# illustrate how to visualize “hsa04110” pathway, which was enriched in KEGG.
pathview(gene.data  = geneID,pathway.id = "hsa05146", species    = "hsa")

#================================================================

cowplot::plot_grid(GOALL, KEGG, ncol=2, labels=LETTERS[1:2])

#===============[gene set enrichment analysis]=================


kk3 <- gseKEGG(gene    = geneID,
               keyType = "kegg",
               organism     = 'hsa',
               pAdjustMethod = "BH",
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk)
gseaplot(kk2, geneSetID = "hsa05146")

#=====================================================================
gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb        = org.Hs.eg.db, 
             pAdjustMethod = "none")
gseaplot(kk2, geneSetID = "hsa05146")
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


remove(kk2)


BiocManager::install("hgu95av2.db")
library(hgu95av2.db)
library(ALL)
library(GOplot) # Load the dataset 
data(ALL)
data(geneList)
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)
sum(topDiffGenes(geneList))
sampleGOdata <- new("topGOdata",
                    + "Simple session", ontology = "BP",
                    + geneList, geneSel = topDiffGenes,
                    + 10,
                    + annFUN.db, affyLib = affyLib)
