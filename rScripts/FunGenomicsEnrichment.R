# https://www.youtube.com/watch?v=jqJ3xnTaGUU
# https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofiler/
  
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)



enrichment <- function (x, y){
  plot = enrichGO(
    x,
    org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = y,
    pvalueCutoff = 0.05,
    p.adjust.methods = "BH",
    qvalueCutoff = 0.2,
    minG
  )
}
enrichment.Mm <- function (x, y){
  plot = enrichGO(
    x,
    org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = y,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  dotplot(plot)
}
enrichment.Hs <- function (x, y){
  plot = enrichGO(
    x,
    org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = y,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
  )
  dotplot(plot)
}

data(geneList, package = "DOSE")
de <- names(geneList)[1:100]
yy <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
head(yy)
dotplot(yy)


enrichment.Mm (DEGs, "MF")

