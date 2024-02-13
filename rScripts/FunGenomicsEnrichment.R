# https://www.youtube.com/watch?v=jqJ3xnTaGUU
# https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofiler/
  
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)

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