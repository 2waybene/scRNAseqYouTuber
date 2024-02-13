## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7859841/
## enrichment with gprofiler2
library(DESeq2)
library(airway)
library(gprofiler2)


# load the airway data
data(airway)
# construct the DESeqDataSet object
ddsMat = DESeqDataSetFromMatrix(countData = assay(airway),
                                colData = colData(airway),
                                design = ~ cell + dex)
# run DESeq2 pipeline
dds = DESeq(ddsMat)
# get the results
results = results(dds, contrast = c("dex", "trt", "untrt"),
                  alpha = 0.05, lfcThreshold = 1)

# keep only the significant genes
results_sig = subset(results, padj < 0.05)
# get the significant up-regulated genes
up = subset(results_sig, log2FoldChange > 0)
# get the significant down-regulated genes
down = subset(results_sig, log2FoldChange < 0)
# enrichment analysis
gp_up = gost(row.names(up), organism = "hsapiens")
gp_down = gost(row.names(down), organism = "hsapiens") 


table(gp_up$result$source)


# GO:BP GO:CC  KEGG  REAC    TF    WP 
# 28     3     2     2     1     4 
# order genes by log2FC
up_ordered = up[order(up$log2FoldChange, decreasing = TRUE),]
# ordered enrichment analysis
gp_up_ordered = gost(row.names(up_ordered), organism = "hsapiens",
                     ordered_query = TRUE)
head(gp_up_ordered$result, 8)

table(gp_up_ordered$result$source)

#GO:BP GO:CC GO:MF  KEGG  REAC    TF    WP 
#27     3     6     1     2     1     2 

gostplot(gp_up, interactive = TRUE)

p1 = gostplot(gp_up, interactive = FALSE)
publish_gostplot(p1, highlight_terms = c("GO:0050896", "KEGG:04978",
                                         "REAC:R-HSA-5661231", "WP:WP3286"))


##. Analysing multiple gene lists

multi_gp = gost(list("up-regulated" = row.names(up),
                     "down-regulated" = row.names(down)))

p2 = gostplot(multi_gp, interactive = FALSE)
publish_gostplot(p2, highlight_terms = c("GO:0099699", "GO:0050896", "KEGG:04978",
                                         "REAC:R-HSA-5661231", "WP:WP3286",
                                         "GO:1990169"))
## convert gene id

results_genes = gconvert(row.names(results), organism = "hsapiens",
                         target = "ENTREZGENE_ACC", filter_na = FALSE)
head(results_genes)




