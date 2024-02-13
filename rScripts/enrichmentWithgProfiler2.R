##=========================================================
##  This is sample code in the original paper
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

head(gp_up$result)
dim(up)
dim(gp_up$result)
table(gp_up$result$source)


# GO:BP GO:CC  KEGG  REAC    TF    WP 
# 28     3     2     2     1     4 
# order genes by log2FC
up_ordered = up[order(up$log2FoldChange, decreasing = TRUE),]
# ordered enrichment analysis
gp_up_ordered = gost(row.names(up_ordered), organism = "hsapiens",
                     ordered_query = TRUE)

dim(gp_up_ordered$result)
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


results_df = as.data.frame(results)
results_df$Ensembl_id = row.names(results_df)
results_df = results_df[order(results_df$padj),]

# add the gene names
results_df = merge(results_df,
                   results_genes[,c("input", "target", "name", "description")],
                   by.x = "Ensembl_id", by.y = "input")

# save the results to a tsv file
write.table(results_df, file = "/Users/li11/myGit/scRNAseqYouTuber/sampleAnalysisResults/DESeq2_results.tsv", sep = "\t",
            quote = F, row.names = F)


##.  Integrating with external tools for visualisations

library(clusterProfiler)
library(enrichplot)
library(DOSE) # needed to convert to enrichResult object

up_names = gconvert(row.names(up))
down_names = gconvert(row.names(down))

# enrichment analysis using gene names
multi_gp = gost(list("up-regulated" = up_names$name,
                     "down-regulated" = down_names$name), multi_query = FALSE, evcodes = TRUE)


# modify the g:Profiler data frame
gp_mod = multi_gp$result[,c("query", "source", "term_id",
                            "term_name", "p_value", "query_size", 
                            "intersection_size", "term_size", 
                            "effective_domain_size", "intersection")]

gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                  "query_size", "Count", "term_size", "effective_domain_size", 
                  "geneID", "GeneRatio", "BgRatio")
gp_mod$geneID = gsub(",", "/", gp_mod$geneID)

row.names(gp_mod) = paste0(gp_mod$ID, "", gp_mod$Cluster)

# define as compareClusterResult object
gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)

# define as enrichResult object
gp_mod_enrich  = new("enrichResult", result = gp_mod)
enrichplot::dotplot(gp_mod_cluster)


barplot(gp_mod_enrich, showCategory = 40, font.size = 6) + 
  ggplot2::facet_grid(~Cluster) +
  ggplot2::ylab("Intersection size")


gostres = gost(query = list("up-regulated" = row.names(up)),
               evcodes = TRUE, multi_query = FALSE,
               sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"))
gem = gostres$result[,c("term_id", "term_name", "p_value", "intersection")]
colnames(gem) = c("GO.ID", "Description", "p.Val", "Genes")
gem$FDR = gem$p.Val
gem$Phenotype = "+1"
gem = gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
# saving the GEM file
write.table(gem, file = "/Users/li11/myGit/scRNAseqYouTuber/sampleAnalysisResults/gProfiler_gem.txt", sep = "\t", quote = F, row.names = F)

download.file(
  url = "http://biit.cs.ut.ee/gprofiler/static/gprofiler_full_hsapiens.ENSG.gmt", 
  "/Users/li11/myGit/scRNAseqYouTuber/sampleAnalysisResults/gProfiler_full_hs_ensg.gmt")

