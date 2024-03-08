

####================================================================
## https://biostatsquid.com/pathway-enrichment-analysis-plots/
# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

# Set input path
in_path <- "x:/myGit/scRNAseqYouTuber/testData/" # input path, where your data is located
#out_path <- "x:/myGit/scRNAseqYouTuber/PEA/Results/" # output path, where you want your results exported to
bg_path <- "x:/project2024/PathwayAnalysisShortCourse/GSEAdb/" # folder with your background genes used for PEA
out_path <- 'x:/project2024/PathwayAnalysisShortCourse/' # output path, where you want your results exported to


name_of_comparison <- 'severevshealthy' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS')) # read in the background genes
filename <- paste0(out_path, 'clusterProfiler/', name_of_comparison, '_', background_genes) # filename of our PEA results

# Read in the data
res_df <- read.csv(paste0(out_path, 'clusterProfiler/', 'severevshealthy_reactome_resclusterp.csv'))
bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS'))
# Convert clusterProfiler object to a new "enrichResult" object
# Select only upregulated genes in Severe

res_df %>% filter(diffexpressed == 'Severe') %>% dim()

res_df <- res_df %>% filter(diffexpressed == 'Severe') %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 

head(res_df)
rownames(res_df) <- res_df$ID

# For visualisation purposes, let's shorten the pathway names
res_df$Description <- gsub('(H|h)iv', 'HIV', 
                           gsub('pd 1', 'PD-1',
                                gsub('ecm', 'ECM', 
                                     gsub('(I|i)nterleukin', 'IL', 
                                          gsub('(R|r)na', 'RNA', 
                                               gsub('(D|d)na', 'DNA',
                                                    gsub(' i ', ' I ', 
                                                         gsub('(A|a)tp ', 'ATP ', 
                                                              gsub('(N|n)adh ', 'NADH ', 
                                                                   gsub('(N|n)ad ', 'NAD ',
                                                                        gsub('t cell', 'T cell',
                                                                             gsub('b cell', 'B cell',
                                                                                  gsub('built from .*', ' (...)',
                                                                                       gsub('mhc', 'MHC',
                                                                                            gsub('mhc class i', 'MHC I', 
                                                                                                 gsub('mhc class ii', 'MHC II', 
                                                                                                      stringr::str_to_sentence(
                                                                                                        gsub('_', ' ',  
                                                                                                             gsub('GOBP_|KEGG_|REACTOME_', '', res_df$Description)))))))))))))))))))



enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = res_df,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = res_df$geneID,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)
class(enrichres)


# Barplot
barplot(enrichres, showCategory = 20) 
mutate(enrichres, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")


# Dotplot
dotplot(enrichres, showCategory = 15) + ggtitle("Severe vs Healthy")

# Cnetplot
cnetplot(enrichres)
# Heatplot
heatplot(enrichres, showCategory = 5)
# Treeplot
enrichres2 <- pairwise_termsim(enrichres) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
treeplot(enrichres2)
# Enrichment map 
emapplot(enrichres2)
# Upsetplot
upsetplot(enrichres)
