# install.packages('tidyverse')
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# For visualisation
# install.packages('pheatmap')
# install.packages("DOSE")
# install.packages("enrichplot")
# install.packages("ggupset")


## https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofiler/
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
out_path <- "x:/myGit/scRNAseqYouTuber/PEA/Results/" # output path, where you want your results exported to
bg_path <- "x:/myGit/scRNAseqYouTuber/PEA/Background_genes/" # folder with your background genes used for PEA

# Functions ------------------------------------------------
## Function: Adjacency matrix to list ####
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Read in data ===================================================
list.files(in_path)
df <- read.csv(paste0(in_path, 'severevshealthy_degresults-1.txt'), row.names = 1)
# Annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2fc > 0 & padj < 0.05 ~ 'UP',
  log2fc < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))

head(df)
genes_in_data <- df$gene_symbol
# Get the genes that are present in your dataframe


genes_in_data <- df$gene_symbol

# Read in the .gmt file

file <- "x:/project2024/PathwayAnalysisShortCourse/GSEAdb/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"

##  from statsquid
## file <- "PEA/Background_genes/c2.cp.kegg.v7.5.1.symbols.gmt"

##  GSEA mouse msigDB for David's project
 file <- "x:/project2024/DavidRNAseq/GSEA_msigDB_mm/m5.go.bp.v2023.2.Mm.symbols.gmt"

 
pwl2 <- read.gmt(file) 
dim(pwl2)
# [1] 12797     2
# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_data,] 
# Save the filtered background gene set
dim(pwl2)
# [1] 12551     2
filename <- 'x:/project2024/PathwayAnalysisShortCourse/GSEAdb/kegg.RDS'
saveRDS(pwl2, filename)


# SQUIDTIP! If you want to parse several .gmt files at once, you can use a loop:
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
for (file in gmt_files){
  file <- gmt_files[1]
  pwl2 <- read.gmt(file) 
  pwl2 <- pwl2[pwl2$gene %in% genes_in_data,]
  filename <- paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
  saveRDS(pwl2, filename)
}



## Prepare deg results -----------------------------------------------

# If you forgot to do that before, annotate according to differential expression
df <- df %>% mutate(diffexpressed = case_when(
  log2fc > 0 & padj < 0.05 ~ 'UP',
  log2fc < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
  
))

# Remove non-significant genes
df <- df[df$diffexpressed != 'NO', ]

# Substitute names so they are annotated nicely in the heatmap later
df$diffexpressed <- gsub('DOWN', 'Healthy', gsub('UP', 'Severe', df$diffexpressed))
unique(df$diffexpressed)

# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
deg_results_list <- split(df, df$diffexpressed)

## Run ClusterProfiler -----------------------------------------------

# Settings
name_of_comparison <- 'severevshealthy' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes <- readRDS(paste0(bg_path, 'reactome.RDS')) # read in the background genes

background_genes <- 'KEGG' # for our filename
bg_genes <- readRDS('x:/project2024/PathwayAnalysisShortCourse/GSEAdb/kegg.RDS')

padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways
out_path <- 'x:/project2024/PathwayAnalysisShortCourse/'
filename <- paste0(out_path, 'clusterProfiler/', name_of_comparison, '_', background_genes) # filename of our PEA results

# Run clusterProfiler on each sub-dataframe
res <- lapply(names(deg_results_list),
              function(x) enricher(gene = deg_results_list[[x]]$gene_symbol,
                                   TERM2GENE = bg_genes))
names(res) <- names(deg_results_list)

#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)

res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(res_df)))


# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df <- res_df[res_df$ID %in% target_pws, ]

dim(res_df)
print('Saving clusterprofiler results')

write.csv(res_df, paste0(filename, '_resclusterp.csv'), row.names = FALSE)

