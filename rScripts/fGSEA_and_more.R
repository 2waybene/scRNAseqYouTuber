

##====================================
##  fgsea -- video 15
##  Pathway analysis
##====================================

##  from bioconductor
# https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

library(fgsea)
library(data.table)
library(ggplot2)

data("examplePathways")
data("exampleRanks")
set.seed(1234)


library(Seurat)
library(tidyverse)
library(fgsea)



## the following was record how the matrix data was prepared from scRNAseq project
##  see /Users/li11/myGit/scRNAseqYouTuber/rScripts/scRNAseq_tutorial_script_part_II.R

# NML <- readRDS( "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_integrated.RDS")
# DimPlot (NML, reduction = 'umap', label = TRUE)
# View(NML@meta.data)


# readRDS (Myofibrolast in fibrotic lung)
# NML <- readRDS("../Desktop/GSE132771/IPF_Fibro.combined.RDS")
# DimPlot(NML, reduction = "umap", label = TRUE)

# FeaturePlot(NML, features=c("CTHRC1", "POSTN", "COL1A1", "ASPN"),
#            label = TRUE, label.size = 2, cols = c('lightgrey', 'blue'))

# Identify differentially expressed genes (DEG)
# DefaultAssay(NML) <- "RNA"

# DEG_MyoF <- FindMarkers(NML, ident.1 = "5", 
#                        logfc.threshold = 0.25, min.pct = 0.25) 

# DEG_MyoF$gene <- rownames(DEG_MyoF)
# DEG_MyoF <- DEG_MyoF %>% arrange(desc(avg_log2FC))
# View(DEG_MyoF)
#write.csv(DEG_MyoF, file="/Users/jyli/SingleCellTutorial/workingData/DEG_MyoF.CSV")


##  my work MacPro
DEG_MyoF <- read.csv (file = "/Users/li11/myGit/scRNAseqYouTuber/testData/DEG_MyoF.CSV")
View(DEG_MyoF)


fold_changes <- DEG_MyoF$avg_log2FC
names(fold_changes) <- DEG_MyoF$gene

head(fold_changes)

# Load GSEA gene sets: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
# MonsterMac
# Reactome <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
# hallmark <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/h.all.v2023.2.Hs.symbols.gmt")
# KEGG <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
#GO <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/c5.go.cc.v2023.2.Hs.symbols.gmt")

##  my work MacPro
Reactome <- fgsea::gmtPathways("/Users/li11/projectsOnMac/project2024/pathwayAnalysis/GSEAdb/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
hallmark <- fgsea::gmtPathways("/Users/li11/projectsOnMac/project2024/pathwayAnalysis/GSEAdb/h.all.v2023.2.Hs.symbols.gmt")
KEGG <- fgsea::gmtPathways("/Users/li11/projectsOnMac/project2024/pathwayAnalysis/GSEAdb/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
GO <- fgsea::gmtPathways("/Users/li11/projectsOnMac/project2024/pathwayAnalysis/GSEAdb/c5.go.cc.v2023.2.Hs.symbols.gmt")

Four_gene_sets <- c(hallmark, KEGG, GO, Reactome)

## https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/
Allgenesets <- fgsea::gmtPathways("/Users/li11/projectsOnMac/project2024/pathwayAnalysis/GSEAdb/msigdb.v2023.2.Hs.symbols.gmt")

# GSEA analysis
gsea_Myo <- fgsea(pathways = Reactome,
                  stats = fold_changes,
                  eps = 0.0,
                  minSize=15,
                  maxSize=500)
View(gsea_Myo)

head(gsea_Myo[order(pval), ])

# make a table plot for a bunch of selected pathways:
topPathwaysUp <- gsea_Myo[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- gsea_Myo[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(Reactome [topPathways], fold_changes, gsea_Myo, 
              gseaParam=0.5)

# make an enrichment plot for a pathway:
plotEnrichment(Reactome[["REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"]],
               fold_changes) + labs(title="ECM_ORGANIZATION")

# Save gsea_Myo
gsea_Myo <- apply(gsea_Myo, 2, as.character)
write.table(gsea_Myo, file="/Users/li11/myGit/scRNAseqYouTuber/sampleAnalysisResults/gsea_MyoFReactome.CSV",
            append = FALSE, quote = TRUE, sep = ",", row.names = F, col.names = TRUE)


