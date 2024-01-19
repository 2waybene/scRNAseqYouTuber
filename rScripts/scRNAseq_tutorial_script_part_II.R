##========================================
## Seurat video tutorial -- video 10
## Find markers
##========================================

library(Seurat)
library(tidyverse)
library(patchwork)

IntegratdNML <- readRDS( "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_integrated.RDS")
IntegratdNML <- FindClusters(IntegratdNML, resolution = 0.01)
DimPlot (IntegratdNML, reduction = 'umap', label = TRUE)


FeaturePlot(IntegratdNML, features = c ("EPCAM", "CLDN5", "COL1A2", "PTPRC"),
            cols = c('lightgrey', 'blue'))

DefaultAssay(IntegratdNML) <-"integrated"
DefaultAssay(IntegratdNML) <- 'RNA'

cluster2.markers <- FindMarkers(IntegratdNML, ident.1 = 2, min.pct = 0.5,
                                only.pos = TRUE, logfc.threshold = 0.25)
ALL.makers <- FindAllMarkers(IntegratdNML, min.pct = 0.5, only.pos = TRUE,
                             logfc.threshold = 0.25)

View(cluster2.markers)
FeaturePlot(IntegratdNML, features = c ("SFTPC", "SFTPB", "NAPSA", "SCGB3A2"),
            cols = c('lightgrey', 'blue'))
DotPlot(IntegratdNML, features = c("SFTPC", "SFTPB", "NAPSA", "SCGB3A2", "AGER"))
VlnPlot(IntegratdNML, features = c("SFTPC", "SFTPB", "NAPSA", "SCGB3A2", "AGER"))

cluster0_3.markers <- FindMarkers(IntegratdNML, ident.1 = 0,  ident.2 = 3, min.pct = 0.25)

cluster0_3and6.markers <- FindMarkers(IntegratdNML, ident.1 = 0,  ident.2 = c(3,6) , min.pct = 0.25)

View(IntegratdNML@commands)

##===================================================================
## Rename cell cluster and use your favorite colors -- video 11
##====================================================================


library(Seurat)
library(tidyverse)
library(patchwork)

IntegratdNML <- readRDS( "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_integrated.RDS")
IntegratdNML <- FindClusters(IntegratdNML, resolution = 0.01)
DimPlot (IntegratdNML, reduction = 'umap', label = TRUE)
FeaturePlot(IntegratdNML, features = c ("EPCAM", "CLDN5", "COL1A2", "PTPRC"),
            cols = c('lightgrey', 'blue'))

cluster5.markers <- FindMarkers(IntegratdNML, ident.1 = 5, min.pct = 0.5,
                                only.pos = TRUE, logfc.threshold = 0.25)
View(cluster5.markers)
FeaturePlot(IntegratdNML, features = c ("MS4A2", "CPA3", "TPSAB1", "TPSAB2"),
            cols = c('grey', 'red'))


cluster6.markers <- FindMarkers(IntegratdNML, ident.1 = 6, min.pct = 0.5,
                                only.pos = TRUE, logfc.threshold = 0.25)
View(cluster6.markers)
FeaturePlot(IntegratdNML, features = c ("CCL21", "NNMT", "TFF3", "MMRN1"),
            cols = c('grey', 'red'))


IntegratdNML <- RenameIdents(IntegratdNML, '0'="Immune", '3'="Immune", '5' = "Immune",
                             '1'="Endothelial", '6'="Lymphatic_Endo", '2' = "Epithelial",
                             '4'="Mensenchymal")
DimPlot (IntegratdNML, reduction = 'umap', label = TRUE)

## chose your favorite color for the cell clusters
## https://htmlcolorcodes.com/color-names/

p1 <- DimPlot (IntegratdNML, reduction = 'umap', label = TRUE, label.size = 3, cols = 
                 c("Red", "Green", "DeepPink", "Blue", "Orange", "DeepSkyBlue"))
p2 <- DimPlot (IntegratdNML, reduction = 'umap', label = TRUE, label.size = 3)
p1 + p2


##====================================
##  About subsetting -- video 12
##====================================
library(Seurat)
library(tidyverse)
library(patchwork)

MergeDNML<- readRDS("/Users/jyli/SingleCellTutorial/workingData/MergeDNML_QC.RDS")
MergeDNML <- subset (MergeDNML, subset =  nFeature_RNA > 500 & nFeature_RNA < 4000
                     & nCount_RNA < 20000 & percent.mt < 10 )


IntegratdNML <- readRDS( "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_integrated.RDS")
IntegratdNML <- FindClusters(IntegratdNML, resolution = 0.01)
DimPlot (IntegratdNML, reduction = 'umap', label = TRUE)

## subset on the expression level of a gene/feature
TypeII_Epi_cells <- subset (x = IntegratdNML, subset = SFTPB > 4)
DimPlot (TypeII_Epi_cells, reduction = 'umap', label = TRUE)

Immune_cells <- subset (x = IntegratdNML, subset = CD52 > 2.8)
DimPlot (Immune_cells, reduction = 'umap', label = TRUE)

Immune_clusters <- subset (x= IntegratdNML, idents = c('0', '3', '5'))
DimPlot (Immune_clusters, reduction = 'umap', label = TRUE)

Mesenchymal_cells <- subset (x = IntegratdNML, idents = '4')
DimPlot (Mesenchymal_cells, reduction = "umap", label = TRUE)

# Run the standard workflow for visualization and clustering
Mesenchymal_cells <- ScaleData(Mesenchymal_cells, verbose = FALSE)
Mesenchymal_cells <- RunPCA(Mesenchymal_cells, npcs = 50,  verbose = FALSE)
Mesenchymal_cells <- FindNeighbors(Mesenchymal_cells, 
                                      reduction = 'pca', dims = 1:30)
Mesenchymal_cells <- FindClusters(Mesenchymal_cells, resolution = 0.1)
Mesenchymal_cells <- RunUMAP(Mesenchymal_cells, reduction = 'pca', dims = 1:30)
DimPlot (Mesenchymal_cells, reduction = "umap", label = TRUE)

FeaturePlot(Mesenchymal_cells, features  = c('PDGFRA', 'LUM', "DES", "RGS5"),
            label = TRUE, cols = c('lightgrey', 'blue'))
saveRDS(Mesenchymal_cells, file = "/Users/jyli/SingleCellTutorial/workingData/Mesenchymal_cells,.RDS")


##============================================
##  How to load TXT scRNA-seq data into R
##. Video 13
##============================================
library(Seurat)
library(tidyverse)

Sciatic_nerve <- read.delim ("/Users/jyli/SingleCellTutorial/dataDownload/GSE147285_RAW/GSM4423509_Uninj_Sciatic.txt",
                             header = TRUE, sep="\t")
View(Sciatic_nerve)
GENE <- Sciatic_nerve[["GENE"]]
rownames(Sciatic_nerve) = GENE
Sciatic_nerve[["GENE"]] <- NULL
class(Sciatic_nerve)

Sciatic_nerve <- as.matrix (Sciatic_nerve)
class(Sciatic_nerve)

Sciatic_nerve <- CreateAssayObject(counts = Sciatic_nerve,
                                   project = "CON_GSE147285", min.cells = 3, min.features = 200)

##====================================
##  About subsetting -- video 14
##  Calculate cluster frequency
##====================================
library(Seurat)
library(tidyverse)
library(patchwork)

NML <- readRDS( "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_integrated.RDS")
DimPlot (NML, reduction = 'umap', label = TRUE)
View(NML@meta.data)

##  cluseter frenquecy
cluster_frequency <- NML@meta.data %>% 
  group_by(seurat_clusters) %>% 
  summarise(count = n()) %>%
  mutate(relative_frequency = count/sum(count)) %>%
  mutate(data_set = "NML")

write.csv(cluster_frequency, file = "/Users/jyli/SingleCellTutorial/workingData/cluster_frenquency.csv")

ggplot(cluster_frequency, aes(x=data_set, y=relative_frequency,
                              fill=seurat_clusters)) + geom_col()


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

NML <- readRDS( "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_integrated.RDS")
DimPlot (NML, reduction = 'umap', label = TRUE)
View(NML@meta.data)




# readRDS (Myofibrolast in fibrotic lung)
# NML <- readRDS("../Desktop/GSE132771/IPF_Fibro.combined.RDS")
DimPlot(NML, reduction = "umap", label = TRUE)

FeaturePlot(NML, features=c("CTHRC1", "POSTN", "COL1A1", "ASPN"),
            label = TRUE, label.size = 2, cols = c('lightgrey', 'blue'))

# Identify differentially expressed genes (DEG)
DefaultAssay(NML) <- "RNA"

DEG_MyoF <- FindMarkers(NML, ident.1 = "5", 
                        logfc.threshold = 0.25, min.pct = 0.25) 

DEG_MyoF$gene <- rownames(DEG_MyoF)
DEG_MyoF <- DEG_MyoF %>% arrange(desc(avg_log2FC))
View(DEG_MyoF)
write.csv(DEG_MyoF, file="/Users/jyli/SingleCellTutorial/workingData/DEG_MyoF.CSV")

fold_changes <- DEG_MyoF$avg_log2FC
names(fold_changes) <- DEG_MyoF$gene

head(fold_changes)

# Load GSEA gene sets: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
Reactome <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
hallmark <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/h.all.v2023.2.Hs.symbols.gmt")
KEGG <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
GO <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/c5.go.cc.v2023.2.Hs.symbols.gmt")

Four_gene_sets <- c(hallmark, KEGG, GO, Reactome)
Allgenesets <- fgsea::gmtPathways("/Users/jyli/SingleCellTutorial/dataDownload/GSEAdb/msigdb.v2023.2.Hs.symbols.gmt")

# GSEA analysis
gsea_Myo <- fgsea(pathways = Reactome,
                  stats = fold_changes,
                  eps = 0.0,
                  minSize=15,
                  maxSize=500)

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
write.table(gsea_Myo, file="/Users/jyli/SingleCellTutorial/sampleAnalysisResults/gsea_MyoFReactome.CSV",
            append = FALSE, quote = TRUE, sep = ",", row.names = F, col.names = TRUE)


