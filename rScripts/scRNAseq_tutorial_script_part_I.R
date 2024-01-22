##========================================================================================
##  scRNAseq_tutorial_script_part_I.R
##  YouTube video tutorial 1- 7
##  https://www.youtube.com/watch?v=43Z13DS_emQ&list=PLOLdjuxsfI4N1SdaQQYXGoa5Z93hPxWVY
##========================================================================================


library(Seurat)

## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771

##====================================
##  Read in raw data -- video 1
##====================================


##  MacMonster
NML_1 <- Read10X (data.dir = "/Users/jyli/SingleCellTutorial/dataDownload/GSE132771_RAW/NML1")
NML_2 <- Read10X (data.dir = "/Users/jyli/SingleCellTutorial/dataDownload/GSE132771_RAW/NML2")
NML_3 <- Read10X (data.dir = "/Users/jyli/SingleCellTutorial/dataDownload/GSE132771_RAW/NML3")

ls()

##  create seurate objectives

##=========================================
##  create seurate objectives-- video 2
##=========================================

NML_I <- CreateSeuratObject(counts = NML_1, project = "NML_I", min.cells = 3, min.features = 200)
class(NML_I)
colnames(NML_I)
rownames(NML_I)
View(NML_I)
View(NML_I@meta.data)

NML_II <- CreateSeuratObject(counts = NML_2, project = "NML_II", min.cells = 3, min.features = 200)
NML_III <- CreateSeuratObject(counts = NML_3, project = "NML_III", min.cells = 3, min.features = 200)

##  save seurate objectives

saveRDS(NML_I, file = "/Users/jyli/SingleCellTutorial/workingData/NML_I.RDS")
saveRDS(NML_II, file = "/Users/jyli/SingleCellTutorial/workingData/NML_II.RDS")
saveRDS(NML_III, file = "/Users/jyli/SingleCellTutorial/workingData/NML_III.RDS")


##=========================================
##  merge three seurate objectives video 3
##==========================================

library(Seurat)
library(tidyverse)
NML_I <- readRDS("/Users/jyli/SingleCellTutorial/workingData/NML_I.RDS")
NML_II <- readRDS("/Users/jyli/SingleCellTutorial/workingData/NML_II.RDS")
NML_III <- readRDS("/Users/jyli/SingleCellTutorial/workingData/NML_III.RDS")

MergeDNML <- merge (NML_I, y = c(NML_II, NML_III),
                   add.cell.ids = ls()[1:3],
                   project = "MergeDNML")

ls()
MergeDNML
View(MergeDNML@meta.data)
saveRDS(MergeDNML, file = "/Users/jyli/SingleCellTutorial/workingData/MergeDNML.RDS")


##====================================
##  Data QC -- video 4
##====================================
library(Seurat)
library(tidyverse)
MergeDNML<- readRDS("/Users/jyli/SingleCellTutorial/workingData/MergeDNML.RDS")

View(MergeDNML@meta.data)
range(MergeDNML$nFeature_RNA)
range(MergeDNML$nCount_RNA)


MergeDNML <- PercentageFeatureSet(MergeDNML, pattern = "^MT-", col.name = "percent.mt")
View(MergeDNML@meta.data)
range(MergeDNML$percent.mt)

VlnPlot(MergeDNML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3)
MergeDNML <- subset (MergeDNML, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
MergeDNML
saveRDS(MergeDNML, file = "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_QC.RDS")


##====================================
##  normalize the data -- video 5
##====================================
library(Seurat)
library(tidyverse)
MergeDNML<- readRDS("/Users/jyli/SingleCellTutorial/workingData/MergeDNML_QC.RDS")

View(MergeDNML)
MergeDNML@assays[["RNA"]]@data@x

MergeDNML <- NormalizeData(MergeDNML,
                           normalization.method = "LogNormalize", scale.factor = 10000)

MergeDNML@assays[["RNA"]]@data@x
MergeDNML

MergeDNML <- FindVariableFeatures(MergeDNML,
                           selection.method = "vst", scale.factor = 10000)
MergeDNML
MergeDNML@assays[["RNA"]]@var.features
VariableFeaturePlot(MergeDNML)

MergeDNML <- ScaleData(MergeDNML) #2000 identified variable features

## Not run
all.genes <- rownames(MergeDNML)
MergeDNML <- ScaleData(MergeDNML, features = all.genes)


MergeDNML@assays[["RNA"]]@scale.data
View(MergeDNML@commands)

saveRDS(MergeDNML, file = "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_Normed.RDS")


##=====================================
##  Dimension reduction -- video 6
##=====================================
library(Seurat)
library(tidyverse)
MergeDNML<- readRDS("/Users/jyli/SingleCellTutorial/workingData/MergeDNML_Normed.RDS")

MergeDNML <- RunPCA(MergeDNML)
DimPlot(MergeDNML, reduction = "pca", dim = c(1,2))
DimPlot(MergeDNML, reduction = "pca", dim = c(1,10))
DimPlot(MergeDNML, reduction = "pca", dim = c(1,50))

MergeDNML <- JackStraw(MergeDNML, num.replicate = 100)
MergeDNML <- ScoreJackStraw(MergeDNML, dims = 1:20)
## somehow failed here!!
JackStrawPlot(MergeDNML, dimes =1:20)

ElbowPlot(MergeDNML)
ElbowPlot(MergeDNML, ndims = 50, reduction = "pca")

MergeDNML <- FindNeighbors(MergeDNML, dims = 1:20)
MergeDNML <- FindClusters(MergeDNML, resolution = 0.1)
MergeDNML <- FindClusters(MergeDNML, resolution = 0.3)

##  seurate suggests that using the same dims in UMAP as in FindNeighbors
MergeDNML <- RunUMAP(MergeDNML, dims = 1:20)
DimPlot(MergeDNML, reduction = 'umap', label = TRUE , repel = TRUE)


MergeDNML <- RunTSNE(object = MergeDNML)
DimPlot(object = MergeDNML, reduction = 'tsne')
View(MergeDNML)

saveRDS(MergeDNML, file = "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_dimDuc.RDS")



##=====================================
##  Data integration -- video 7
##=====================================
library(Seurat)
library(tidyverse)
MergeDNML_SWF <- readRDS("/Users/jyli/SingleCellTutorial/workingData/MergeDNML_dimDuc.RDS")
DimPlot(MergeDNML_SWF, reduction = 'umap', label = TRUE)

FeaturePlot(MergeDNML_SWF, features = c ("EPCAM", "CLDN5", "COL1A2", "PTPRC"),
            cols = c('lightgrey', 'blue'))

view(MergeDNML_SWF@meta.data)
DimPlot (MergeDNML_SWF, reduction = 'umap', label = TRUE, group.by = "orig.ident")

MergeDNML<- readRDS("/Users/jyli/SingleCellTutorial/workingData/MergeDNML_QC.RDS")
View(MergeDNML@meta.data)
View(MergeDNML)

MergeDNML.list <- SplitObject(MergeDNML, split.by = 'orig.ident')

MergeDNML.list <- lapply (X = MergeDNML.list, FUN = function(x) {
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method = "vst", scale.factor = 10000)
})

features <- SelectIntegrationFeatures(object.list = MergeDNML.list)

MergeDNML.anchors <- FindIntegrationAnchors(object.list = MergeDNML.list,
                                            anchor.features = features)

MergedNML.integrated <- IntegrateData (anchorset = MergeDNML.anchors)

### Perform an integrated analysis
# specify that we will perform downstream analysis on the integrated data
DefaultAssay(MergedNML.integrated) <- 'integrated'

# Run the standard workflow for visualization and clustering
MergedNML.integrated <- ScaleData(MergedNML.integrated, verbose = FALSE)
MergedNML.integrated <- RunPCA(MergedNML.integrated, npcs = 50,  verbose = FALSE)
MergedNML.integrated <- FindNeighbors(MergedNML.integrated, 
                                      reduction = 'pca', dims = 1:30)
MergedNML.integrated <- FindClusters(MergedNML.integrated, resolution = 0.3)
MergedNML.integrated <- RunUMAP(MergedNML.integrated, reduction = 'pca', dims = 1:30)

DimPlot(MergedNML.integrated, reduction = 'umap', lqbel = TRUE)

plot1 <- DimPlot(MergeDNML_SWF, reduction = 'umap', group.by = 'orig.ident')
plot2 <- DimPlot(MergedNML.integrated, reduction = 'umap', group.by = 'orig.ident')
plot1 + plot2

saveRDS(MergedNML.integrated, file = "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_integrated.RDS")




