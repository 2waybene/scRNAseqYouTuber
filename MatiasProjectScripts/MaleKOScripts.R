
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(monocle3)

# 1. Analyze scRNA-seq data using Seurat------
# Read in 10x raw file
MG9 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG9-121-1/")
MG10 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG10-122-1/")
MG12 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG12-124-1/")
MG14 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG14-125-2/")


# 2. Create Seurat object -----
# 
MG9 <- CreateSeuratObject(counts = MG9, project = "MG9", 
                          min.cells = 3, min.features = 200) 
MG9 <- PercentageFeatureSet(MG9, pattern = "^mt-", col.name = "percent.mt")

MG10 <- CreateSeuratObject(counts = MG10, project = "MG10", 
                          min.cells = 3, min.features = 200) 
MG10 <- PercentageFeatureSet(MG10, pattern = "^mt-", col.name = "percent.mt")

MG12 <- CreateSeuratObject(counts = MG12, project = "MG12", 
                          min.cells = 3, min.features = 200) 
MG12 <- PercentageFeatureSet(MG12, pattern = "^mt-", col.name = "percent.mt")

MG14 <- CreateSeuratObject(counts = MG14, project = "MG14", 
                          min.cells = 3, min.features = 200) 
MG14 <- PercentageFeatureSet(MG14, pattern = "^mt-", col.name = "percent.mt")

##==================================================================================
# 3. Merge into a single Seurat object
# and start standard analysis with a few key steps
#
##  A. split into separate object, perform normalization, and integrate them with cca 
##  B. Normalize the integrate data, dimension reduction, and clustering 
##  C. Feature plot to identify megakaryocets celll populations
##  split into separate oject, perform normalization, and integrate them with cca 
##===================================================================================


##  A. split into separate object, perform normalization, and integrate them with cca 

mKO <- merge(MG9, y = c(MG10, MG12, MG14))



VlnPlot(mKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE)
mKO <- subset(mKO, nFeature_RNA<7500 & nCount_RNA<50000 & percent.mt< 7.5)
dim(mKO@meta.data)
table(mKO@meta.data$orig.ident)
View(mKO@meta.data)



rm (MG9, MG10, MG12, MG14)

mKO.list <- SplitObject(mKO, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
mKO.list <- lapply(X = mKO.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = mKO.list)
mKO.anchors <- FindIntegrationAnchors(object.list = mKO.list, 
                                        anchor.features = features)

# creates an 'integrated' data assay
mKO.combined <- IntegrateData(anchorset = mKO.anchors)

##  B. Normalize the integrate data, dimension reduction, and clustering 
# specify that we will perform downstream analysis on the integrated data
DefaultAssay(mKO.combined) <- "integrated"

# Run the standard workflow for cell clustering
mKO.combined <- ScaleData(mKO.combined, verbose = FALSE)
mKO.combined <- RunPCA(mKO.combined, npcs = 30, verbose = FALSE)
mKO.combined <- RunUMAP(mKO.combined, reduction = "pca", dims = 1:30)
mKO.combined <- FindNeighbors(mKO.combined, reduction = "pca", dims = 1:30)
mKO.combined <- FindClusters(mKO.combined, resolution = 0.5)

# Visualization
DimPlot(mKO.combined, reduction = "umap", label = TRUE)
DimPlot(mKO.combined, reduction = "umap", label = TRUE, split.by = "orig.ident")
##  C. Feature plot to identify megakaryocets celll populations 
saveRDS (mKO.combined, file = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/processedData/maleKO4samples_standardProcessed.rds")


Sun.markers <- c("Itga2b", "Pf4", "Itgb3", "Gp9")
MK.markers.minus.sun <- c("Gata1", "Mpl", "Kit", "Cd34", "Vwf", "Tubb1", "Gp1bb")
MK.markers.extra <- c("Ppbp", "Selp", "Gp6", "Gp1ba")


mKO.combined <- FindClusters(mKO.combined, resolution = 0.1)
FeaturePlot(mKO.combined,, features = Sun.markers, label = TRUE)
FeaturePlot(mKO.combined,, features = MK.markers.minus.sun, label = TRUE)
FeaturePlot(mKO.combined,, features = MK.markers.extra , label = TRUE)

