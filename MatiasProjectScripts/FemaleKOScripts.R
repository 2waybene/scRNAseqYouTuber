
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(monocle3)

# 1. Analyze scRNA-seq data using Seurat------
# Read in 10x raw file
MG1 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG1-120-1/")
MG3 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG3-121-6/")
MG5 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG5-123-3/")
MG7 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG7-124-3/")


# 2. Create Seurat object -----
# 
MG1 <- CreateSeuratObject(counts = MG1, project = "MG1", 
                          min.cells = 3, min.features = 200) 
MG1 <- PercentageFeatureSet(MG1, pattern = "^mt-", col.name = "percent.mt")

MG3 <- CreateSeuratObject(counts = MG3, project = "MG3", 
                          min.cells = 3, min.features = 200) 
MG3 <- PercentageFeatureSet(MG3, pattern = "^mt-", col.name = "percent.mt")

MG5 <- CreateSeuratObject(counts = MG5, project = "MG5", 
                          min.cells = 3, min.features = 200) 
MG5 <- PercentageFeatureSet(MG5, pattern = "^mt-", col.name = "percent.mt")

MG7 <- CreateSeuratObject(counts = MG7, project = "MG7", 
                          min.cells = 3, min.features = 200) 
MG7 <- PercentageFeatureSet(MG7, pattern = "^mt-", col.name = "percent.mt")

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

fKO <- merge(MG1, y = c(MG3, MG5, MG7))



VlnPlot(fKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE)
fKO <- subset(fKO, nFeature_RNA<7500 & nCount_RNA<50000 & percent.mt< 7.5)
dim(fKO@meta.data)
table(fKO@meta.data$orig.ident)
View(fKO@meta.data)



rm (MG1, MG3, MG5, MG7)

fKO.list <- SplitObject(fKO, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
fKO.list <- lapply(X = fKO.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = fKO.list)
fKO.anchors <- FindIntegrationAnchors(object.list = fKO.list, 
                                        anchor.features = features)

# creates an 'integrated' data assay
fKO.combined <- IntegrateData(anchorset = fKO.anchors)

##  B. Normalize the integrate data, dimension reduction, and clustering 
# specify that we will perform downstream analysis on the integrated data
DefaultAssay(fKO.combined) <- "integrated"

# Run the standard workflow for cell clustering
fKO.combined <- ScaleData(fKO.combined, verbose = FALSE)
fKO.combined <- RunPCA(fKO.combined, npcs = 30, verbose = FALSE)
fKO.combined <- RunUMAP(fKO.combined, reduction = "pca", dims = 1:30)
fKO.combined <- FindNeighbors(fKO.combined, reduction = "pca", dims = 1:30)
fKO.combined <- FindClusters(fKO.combined, resolution = 0.5)

# Visualization
DimPlot(fKO.combined, reduction = "umap", label = TRUE)
DimPlot(fKO.combined, reduction = "umap", label = TRUE, split.by = "orig.ident")
##  C. Feature plot to identify megakaryocets celll populations 
saveRDS (fKO.combined, file = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/processedData/femaleKO4samples_standardProcessed.rds")
