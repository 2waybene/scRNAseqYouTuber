##=========================================================
##  Tutorial video 3
##  Convert Seurat object to cds for trajectory analysis
##=========================================================

library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(monocle3)

# 1. Analyze scRNA-seq data using Seurat------
# Read in 10x raw file
MG2 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG2-120-2/")
MG4 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG4-121-7/")
MG6 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG6-123-4/")
MG8 <- Read10X(data.dir = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/CellRanger/MG8-124-4/")


# 2. Create Seurat object -----
# 
MG2 <- CreateSeuratObject(counts = MG2, project = "MG2", 
                          min.cells = 3, min.features = 200) 
MG2 <- PercentageFeatureSet(MG2, pattern = "^mt-", col.name = "percent.mt")

MG4 <- CreateSeuratObject(counts = MG4, project = "MG4", 
                          min.cells = 3, min.features = 200) 
MG4 <- PercentageFeatureSet(MG4, pattern = "^mt-", col.name = "percent.mt")

MG6 <- CreateSeuratObject(counts = MG6, project = "MG6", 
                          min.cells = 3, min.features = 200) 
MG6 <- PercentageFeatureSet(MG6, pattern = "^mt-", col.name = "percent.mt")

MG8 <- CreateSeuratObject(counts = MG8, project = "MG8", 
                          min.cells = 3, min.features = 200) 
MG8 <- PercentageFeatureSet(MG8, pattern = "^mt-", col.name = "percent.mt")

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

fFlox <- merge(MG2, y = c(MG4, MG6, MG8))

VlnPlot(fFlox, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fFlox <- subset(fFlox, nFeature_RNA<7500 & nCount_RNA<50000 & percent.mt< 7.5)
dim(fFlox@meta.data)
table(fFlox@meta.data$orig.ident)
View(fFlox@meta.data)

fFlox.list <- SplitObject(fFlox, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
fFlox.list <- lapply(X = fFlox.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = fFlox.list)
fFlox.anchors <- FindIntegrationAnchors(object.list = fFlox.list, 
                                      anchor.features = features)

# creates an 'integrated' data assay
fFlox.combined <- IntegrateData(anchorset = fFlox.anchors)

##  B. Normalize the integrate data, dimension reduction, and clustering 
# specify that we will perform downstream analysis on the integrated data
DefaultAssay(fFlox.combined) <- "integrated"

# Run the standard workflow for cell clustering
fFlox.combined <- ScaleData(fFlox.combined, verbose = FALSE)
fFlox.combined <- RunPCA(fFlox.combined, npcs = 30, verbose = FALSE)
fFlox.combined <- RunUMAP(fFlox.combined, reduction = "pca", dims = 1:30)
fFlox.combined <- FindNeighbors(fFlox.combined, reduction = "pca", dims = 1:30)
fFlox.combined <- FindClusters(fFlox.combined, resolution = 0.5)

# Visualization
DimPlot(fFlox.combined, reduction = "umap", label = TRUE)
DimPlot(fFlox.combined, reduction = "umap", label = TRUE, split.by = "orig.ident")
##  C. Feature plot to identify megakaryocets celll populations 
saveRDS (fFlox.combined, file = "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/processedData/femaleFlox4samples_standardProcessed.rds")

##  C. Feature plot to identify megakaryocets celll populations

fFlox.combined <- readRDS ( "/Volumes/li11/project2023/CidlowskiRNAseq/MatiasProject/combinedData/processedData/femaleFlox4samples_standardProcessed.rds")
DimPlot(fFlox.combined, reduction = "umap", label = TRUE)

fFlox.combined <- FindClusters(fFlox.combined, resolution = 0.1)
DimPlot(fFlox.combined, reduction = "umap", label = TRUE)


