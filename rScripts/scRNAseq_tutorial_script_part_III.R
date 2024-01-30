# load libraries
library(Seurat)
library(tidyverse)

# Read
CON_I<- Read10X(data.dir = "E:/SingleCellPublicData/dataDownload/GSE132771_RAW/NMLAll1/")
CON_II<- Read10X(data.dir = "E:/SingleCellPublicData/dataDownload/GSE132771_RAW/NMLAll2/")
CON_III<- Read10X(data.dir = "E:/SingleCellPublicData/dataDownload/GSE132771_RAW/NMLAll3/")
IPF_I<- Read10X(data.dir = "E:/SingleCellPublicData/dataDownload/GSE132771_RAW/IPFAll1/")
IPF_II<- Read10X(data.dir = "E:/SingleCellPublicData/dataDownload/GSE132771_RAW/IPFAll2/")
IPF_III<- Read10X(data.dir = "E:/SingleCellPublicData/dataDownload/GSE132771_RAW/IPFAll3/")

# Create Seurat Object
CON_I <- CreateSeuratObject(counts = CON_I, project = "CON", min.cells = 3, min.features = 200)
CON_I <- PercentageFeatureSet(CON_I, pattern = "^MT-", col.name = "percent.mt")

CON_II <- CreateSeuratObject(counts = CON_II, project = "CON", min.cells = 3, min.features = 200)
CON_II <- PercentageFeatureSet(CON_II, pattern = "^MT-", col.name = "percent.mt")

CON_III <- CreateSeuratObject(counts = CON_III, project = "CON", min.cells = 3, min.features = 200)
CON_III <- PercentageFeatureSet(CON_III, pattern = "^MT-", col.name = "percent.mt")

IPF_I <- CreateSeuratObject(counts = IPF_I, project = "IPF", min.cells = 3, min.features = 200)
IPF_I <- PercentageFeatureSet(IPF_I, pattern = "^MT-", col.name = "percent.mt")

IPF_II <- CreateSeuratObject(counts = IPF_II, project = "IPF", min.cells = 3, min.features = 200)
IPF_II <- PercentageFeatureSet(IPF_II, pattern = "^MT-", col.name = "percent.mt")

IPF_III <- CreateSeuratObject(counts = IPF_III, project = "IPF", min.cells = 3, min.features = 200)
IPF_III <- PercentageFeatureSet(IPF_III, pattern = "^MT-", col.name = "percent.mt")

# Filter low quality cells
VlnPlot(CON_I, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CON_I <- subset(CON_I, subset = nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)

VlnPlot(CON_II, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CON_II <- subset(CON_II, subset = nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 10)

VlnPlot(CON_III, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CON_III <- subset(CON_III, subset = nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)

VlnPlot(IPF_I, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
IPF_I <- subset(IPF_I, subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & percent.mt < 10)

VlnPlot(IPF_II, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
IPF_II <- subset(IPF_II, subset = nFeature_RNA < 3000 & nCount_RNA < 15000 & percent.mt < 5)

VlnPlot(IPF_III, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
IPF_III <- subset(IPF_III, subset = nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)

view(CON_I@meta.data)
view(IPF_I@meta.data)

CONIPF.list <- list(CON_I = CON_I, CON_II = CON_II, CON_III = CON_III,
                    IPF_I = IPF_I, IPF_II = IPF_II, IPF_III = IPF_III)

rm(CON_I, CON_II, CON_III, IPF_I, IPF_II, IPF_III)

# normalize and identify variable features for each dataset independently
CONIPF.list <- lapply(X = CONIPF.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = CONIPF.list)
CONIPF.anchors <- FindIntegrationAnchors(object.list = CONIPF.list, anchor.features = features)

# creates an 'integrated' data assay
CONIPF.combined <- IntegrateData(anchorset = CONIPF.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(CONIPF.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
CONIPF.combined <- ScaleData(CONIPF.combined, verbose = FALSE)
CONIPF.combined <- RunPCA(CONIPF.combined, npcs = 50, verbose = FALSE)
CONIPF.combined <- RunUMAP(CONIPF.combined, reduction = "pca", dims = 1:30)
CONIPF.combined <- FindNeighbors(CONIPF.combined, reduction = "pca", dims = 1:30)
CONIPF.combined <- FindClusters(CONIPF.combined, resolution = 0.1)
DimPlot(CONIPF.combined, reduction = "umap", label = TRUE)

VlnPlot(CONIPF.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

view(CONIPF.combined@meta.data)

DimPlot(CONIPF.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(CONIPF.combined, reduction = "umap", split.by = "orig.ident")

FeaturePlot(CONIPF.combined, features=c("EPCAM", "SFTPC", "AGER", "SCGB1A1"),
            cols = c('lightgrey', 'blue'))
FeaturePlot(CONIPF.combined, features=c("CLDN5", "CCL21", "PECAM1", "EMCN"),label = TRUE,
            cols = c('lightgrey', 'blue'))
FeaturePlot(CONIPF.combined, features=c("COL1A2", "LUM", "PDGFRA","PDGFRB"),label = TRUE,
            cols = c('lightgrey', 'blue'))
FeaturePlot(CONIPF.combined, features=c("PTPRC", "CD52", "AIF1", "TRBC2"),label = TRUE,
            cols = c('lightgrey', 'blue'))
FeaturePlot(CONIPF.combined, features=c("MSLN", "CALB2", "HP", "PRG4"),label = TRUE,
            cols = c('lightgrey', 'blue'))

CONIPF.combined <- RenameIdents(CONIPF.combined, 
                                `4` = "EPI", `9` = "EPI", `1` = "ENDO",`12` = "ENDO", 
                                `5` = "MESN", `6` = "MESN", `14` = "MESN", 
                                `0` = "IM", `2` = "IM", `3` = "IM", `7` = "IM",
                                `8` = "IM",`10` = "IM",`11` = "IM",`13` = "IM")

# Save CONIPF.combined.RDS
saveRDS(CONIPF.combined, file="../Desktop/Seurat Video Tutorials/Data/CONIPF.combined.RDS")
CONIPF.combined <- readRDS("../Desktop/Seurat Video Tutorials/Data/CONIPF.combined.RDS")