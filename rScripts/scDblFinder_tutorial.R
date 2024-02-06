# library(scDblFinder)
library(Seurat)
suppressPackageStartupMessages(library(scDblFinder))
library(tidyverse)

set.seed(123)


##  Need data path for sample run
# 1. Create Seurat Objects

# 1. Analyze scRNA-seq data using Seurat------
# Read

##  MyWorkMac
NML1 <- Read10X(data.dir = "/Users/li11/myGit/SingleCellDataDownload/GSE132771_RAW/NML1Lin/") # GSE132771

NML1 <- CreateSeuratObject(counts = NML1, project = "NML1", min.cells = 3, 
                           min.features = 200)

# 2. QC
NML1 <- PercentageFeatureSet(NML1, pattern = "^MT-", col.name = "percent.mt")

VlnPlot(NML1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

NML1 <- subset(NML1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & 
                 nCount_RNA < 10000 & percent.mt < 10)

# 3. Seurat standard workflow
NML1 <- NormalizeData(object = NML1)
NML1 <- FindVariableFeatures(object = NML1)
NML1 <- ScaleData(object = NML1)
NML1 <- RunPCA(object = NML1)
ElbowPlot(NML1)
NML1 <- FindNeighbors(object = NML1, dims = 1:20)
NML1 <- FindClusters(object = NML1)
NML1 <- RunUMAP(object = NML1, dims = 1:20)
DimPlot(NML1, reduction = 'umap', label = TRUE)

# 4. run scDblFinde to identify doublets, SingleCellExperiment object
Idents(NML1)

# Seurat v.5
sce <- scDblFinder(GetAssayData(NML1, layer = "counts"), clusters = Idents(NML1))

Seurat.cnt <- GetAssayData(NML1, layer = "counts")
View(Seurat.cnt)


# Seurat v.4
#sce <- scDblFinder(GetAssayData(NML1, slot= "counts"), clusters = Idents(NML1))

# 5. import the scDblFinder.class from sce object to the Seurat object
view(NML1@meta.data)

NML1$scDblFinder.class <- sce$scDblFinder.class

view(NML1@meta.data)

table(NML1@meta.data$scDblFinder.class)

# 6. visualize the doublets
DimPlot(NML1, reduction = 'umap', group.by = "scDblFinder.class")
DimPlot(NML1, reduction = 'umap', split.by = "scDblFinder.class")

Idents(NML1) <- NML1@meta.data$scDblFinder.class
Idents(NML1)
VlnPlot(NML1, features = c("nFeature_RNA", "nCount_RNA"), ncol=2)

# 7. subset the object: 
# NML1 <- subset(NML1, subset = scDblFinder.class == "singlet")
NML1 <- subset(NML1, idents = "singlet")


DimPlot(NML1, reduction = 'umap', label = TRUE)

# 8. re-cluster the cells
NML1 <- NormalizeData(object = NML1)
NML1 <- FindVariableFeatures(object = NML1)
NML1 <- ScaleData(object = NML1)
NML1 <- RunPCA(object = NML1)
NML1 <- FindNeighbors(object = NML1, dims = 1:20)
NML1 <- FindClusters(object = NML1)
NML1 <- RunUMAP(object = NML1, dims = 1:20)
DimPlot(NML1, reduction = 'umap', label = TRUE)

