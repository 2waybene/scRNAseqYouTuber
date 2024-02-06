
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(Seurat)

##  MyMacMonster
hdf5_obj <- Read10X_h5(filename = "~/SingleCellTutorial/dataDownload/HumanPBMCv3dot120K/20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5",
                       use.name=TRUE,
                       unique.features=TRUE)

pbmc.seurat <- CreateSeuratObject(counts  = hdf5_obj)

dim(pbmc.seurat@meta.data)

pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")

pbmc.seurat.filtered <- subset (pbmc.seurat, subset = nCount_RNA > 800 & 
                                  nFeature_RNA > 500 &
                                  mitoPercent < 10 )


dim(pbmc.seurat.filtered)

pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

View(pbmc.seurat.filtered)
DimPlot(pbmc.seurat.filtered, reduction = "umap")

saveRDS(pbmc.seurat.filtered, file = "~/SingleCellTutorial/singleR/pbmc_filtered_processed.RDS")
pbmc.seurat.filtered <- readRDS("~/SingleCellTutorial/singleR/pbmc_filtered_processed.RDS")

ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))
     
pbmc_count <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')
pred <- SingleR (test = pbmc_count,
        ref=ref,
        labels = ref$label.main)

pred

pbmc.seurat.filtered$singleR.labels <- pred$labels[match (rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
DimPlot (pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')

library(viridis)
grid::current.viewport()
plotScoreHeatmap(pred)
plotDeltaDistribution(pred)
tab <- table (Assigned=pred$labels, Clusters = pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), colorRampPalette(c('white','blue'))(10))
