##========================================
## Convert Seurat Objects into H5AD Files
## -- video 20
##========================================
library(Seurat)
library(tidyverse)
library(Matrix)


NML <- readRDS( "/Users/jyli/SingleCellTutorial/workingData/MergeDNML_integrated.RDS")
DimPlot (NML, reduction = 'umap', label = TRUE)
View(NML@meta.data)


##  Seurat v4
counts_matrix <- GetAssayData(NML, assay = 'RNA', slot = 'counts')

##  Seurat v5
counts_matrix <- GetAssayData(NML, assay = 'RNA', layer  = 'counts')

View(counts_matrix)
writeMM(counts_matrix, file = paste0(file = "/Users/jyli/SingleCellTutorial/SeuratToH5AD/matrix.mtx"))


write.csv (NML@reductions$pca@cell.embeddings,
           file = "/Users/jyli/SingleCellTutorial/SeuratToH5AD/pca.csv", quote = F, row.names = F)

write.table (data.frame('gene'=rownames(counts_matrix)),
             file = "/Users/jyli/SingleCellTutorial/SeuratToH5AD/gene_names.csv",
             quote=F, row.names = F, col.names = F)

View(NML@meta.data)
NML$barcode=colnames(NML)
NML$UMAP_1 <- NML@reductions$umap@cell.embeddings[,1]
NML$UMAP_2 <- NML@reductions$umap@cell.embeddings[,2]

write.csv(NML@meta.data, file = "/Users/jyli/SingleCellTutorial/SeuratToH5AD/metadata.csv",
          quote = F, row.names = F)


## continue with python code
##  /Users/jyli/SingleCellTutorial/SeuratToH5AD/convertSeurat.ipynb
## In there, I am having trouble with package scvelo and have to choose a temporary way around it
## ## https://stackoverflow.com/questions/53014306/error-15-initializing-libiomp5-dylib-but-found-libiomp5-dylib-already-initial
## have to temporarily use the following two lines of code to get around with this problem

## import os
## os.environ['KMP_DUPLICATE_LIB_OK']='True'





