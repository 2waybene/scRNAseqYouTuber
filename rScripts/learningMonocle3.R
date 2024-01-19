##=======================================================================================
##  learningMonocle3.R
##  https://www.youtube.com/watch?v=YQqVsXdwFNU&list=PLOLdjuxsfI4P48Sl_N0hOFs7WYOv23yUr
##  https://cole-trapnell-lab.github.io/monocle3/docs/installation/
##=======================================================================================

##  dependencies for monocle3
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                       'terra', 'ggrastr'))



##============================================
##  Tutorial video 1
##  Read in raw data and create cds object
##===========================================

library(monocle3)
library(tidyverse)
library(patchwork)
set.seed(1234)

cds1 <- load_mm_data(mat_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF1/GSM3891626_IPF1_Lin_matrix.mtx",
                     feature_anno_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF1/GSM3891626_IPF1_Lin_genes.tsv",
                     cell_anno_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF1/GSM3891626_IPF1_Lin_barcodes.tsv")

colData(cds1)
colData(cds1)$samples <- paste("IPF1")

cds2 <- load_mm_data(mat_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF2/GSM3891628_IPF2_Lin_matrix.mtx",
                     feature_anno_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF2/GSM3891628_IPF2_Lin_genes.tsv",
                     cell_anno_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF2/GSM3891628_IPF2_Lin_barcodes.tsv")
colData(cds2)$samples <- paste("IPF2")
colData(cds2)

cds3 <- load_mm_data(mat_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF3/GSM3891630_IPF3_Lin_matrix.mtx",
                     feature_anno_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF3/GSM3891630_IPF3_Lin_features.tsv",
                     cell_anno_path = "/Users/li11/myGit/SingleCellDataDownload/GSE132771/IPF3/GSM3891630_IPF3_Lin_barcodes.tsv")
colData(cds3)$samples <- paste("IPF3")
colData(cds3)

View(cds1)
cds <- combine_cds(list(cds1,cds2,cds3))

rm (cds1,cds2,cds3)

fData(cds)
names(fData(cds)) <- "gene_short_name"
fData(cds)


cds <- preprocess_cds(cds,num_dim = 100)

plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)
colData(cds)
plot_cells(cds, color_cells_by = "samples", group_label_size = 6)

##  failed to perform batch effect!!
##  https://github.com/cole-trapnell-lab/monocle3/issues/690
cds <- align_cds(cds, alignment_group = "samples", useNames=TRUE)
cds <- reduce_dimension(cds) 
plot_cells(cds, color_cells_by = "samples", group_label_size = 6)


cds <- cluster_cells(cds, resolution = 1e-5)

plot_cells(cds, group_label_size = 6)

marker_test_res <- top_markers(cds, group_cells_by = "cluster", verbose = T)

top_specific_makers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, marker_score) #3

top_specific_makers_ids <- unique(top_specific_makers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_makers_ids,
                    group_cells_by = "cluster",
                    ordering_type = "maximal_on_diag",
                    max.size = 3)

P1 <- plot_cells(cds, color_cells_by = "cluster", group_cells_by = "cluster", group_label_size = 6)
P2 <- plot_genes_by_group(cds,
                          top_specific_makers_ids,
                          group_cells_by = "cluster",
                          max.size = 3)
P1|P2
clusters(cds)
colData(cds)
colData(cds)$assigned_cell_type <- as.character(clusters(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1" = "ALV", "2" = "Perc", "3"= "SMCs", "4" = "Airway",
                                                 "5" = "Mes", "6" = "IM", "7" = "AT2", "8" = "Basal",
                                                 "9" = "_ALV", "10" = "_Perc", "11"= "_SMCs", "12" = "else")
colData(cds)
plot_cells(cds, group_cells_by = "cluster", color_cells_by = "assigned_cell_type",
           group_label_size = 6)
saveRDS (cds, file = "/Users/li11/myGit/SingleCellDataDownload/processData/cds.RDS")
   
rm (cds)

save_monocle_objects(cds=cds, directory_path='/Users/li11/myGit/SingleCellDataDownload/processData/')

## cds_object.rds  (full_cds  RDS  from  cell_data_set)
## rdd_umap_transform_model_umap.idx  (UMAP  UMAP_NN_index  from  reduce_dimension)

cds <- load_monocle_objects(directory_path='/Users/li11/myGit/SingleCellDataDownload/processData/')



##=========================================================================================
## my comments
## I liked your videos a lot and deeply appreciated for your enthusiasm! 
## I found the monocle3 more challenging!! The installation has problem,  
## and functions have changed. One hurdle that I can't overcome is the align_cds, 
## it does not seem work for me. I used your data but failed to "integrate" 
## three samples well. Therefore, the following analysis does not make a whole lot sense.
##=========================================================================================

##===============================
##  Tutorial video 2
##  Constructing Trajectories
##==============================



library(monocle3)

# 1. Load the analysed data
cds <- readRDS("/Users/li11/myGit/SingleCellDataDownload/processData/cds.RDS")
plot_cells(cds, group_cells_by="cluster", color_cells_by="assigned_cell_type",
           group_label_size = 6)

# 2. subset SMCs, Pericytes and fibroblasts
cds_subset <- choose_cells(cds)
cds_subset <- reduce_dimension(cds_subset)
cds_subset <- cluster_cells(cds_subset, resolution = 1e-5)
plot_cells(cds_subset, color_cells_by="assigned_cell_type", group_label_size = 6)

saveRDS(cds_subset, file="/Users/li11/myGit/SingleCellDataDownload/processData/cds_subset.RDS")

plot_cells(cds_subset, genes=c("PDGFRB", "ACTA2", "RGS5", "CTHRC1"))
plot_cells(cds_subset, genes=c("COL1A1", "COL3A1", "CTHRC1", "POSTN"))
plot_cells(cds_subset, genes=c("CD34", "PI16", "SCARA5", "MFAP5"))

# 3. Constructing Single Cell Trajectories 
cds_subset <- learn_graph(cds_subset)

plot_cells(cds_subset, color_cells_by = "assigned_cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size = 6)

cds_subset <- order_cells(cds_subset)

plot_cells(cds_subset, color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size= 3)



##=========================================================
##  Tutorial video 3
##  Convert Seurat object to cds for trajectory analysis
##=========================================================

library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library(monocle3)

# 1. Analyze scRNA-seq data using Seurat------
# Read
IPF_I<- Read10X(data.dir = "/Users/li11/myGit/SingleCellDataDownload/GSE132771_RAW/IPF1Lin/")
IPF_II<- Read10X(data.dir = "/Users/li11/myGit/SingleCellDataDownload/GSE132771_RAW/IPF2Lin/")
IPF_III<- Read10X(data.dir = "/Users/li11/myGit/SingleCellDataDownload/GSE132771_RAW/IPF3Lin/")

# Create Seurat Object
IPF_I <- CreateSeuratObject(counts = IPF_I, project = "IPF_I", 
                            min.cells = 3, min.features = 200) 
IPF_I <- PercentageFeatureSet(IPF_I, pattern = "^MT-", col.name = "percent.mt")

IPF_II <- CreateSeuratObject(counts = IPF_II, project = "IPF_II", 
                             min.cells = 3, min.features = 200)
IPF_II <- PercentageFeatureSet(IPF_II, pattern = "^MT-", col.name = "percent.mt")

IPF_III <- CreateSeuratObject(counts = IPF_III, project = "IPF_III", 
                              min.cells = 3, min.features = 200)
IPF_III <- PercentageFeatureSet(IPF_III,pattern = "^MT-",col.name = "percent.mt")


IPF <- merge(IPF_I, y = c(IPF_II, IPF_III))
VlnPlot(IPF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
IPF <- subset(IPF, nFeature_RNA<5000 & nCount_RNA<20000 & percent.mt<10)

View(IPF@meta.data)

IPF.list <- SplitObject(IPF, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
IPF.list <- lapply(X = IPF.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = IPF.list)
IPF.anchors <- FindIntegrationAnchors(object.list = IPF.list, 
                                      anchor.features = features)

# creates an 'integrated' data assay
IPF.combined <- IntegrateData(anchorset = IPF.anchors)

# specify that we will perform downstream analysis on the integrated data
DefaultAssay(IPF.combined) <- "integrated"

# Run the standard workflow for cell clustering
IPF.combined <- ScaleData(IPF.combined, verbose = FALSE)
IPF.combined <- RunPCA(IPF.combined, npcs = 30, verbose = FALSE)
IPF.combined <- RunUMAP(IPF.combined, reduction = "pca", dims = 1:30)
IPF.combined <- FindNeighbors(IPF.combined, reduction = "pca", dims = 1:30)
IPF.combined <- FindClusters(IPF.combined, resolution = 0.5)

# Visualization
DimPlot(IPF.combined, reduction = "umap", label = TRUE)

# four markers
FeaturePlot(IPF.combined, features=c("COL1A2", "LUM", "PDGFRA","PDGFRB"),
            label = T, cols = c('lightgrey', 'blue'))
FeaturePlot(IPF.combined, features=c("EPCAM", "SFTPC", "AGER", "SCGB1A1"),
            label = T, cols = c('lightgrey', 'blue'))
FeaturePlot(IPF.combined, features=c("CLDN5", "CCL21", "PECAM1", "EMCN"),
            label = T, cols = c('lightgrey', 'blue'))
FeaturePlot(IPF.combined, features=c("PTPRC", "CD52", "AIF1", "TRBC2"),
            label = T,cols = c('lightgrey', 'blue'))
FeaturePlot(IPF.combined, features=c("MSLN", "CALB2", "HP", "PRG4"),
            label = T, cols = c('lightgrey', 'blue'))

IPF_Mesenchymal <- subset(x = IPF.combined, idents = c("0", "1", "2", 
                                                       "3", "4", "5", "6", "7","8", "9"))
DimPlot(IPF_Mesenchymal, reduction = "umap", label = TRUE)

# Run the standard workflow for cell clustering
IPF_Mesenchymal <- ScaleData(IPF_Mesenchymal, verbose = FALSE)
IPF_Mesenchymal <- RunPCA(IPF_Mesenchymal, npcs = 30, verbose = FALSE)
IPF_Mesenchymal <- RunUMAP(IPF_Mesenchymal, reduction = "pca", dims = 1:30)
IPF_Mesenchymal <- FindNeighbors(IPF_Mesenchymal, reduction = "pca", dims = 1:30)
IPF_Mesenchymal <- FindClusters(IPF_Mesenchymal, resolution = 0.3)

# Visualization
DimPlot(IPF_Mesenchymal, reduction = "umap", label = TRUE)

IPF_Mesenchymal <- subset(x = IPF_Mesenchymal, idents = c("0", "1", "2", 
                                                          "3", "4", "5", "6"))

FeaturePlot(IPF_Mesenchymal, features=c("RGS5", "ACTA2", "LUM", "ACTG2"),
            label = T, cols = c('lightgrey', 'blue'))
FeaturePlot(IPF_Mesenchymal, features=c("CTHRC1", "COL1A1", "COL3A1", "POSTN"),
            label = T, cols = c('lightgrey', 'blue'))

IPF_Mesenchymal <- RenameIdents(IPF_Mesenchymal,
                                `0` = "Pericytes", `4` = "Pericytes", `1` = "SMCs",
                                `3` = "MyoF", `2` = "AlvF", `5` = "AlvF", `6` = "AirwayF")

DimPlot(IPF_Mesenchymal , reduction = "umap", label = TRUE)

saveRDS(IPF_Mesenchymal, 
       file="/Users/li11/myGit/SingleCellDataDownload/processData/IPF_Mesenchymal.RDS")

rm("features","IPF.anchors", "IPF.combined","IPF", "IPF.list","IPF_I", 
   "IPF_II", "IPF_III")

# 2. Convert to Seurat object to cell_data_set object ------
cds <- as.cell_data_set(IPF_Mesenchymal, group.by='ident') # SeuratWrappers
DefaultAssay(IPF_Mesenchymal)
DefaultAssay(IPF_Mesenchymal) <- "RNA"

# plot
colData(cds)
plot_cells(cds, color_cells_by = 'seurat_clusters',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")

plot_cells(cds, color_cells_by = 'ident', 
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")

##==============================================================================
# Learn the trajectory graph
# Failed because: Error: names(cds@clusters[[reduction_method]]$partitions) == NULL
# Temporarily comment the following line

##  cds <- learn_graph(cds, use_partition = FALSE)
##==============================================================================

##  we do NOT need to re-run UMAP here, just cluster_cells

cds <- cluster_cells(cds, resolution = 1e-5)

plot_cells(cds, color_cells_by = 'ident', 
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")

# 3. trajectory analysis------
cds <- learn_graph(cds, use_partition = FALSE)

# plot
plot_cells(cds, color_cells_by = 'ident', 
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")

##  Error again here!!
## Error in plot_cells(cds, genes = c("CD34", "PI16", "SCARA5", "MFAP5")) : 
##   None of the provided genes were found in the cds
plot_cells(cds, genes=c("CD34", "PI16", "SCARA5", "MFAP5"))

# Create gene annotation file
fData(cds) # gene_annotation
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
fData(cds)

# Order the cells by pseudotime 
cds <- order_cells(cds)

#plot_cells in pseudotime
plot_cells(cds, color_cells_by = 'pseudotime')



