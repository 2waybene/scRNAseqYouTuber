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



##===========================
##  Tutorial video 1
##===========================

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

##  faile to perform batch effect!!
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

##===========================
##  Tutorial video 2
##===========================



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






