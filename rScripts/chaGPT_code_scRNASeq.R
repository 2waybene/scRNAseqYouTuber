# Install required packages if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Seurat")

# Load required libraries
library(Seurat)

# Load the SingleCellExperiment package for reading the data
BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

# Load the GSE132771 dataset
library(GEOquery)
gse <- getGEO("GSE132771", GSEMatrix = TRUE)
sce <- SingleCellExperiment(assays = list(counts = as.matrix(exprs(gse))))

# Preprocess the data
sce <- sce[, rowSums(counts(sce)) > 0]  # Remove cells with zero counts
sce <- log1p(sce)  # Log-transform the data

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = as.matrix(sce), project = "GSE132771")

# Perform quality control
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# Cluster cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Visualize the results
plot1 <- DimPlot(seurat_obj, group.by = "orig.ident", label = TRUE)
plot2 <- FeaturePlot(seurat_obj, features = c("gene_of_interest"))

# Save plots (optional)
pdf("dimplot.pdf")
print(plot1)
dev.off()

pdf("featureplot.pdf")
print(plot2)
dev.off()
