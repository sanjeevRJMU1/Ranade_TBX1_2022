
## Script 1: Clustering

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)
library(URD)

# load neural crest subset object
Tbx1_CNC <- readRDS("~/scRNA-seq/2020_Tbx1/all_aggr_210429/Tbx1_CNC_annotated.RDS")

DimPlot(Tbx1_CNC, reduction = "umap", label = T)
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

# write seurat3ToURD, use SCT instead of RNA in assay 
# https://rdrr.io/github/farrellja/URD/src/R/import.R
# https://github.com/farrellja/URD/issues/29
seurat3ToURD <- function(seurat.object) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    
    # Copy over data (TODO: do we want raw RNA data or SCT data?)
    ds@logupx.data <- as(as.matrix(seurat.object@assays$SCT@data), "dgCMatrix")
    if(!any(dim(seurat.object@assays$SCT@counts) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$SCT@counts[rownames(seurat.object@assays$SCT@data), colnames(seurat.object@assays$SCT@data)]), "dgCMatrix")
    
    # Copy over metadata
    ## TO DO - grab kmeans clustering info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) { 
      get.data <- as.data.frame(seurat.object@assays$SCT@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) { 
      get.data <- as.data.frame(seurat.object@meta.data) 
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    
    # Copy over var.genes (removed if statement from original code, using scale data column to pull variable genes)
    ds@var.genes <- row.names(seurat.object@assays$SCT@scale.data)
    
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    }
    
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("pca" %in% names(seurat.object@reductions)) && !any(dim(Loadings(seurat.object, reduction = "pca")) == 0)) {
        ds@pca.load <- as.data.frame(Loadings(seurat.object, reduction = "pca"))
        ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
        ds@pca.sdev <- seurat.object@reductions$pca@stdev
        ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
}

Tbx1_CNC_URD <- seurat3ToURD(Tbx1_CNC)
Tbx1_CNC_URD
#URD object: 22938 genes x 10231 cells.

#calculate diffusion map
#Use sqroot of total number cells in dataset as KNN: sqroot of 35792 = 101
#run "##" line of code, after URD calculates global sigma, STOP and record sigma
Tbx1_CNC_URD <- calcDM(Tbx1_CNC_URD, knn=101, sigma.use=24)
#[1] "destiny determined an optimal global sigma of 23.943"

#plot diffusion maps using pairs of dimensions, for small datasets, structure of diff may already be apparent
plotDimArray(Tbx1_CNC_URD, reduction.use = "dm", dims.to.plot = 1:8, 
             outer.title = "Diffusion Map (Sigma 10, 189 NNs): Age", label="clusters", plot.title="", legend=T)
#define roots 
#here we chose cluster 1 from res 1.5
root.cells <- cellsInCluster(Tbx1_CNC_URD, "WT_E925", cluster=1)

# Then we run 'flood' simulations
combined.floods <- floodPseudotime(Tbx1_CNC_URD, 
                                   root.cells = root.cells, n=80, minimum.cells.flooded = 10, verbose=T)







