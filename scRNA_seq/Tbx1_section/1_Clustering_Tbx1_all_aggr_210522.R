####### Tbx1 all aggr: E9.5, E10.5, E11.5

## Script 1: Clustering

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

### 1. Initialize Seurat object
Tbx1_all_aggr.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
Tbx1_all_aggr_sct <- CreateSeuratObject(counts = Tbx1_all_aggr.data, project = "Tbx1_all_aggr", min.cells = 3, min.features = 200)
Tbx1_all_aggr_sct
# An object of class Seurat 
# 27627 features across 104897 samples within 1 assay 
# Active assay: RNA (27627 features, 0 variable features)

### 2. Add mito percentage and gem group to meta data
# mito percentage 
Tbx1_all_aggr_sct[["percent.mt"]] <- PercentageFeatureSet(Tbx1_all_aggr_sct, pattern = "^mt-")

# gem group naming
classification.vec <- as.numeric(gsub(".*-","", (colnames(x = Tbx1_all_aggr_sct))))
names(classification.vec) <- colnames(x = Tbx1_all_aggr_sct)
classification.vec[classification.vec=="1"] <- "WT_E925"
classification.vec[classification.vec=="2"] <- "KO_E925"
classification.vec[classification.vec=="3"] <- "WT_E925"
classification.vec[classification.vec=="4"] <- "KO_E925"
classification.vec[classification.vec=="5"] <- "KO_E105"
classification.vec[classification.vec=="6"] <- "KO_E105"
classification.vec[classification.vec=="7"] <- "WT_E105"
classification.vec[classification.vec=="8"] <- "WT_E105"
classification.vec[classification.vec=="9"] <- "KO_E115"
classification.vec[classification.vec=="10"] <- "WT_E115"
classification.vec[classification.vec=="11"] <- "WT_E115"
classification.vec[classification.vec=="12"] <- "KO_E115"
classification.vec[classification.vec=="13"] <- "KO_E115"
Tbx1_all_aggr_sct$"gem.group" <- classification.vec
head(Tbx1_all_aggr_sct@meta.data)
tail(Tbx1_all_aggr_sct@meta.data)

### 3. Filter poor QC cells --> This is a critical step!!!
#VlnPlot(Tbx1_all_aggr_sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filter cells by percent.mito
Tbx1_all_aggr_sct <- subset(Tbx1_all_aggr_sct, subset = percent.mt < 15)

# filter cells by nFeature_RNA 
Tbx1_all_aggr_sct <- subset(Tbx1_all_aggr_sct, subset = nFeature_RNA > 3500 & nFeature_RNA < 11000)

# filter cells by nCount_RNA (UMI counts)
Tbx1_all_aggr_sct <- subset(Tbx1_all_aggr_sct, subset = nCount_RNA < 100000)

# re-run metrics plot
#VlnPlot(Tbx1_all_aggr_sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# see how many cells are left
Tbx1_all_aggr_sct
# An object of class Seurat 
# 27627 features across 57777 samples within 1 assay 
# Active assay: RNA (27627 features, 0 variable features)

### 4. Assign cell cycle scores
cc.genes <- readLines(con = "/Users/sranade/scRNA-seq/cell_cycle_vignette_files/cell_cycle_ssr.txt")
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:45]
g2m.genes <- cc.genes[46:100]

# assign cell cycle scores
Tbx1_all_aggr_sct <- CellCycleScoring(Tbx1_all_aggr_sct, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = Tbx1_all_aggr_sct@meta.data)

### 5. Run SCTransform
Tbx1_all_aggr_sct <- SCTransform(Tbx1_all_aggr_sct, vars.to.regress = c("S.Score", "G2M.Score"), verbose = TRUE)

### 6. Run PCA
Tbx1_all_aggr_sct <- RunPCA(Tbx1_all_aggr_sct, verbose = TRUE, assay = "SCT")

### 7. Run batch correction, MNN
Tbx1_all_aggr_sct <- RunFastMNN(object.list = SplitObject(Tbx1_all_aggr_sct, split.by = "gem.group"))

### 8a. Cluster cells -- first round
Tbx1_all_aggr_sct <- RunUMAP(Tbx1_all_aggr_sct, dims = 1:30, reduction = "mnn", verbose = TRUE)
Tbx1_all_aggr_sct <- FindNeighbors(Tbx1_all_aggr_sct,reduction = "mnn", dims = 1:30, verbose = TRUE)
Tbx1_all_aggr_sct <- FindClusters(Tbx1_all_aggr_sct, verbose = TRUE, resolution = 1)
p1 <- DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = T)
p2 <- DimPlot(Tbx1_all_aggr_sct, reduction = "umap", group.by = "gem.group")
CombinePlots(plots = list(p1, p2))


### 8b. QC checks
# VlnPlot(Tbx1_all_aggr_sct, c("nFeature_RNA"), pt.size = 0.1, ncol = 1)
# Tbx1_all_aggr_sct <- subset(Tbx1_all_aggr_sct, subset = nFeature_RNA > 4000 & nFeature_RNA < 9000)
# VlnPlot(Tbx1_all_aggr_sct, c("nCount_RNA"), pt.size = 0.1, ncol = 1)
# Tbx1_all_aggr_sct <- subset(Tbx1_all_aggr_sct, subset = nCount_RNA < 80000)

# ### 8c. Re-Cluster cells
# Tbx1_all_aggr_sct <- RunUMAP(Tbx1_all_aggr_sct, dims = 2:30, reduction = "mnn", verbose = TRUE)
# Tbx1_all_aggr_sct <- FindNeighbors(Tbx1_all_aggr_sct,reduction = "mnn", dims = 2:30, verbose = TRUE, force.recalc = T)
# Tbx1_all_aggr_sct <- FindClusters(Tbx1_all_aggr_sct, verbose = TRUE, resolution = 0.7)
# DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = T)
# 
# table(Tbx1_all_aggr_sct@active.ident,Tbx1_all_aggr_sct$gem.group)


### 8d. Markers and save object
Tbx1_all_aggr_sct.markers <- FindAllMarkers(Tbx1_all_aggr_sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_all_aggr_sct.markers, file = "Tbx1_all_aggr_sct_all-markers.csv")
saveRDS(Tbx1_all_aggr_sct, file = "Tbx1_all_aggr_sct.RDS")
Tbx1_all_aggr_sct
# An object of class Seurat 
# 54282 features across 53805 samples within 2 assays 
# Active assay: SCT (26655 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

####

Tbx1_trim <- SubsetData(Tbx1_all_aggr_sct, ident.use = c(0,22,4,29,18,11,25,9,7,1,14,39,36,28,13,8,2,31,12,20,10,21,5))
#
Tbx1_trim <- RunUMAP(Tbx1_trim, dims = 1:20, reduction = "mnn", verbose = TRUE)
Tbx1_trim <- FindNeighbors(Tbx1_trim,reduction = "mnn", dims = 1:20, verbose = TRUE, force.recalc = T)
Tbx1_trim <- FindClusters(Tbx1_trim, verbose = TRUE, resolution = 0.5)
DimPlot(Tbx1_trim, reduction = "umap", label = T)

VlnPlot(Tbx1_trim, c("nFeature_RNA"), pt.size = 0.1, ncol = 1)
Tbx1_trim <- SubsetData(Tbx1_trim, ident.remove = 15)


#
Tbx1_trim <- RunUMAP(Tbx1_trim, dims = 2:25, reduction = "mnn", verbose = TRUE)
Tbx1_trim <- FindNeighbors(Tbx1_trim,reduction = "mnn", dims = 2:25, verbose = TRUE, force.recalc = T)
Tbx1_trim <- FindClusters(Tbx1_trim, verbose = TRUE, resolution = 0.6)
DimPlot(Tbx1_trim, reduction = "umap", label = T)

####
Tbx1_trim.markers <- FindAllMarkers(Tbx1_trim, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_trim.markers, file = "Tbx1_trim_all-markers.csv")
saveRDS(Tbx1_trim, file = "Tbx1_trim.RDS")

table(Tbx1_trim@active.ident,Tbx1_trim$gem.group)


####
Tbx1_CNC <- SubsetData(Tbx1_trim, ident.use = c(19,11,2,9,10,6))
#
Tbx1_CNC <- RunUMAP(Tbx1_CNC, dims = 2:25, reduction = "mnn", verbose = TRUE)
Tbx1_CNC <- FindNeighbors(Tbx1_CNC,reduction = "mnn", dims = 2:25, verbose = TRUE, force.recalc = T)
Tbx1_CNC <- FindClusters(Tbx1_CNC, verbose = TRUE, resolution = 0.3)
DimPlot(Tbx1_CNC, reduction = "umap", label = T)
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)
###
Tbx1_CNC.markers <- FindAllMarkers(Tbx1_CNC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_CNC.markers, file = "Tbx1_CNC_all-markers.csv")
saveRDS(Tbx1_CNC, file = "Tbx1_CNC.RDS")
####
# re-run SCT on Tbx1_CNC
Tbx1_CNC <- CellCycleScoring(Tbx1_CNC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = Tbx1_CNC@meta.data)
Tbx1_CNC <- SCTransform(Tbx1_CNC, vars.to.regress = c("S.Score", "G2M.Score"), verbose = TRUE)
Tbx1_CNC <- RunPCA(Tbx1_CNC, verbose = TRUE, assay = "SCT")
Tbx1_CNC <- RunFastMNN(object.list = SplitObject(Tbx1_CNC, split.by = "gem.group"))

#
Tbx1_CNC <- RunUMAP(Tbx1_CNC, dims = 2:25, reduction = "mnn", verbose = TRUE)
Tbx1_CNC <- FindNeighbors(Tbx1_CNC,reduction = "mnn", dims = 2:25, verbose = TRUE, force.recalc = T)
Tbx1_CNC <- FindClusters(Tbx1_CNC, verbose = TRUE, resolution = 0.3)
p1 <- DimPlot(Tbx1_CNC, reduction = "umap", label = T)
p2 <- DimPlot(Tbx1_CNC, reduction = "umap", group.by = "gem.group")
CombinePlots(plots = list(p1, p2))


Tbx1_CNC <- RunUMAP(Tbx1_CNC, dims = 1:12, reduction = "mnn", verbose = TRUE)
Tbx1_CNC <- FindNeighbors(Tbx1_CNC,reduction = "mnn", dims = 1:12, verbose = TRUE, force.recalc = T)
Tbx1_CNC <- FindClusters(Tbx1_CNC, verbose = TRUE, resolution = 0.3)

DimPlot(Tbx1_CNC, reduction = "umap", label = T)

table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)
Tbx1_CNC <-SubsetData(Tbx1_CNC, ident.remove = 6)

Tbx1_CNC <- RunUMAP(Tbx1_CNC, dims = 1:12, reduction = "mnn", verbose = TRUE)
Tbx1_CNC <- FindNeighbors(Tbx1_CNC,reduction = "mnn", dims = 1:12, verbose = TRUE, force.recalc = T)
Tbx1_CNC <- FindClusters(Tbx1_CNC, verbose = TRUE, resolution = 0.3)
DimPlot(Tbx1_CNC, reduction = "umap", label = T)
VlnPlot(Tbx1_CNC, c("nFeature_RNA"), pt.size = 0.1, ncol = 1)

Tbx1_CNC.markers <- FindAllMarkers(Tbx1_CNC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_CNC.markers, file = "Tbx1_CNC_all-markers.csv")
saveRDS(Tbx1_CNC, file = "Tbx1_CNC.RDS")

### 
# re-load
Tbx1_all_aggr_sct <-readRDS(file = "Tbx1_all_aggr_sct.RDS")
DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = T)

# 
Tbx1_CNC <-readRDS(file = "Tbx1_CNC.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = T)


