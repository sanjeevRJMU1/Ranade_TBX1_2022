###### CM subset at E11.5

#####################


setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

#
Tbx1_mesoderm_cnc <- readRDS(file = "Tbx1_mesoderm_cnc_annotated_210830.RDS")
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
table(Tbx1_mesoderm_cnc@active.ident,Tbx1_mesoderm_cnc$gem.group)

# subset to only E11.5
E115_subset <- SubsetData(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E115","KO_E115"))
table(E115_subset@active.ident,E115_subset$gem.group)
DimPlot(E115_subset, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()

# subset only CMs
E115_subset <- SubsetData(E115_subset, ident.use = "Cardiomyocyte")
table(E115_subset@active.ident,E115_subset$gem.group)
E115_subset@active.ident<-droplevels(E115_subset@active.ident)
table(E115_subset@active.ident,E115_subset$gem.group)
DimPlot(E115_subset, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
DimPlot(E115_subset, reduction = "umap", label = TRUE, pt.size = 0.1, group.by = "gem.group") + NoLegend()

## vln plot marker genes, DGE all clusters
VlnPlot(E115_subset, c("Tnnt2"), group.by = "gem.group", pt.size = 0.1)
VlnPlot(E115_subset, c("Rspo3"), group.by = "gem.group", pt.size = 0.1)
VlnPlot(E115_subset, c("Bmp4"), group.by = "gem.group", pt.size = 0.1)
VlnPlot(E115_subset, c("Vsnl1"), group.by = "gem.group", pt.size = 0.1)
VlnPlot(E115_subset, c("Nr2f2"), group.by = "gem.group", pt.size = 0.1)
VlnPlot(E115_subset, c("Irx4"), group.by = "gem.group", pt.size = 0.1)
VlnPlot(E115_subset, c("Bmp2"), group.by = "gem.group", pt.size = 0.1)


##
# recluster
E115_subset <- RunUMAP(E115_subset, dims = 1:20, reduction = "mnn", verbose = TRUE)
E115_subset <- FindNeighbors(E115_subset,reduction = "mnn", dims = 1:20, verbose = TRUE, force.recalc = T)
E115_subset <- FindClusters(E115_subset, verbose = TRUE, resolution = 0.2)
DimPlot(E115_subset, reduction = "umap", label = T)
E115_subset.markers <- FindAllMarkers(E115_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DimPlot(E115_subset, reduction = "umap", label = T, group.by = "gem.group")
table(E115_subset@active.ident, E115_subset$gem.group)

VlnPlot(E115_subset, c("Irx4"), pt.size = 0.1)
VlnPlot(E115_subset, c("Vsnl1"), pt.size = 0.1)
VlnPlot(E115_subset, c("Rspo3"), pt.size = 0.1)
VlnPlot(E115_subset, c("Bmp4"), pt.size = 0.1)
VlnPlot(E115_subset, c("Nr2f2"), pt.size = 0.1)
VlnPlot(E115_subset, c("Tbx1"), pt.size = 0.1)
VlnPlot(E115_subset, c("Tbx1"), pt.size = 0.1)
VlnPlot(E115_subset, c("Sema3c"), pt.size = 0.1)
VlnPlot(E115_subset, c("Tbx5"), pt.size = 0.1)
VlnPlot(E115_subset, c("Acta2"), pt.size = 0.1)

VlnPlot(E115_subset, c("Irx4","Vsnl1","Rspo3","Nr2f2","Tnc","Mybpc1"), pt.size = 0.1, ncol = 3)
FeaturePlot(E115_subset, c("Irx4","Vsnl1","Rspo3","Nr2f2","Tnc","Mybpc1"), ncol = 3, pt.size = 0.1)


### by clusters. i think cluster 0 = OFT, cluster 1 = atrium and cluster 2 = smooth muscle
### none of the clusters look ventricular due to lack of Irx4 expression

## Findmarkers
clust0 <- SubsetData(E115_subset, ident.use = 0)
clust0_WT <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("WT_E115"))
clust0_KO <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust0, cells= clust0_WT) <- "WT_E115"
Idents(object = clust0, cells= clust0_KO) <- "KO_E115"
clust0_findmark <- FindMarkers(clust0, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust0_findmark, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/findallmarkers_outputs/E115_CM_clust0_DGE_220606.csv")

clust1 <- SubsetData(E115_subset, ident.use = 1)
clust1_WT <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("WT_E115"))
clust1_KO <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust1, cells= clust1_WT) <- "WT_E115"
Idents(object = clust1, cells= clust1_KO) <- "KO_E115"
clust1_findmark <- FindMarkers(clust1, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust1_findmark, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/findallmarkers_outputs/E115_CM_clust1_DGE_220606.csv")


clust2 <- SubsetData(E115_subset, ident.use = 2)
clust2_WT <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("WT_E115"))
clust2_KO <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust2, cells= clust2_WT) <- "WT_E115"
Idents(object = clust2, cells= clust2_KO) <- "KO_E115"
clust2_findmark <- FindMarkers(clust2, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust2_findmark, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/findallmarkers_outputs/E115_CM_clust2_DGE_220606.csv")



## look at specific genes per sub-cluster to confirm
VlnPlot(E115_subset, c("Vsnl1"), group.by = "gem.group",idents = 0, pt.size = 0.1)
VlnPlot(E115_subset, c("Vsnl1"), group.by = "gem.group",idents = 1, pt.size = 0.1)
VlnPlot(E115_subset, c("Vsnl1"), group.by = "gem.group",idents = 2, pt.size = 0.1)

VlnPlot(E115_subset, c("Nr2f2"), group.by = "gem.group",idents = 0, pt.size = 0.1)
VlnPlot(E115_subset, c("Nr2f2"), group.by = "gem.group",idents = 1, pt.size = 0.1)
VlnPlot(E115_subset, c("Nr2f2"), group.by = "gem.group",idents = 2, pt.size = 0.1)

VlnPlot(E115_subset, c("Tbx5"), group.by = "gem.group",idents = 0, pt.size = 0.1)
VlnPlot(E115_subset, c("Tbx5"), group.by = "gem.group",idents = 1, pt.size = 0.1)
VlnPlot(E115_subset, c("Tbx5"), group.by = "gem.group",idents = 2, pt.size = 0.1)

####
# plot
my_order <- c("WT_E115","KO_E115")
clust0@active.ident <- factor(x = clust0@active.ident, levels = my_order)
table(clust0@active.ident)

VlnPlot(clust0, c("Rspo3","Bmp4","Sema3c","Isl1","Tbx5","Nr2f1","Nr2f2","Wnt2"), ncol = 4, pt.size = 0.1)

VlnPlot(clust1, c("Tbx5","Nr2f1","Nr2f2","Wnt2"), ncol = 4, pt.size = 0.1)




