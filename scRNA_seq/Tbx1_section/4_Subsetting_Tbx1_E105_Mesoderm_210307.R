####### E10.5 Tbx1 scRNA-seq
## SR38_1 and SR38_2 = KO
## SR38_3 and SR38_4 = WT
## Script 4: Cluster subsetting


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/E105_SR40")
#
Mesoderm_annotated <- readRDS("Mesoderm_annotated.RDS")
DimPlot(Mesoderm_annotated, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()



# Mesoderm Subset --> not including phox2b/myt1l clusters. do that separately
Mesoderm_subset <- SubsetData(Mesoderm_annotated, ident.use = c(20,5,24,8,26,4,13,18,1,3))
Mesoderm_subset <- RunUMAP(Mesoderm_subset, dims = 1:30, reduction = "mnn", verbose = TRUE)
Mesoderm_subset <- FindNeighbors(Mesoderm_subset,reduction = "mnn", dims = 1:30, verbose = TRUE, force.recalc = T)
Mesoderm_subset <- FindClusters(Mesoderm_subset, verbose = TRUE, resolution = 1)
#p3 <- DimPlot(Mesoderm_subset, reduction = "umap", label = T, pt.size = 0.2)
#p4 <- DimPlot(Mesoderm_subset, reduction = "umap", group.by = "gem.group", pt.size = 0.2)
#CombinePlots(plots = list(p3, p4))
DimPlot(Mesoderm_subset, reduction = "umap", label = T, pt.size = 0.2)
VlnPlot(Mesoderm_subset, c("nFeature_RNA","nCount_RNA"), pt.size = 0.1, ncol = 1)

Mesoderm_subset <- RunUMAP(Mesoderm_subset, dims = 1:12, reduction = "mnn", verbose = TRUE)
Mesoderm_subset <- FindNeighbors(Mesoderm_subset,reduction = "mnn", dims = 1:12, verbose = TRUE, force.recalc = T)
Mesoderm_subset <- FindClusters(Mesoderm_subset, verbose = TRUE, resolution = 0.5)
DimPlot(Mesoderm_subset, reduction = "umap", label = T, pt.size = 0.3)

DimPlot(Mesoderm_subset, reduction = "umap", label = F, pt.size = 0.3, group.by = "gem.group")

table(Mesoderm_subset@active.ident, Mesoderm_subset$gem.group)
Mesoderm_subset.markers <- FindAllMarkers(Mesoderm_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Mesoderm_subset.markers, file = "Mesoderm_subset.markers.csv")
saveRDS(Mesoderm_subset, file = "Mesoderm_subset.RDS")
Mesoderm_subset <- readRDS("~/scRNA-seq/2020_Tbx1/E105_SR40/Mesoderm_subset.RDS")

### DGE

### DGE
# subset each cluster
table(Mesoderm_subset@active.ident, Mesoderm_subset$gem.group)
c0 <- SubsetData(Mesoderm_subset, ident.use = 0)
c1 <- SubsetData(Mesoderm_subset, ident.use = 1)
c2 <- SubsetData(Mesoderm_subset, ident.use = 2)
c3 <- SubsetData(Mesoderm_subset, ident.use = 3)
c4 <- SubsetData(Mesoderm_subset, ident.use = 4)
c5 <- SubsetData(Mesoderm_subset, ident.use = 5)
c6 <- SubsetData(Mesoderm_subset, ident.use = 6)
c7 <- SubsetData(Mesoderm_subset, ident.use = 7)
c8 <- SubsetData(Mesoderm_subset, ident.use = 8)
c9 <- SubsetData(Mesoderm_subset, ident.use = 9)
c10 <- SubsetData(Mesoderm_subset, ident.use = 10)
c11 <- SubsetData(Mesoderm_subset, ident.use = 11)


# 
c0_WT <- OldWhichCells(c0, subset.name = "gem.group", accept.value = c("WT"))
c0_KO <- OldWhichCells(c0, subset.name = "gem.group", accept.value = c("KO"))

c1_WT <- OldWhichCells(c1, subset.name = "gem.group", accept.value = c("WT"))
c1_KO <- OldWhichCells(c1, subset.name = "gem.group", accept.value = c("KO"))

c2_WT <- OldWhichCells(c2, subset.name = "gem.group", accept.value = c("WT"))
c2_KO <- OldWhichCells(c2, subset.name = "gem.group", accept.value = c("KO"))

c3_WT <- OldWhichCells(c3, subset.name = "gem.group", accept.value = c("WT"))
c3_KO <- OldWhichCells(c3, subset.name = "gem.group", accept.value = c("KO"))

c4_WT <- OldWhichCells(c4, subset.name = "gem.group", accept.value = c("WT"))
c4_KO <- OldWhichCells(c4, subset.name = "gem.group", accept.value = c("KO"))

c5_WT <- OldWhichCells(c5, subset.name = "gem.group", accept.value = c("WT"))
c5_KO <- OldWhichCells(c5, subset.name = "gem.group", accept.value = c("KO"))

c6_WT <- OldWhichCells(c6, subset.name = "gem.group", accept.value = c("WT"))
c6_KO <- OldWhichCells(c6, subset.name = "gem.group", accept.value = c("KO"))

c7_WT <- OldWhichCells(c7, subset.name = "gem.group", accept.value = c("WT"))
c7_KO <- OldWhichCells(c7, subset.name = "gem.group", accept.value = c("KO"))

c8_WT <- OldWhichCells(c8, subset.name = "gem.group", accept.value = c("WT"))
c8_KO <- OldWhichCells(c8, subset.name = "gem.group", accept.value = c("KO"))

c9_WT <- OldWhichCells(c9, subset.name = "gem.group", accept.value = c("WT"))
c9_KO <- OldWhichCells(c9, subset.name = "gem.group", accept.value = c("KO"))

c10_WT <- OldWhichCells(c10, subset.name = "gem.group", accept.value = c("WT"))
c10_KO <- OldWhichCells(c10, subset.name = "gem.group", accept.value = c("KO"))

c11_WT <- OldWhichCells(c11, subset.name = "gem.group", accept.value = c("WT"))
c11_KO <- OldWhichCells(c11, subset.name = "gem.group", accept.value = c("KO"))



# b.
Idents(object = c0, cells= c0_WT) <- "c0_WT"
Idents(object = c0, cells= c0_KO) <- "c0_KO"

Idents(object = c1, cells= c1_WT) <- "c1_WT"
Idents(object = c1, cells= c1_KO) <- "c1_KO"

Idents(object = c2, cells= c2_WT) <- "c2_WT"
Idents(object = c2, cells= c2_KO) <- "c2_KO"

Idents(object = c3, cells= c3_WT) <- "c3_WT"
Idents(object = c3, cells= c3_KO) <- "c3_KO"

Idents(object = c4, cells= c4_WT) <- "c4_WT"
Idents(object = c4, cells= c4_KO) <- "c4_KO"

Idents(object = c5, cells= c5_WT) <- "c5_WT"
Idents(object = c5, cells= c5_KO) <- "c5_KO"

Idents(object = c6, cells= c6_WT) <- "c6_WT"
Idents(object = c6, cells= c6_KO) <- "c6_KO"

Idents(object = c7, cells= c7_WT) <- "c7_WT"
Idents(object = c7, cells= c7_KO) <- "c7_KO"

Idents(object = c8, cells= c8_WT) <- "c8_WT"
Idents(object = c8, cells= c8_KO) <- "c8_KO"

Idents(object = c9, cells= c9_WT) <- "c9_WT"
Idents(object = c9, cells= c9_KO) <- "c9_KO"

Idents(object = c10, cells= c10_WT) <- "c10_WT"
Idents(object = c10, cells= c10_KO) <- "c10_KO"

Idents(object = c11, cells= c11_WT) <- "c11_WT"
Idents(object = c11, cells= c11_KO) <- "c11_KO"

# 3. Run Findmarkers -- normally you set the lower cell # but since its n=1, just run it
c0_findmark <- FindMarkers(c0, ident.1 = "c0_WT", ident.2 = "c0_KO",max.cells.per.ident = 563,logfc.threshold = 0.3)
c1_findmark <- FindMarkers(c1, ident.1 = "c1_WT", ident.2 = "c1_KO",max.cells.per.ident = 243,logfc.threshold = 0.3)
c2_findmark <- FindMarkers(c2, ident.1 = "c2_WT", ident.2 = "c2_KO",max.cells.per.ident = 321,logfc.threshold = 0.3)
c3_findmark <- FindMarkers(c3, ident.1 = "c3_WT", ident.2 = "c3_KO",max.cells.per.ident = 278,logfc.threshold = 0.3)
c4_findmark <- FindMarkers(c4, ident.1 = "c4_WT", ident.2 = "c4_KO",max.cells.per.ident = 302,logfc.threshold = 0.3)
c5_findmark <- FindMarkers(c5, ident.1 = "c5_WT", ident.2 = "c5_KO",max.cells.per.ident = 270,logfc.threshold = 0.3)
c6_findmark <- FindMarkers(c6, ident.1 = "c6_WT", ident.2 = "c6_KO", max.cells.per.ident = 211,logfc.threshold = 0.3)
c7_findmark <- FindMarkers(c7, ident.1 = "c7_WT", ident.2 = "c7_KO",max.cells.per.ident = 187,logfc.threshold = 0.3)
c8_findmark <- FindMarkers(c8, ident.1 = "c8_WT", ident.2 = "c8_KO",max.cells.per.ident = 160,logfc.threshold = 0.3)
c9_findmark <- FindMarkers(c9, ident.1 = "c9_WT", ident.2 = "c9_KO",max.cells.per.ident = 108,logfc.threshold = 0.3)
c10_findmark <- FindMarkers(c10, ident.1 = "c10_WT", ident.2 = "c10_KO",max.cells.per.ident = 93,logfc.threshold = 0.3)
c11_findmark <- FindMarkers(c11, ident.1 = "c11_WT", ident.2 = "c11_KO",max.cells.per.ident = 52,logfc.threshold = 0.3)



# save
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/E105_SR40/DGE/Mesoderm")
write.csv(c0_findmark, file = "c0_findmark.csv")
write.csv(c1_findmark, file = "c1_findmark.csv")
write.csv(c2_findmark, file = "c2_findmark.csv")
write.csv(c3_findmark, file = "c3_findmark.csv")
write.csv(c4_findmark, file = "c4_findmark.csv")
write.csv(c5_findmark, file = "c5_findmark.csv")
write.csv(c6_findmark, file = "c6_findmark.csv")
write.csv(c7_findmark, file = "c7_findmark.csv")
write.csv(c8_findmark, file = "c8_findmark.csv")
write.csv(c9_findmark, file = "c9_findmark.csv")
write.csv(c10_findmark, file = "c10_findmark.csv")
write.csv(c11_findmark, file = "c11_findmark.csv")


## heatmap marker genes
### Endoderm clusters
features <- c("Foxd1","Ebf1",
              "Foxf1","Osr1",
              "Lix1","Alx3",
              "Six1","Eya4",
              "Wt1","Tbx20","Mab21l2","Hand1",
              "Upk3b","Upk1b",
              "Smpx","Myh6")

DoHeatmap(subset(Mesoderm_subset, downsample=100), features = features,slot = "scale.data",assay = "SCT", label = F, raster = F) + scale_fill_gradientn(colors = c("cornflowerblue", "white", "red")) 


new.cluster.ids <- c("PharyngealMesoderm",
                     "pSHF",
                     "PharyngealMesoderm",
                     "AHF",
                     "AHF",
                     "JCF",
                     "Epicardium",
                     "AHF",
                     "PharyngealMesoderm",
                     "pSHF",
                     "Epicardium",
                     "Cardiomyocyte")
names(new.cluster.ids) <- levels(Mesoderm_subset)
Mesoderm_subset_annotated <- RenameIdents(Mesoderm_subset, new.cluster.ids)
DimPlot(Mesoderm_subset_annotated, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()

c0_findmark <- FindMarkers(c0, ident.1 = "c0_WT", ident.2 = "c0_KO",max.cells.per.ident = 563,logfc.threshold = 0.3)
c1_findmark <- FindMarkers(c1, ident.1 = "c1_WT", ident.2 = "c1_KO",max.cells.per.ident = 243,logfc.threshold = 0.3)
c2_findmark <- FindMarkers(c2, ident.1 = "c2_WT", ident.2 = "c2_KO",max.cells.per.ident = 321,logfc.threshold = 0.3)
c3_findmark <- FindMarkers(c3, ident.1 = "c3_WT", ident.2 = "c3_KO",max.cells.per.ident = 278,logfc.threshold = 0.3)
c4_findmark <- FindMarkers(c4, ident.1 = "c4_WT", ident.2 = "c4_KO",max.cells.per.ident = 302,logfc.threshold = 0.3)
c5_findmark <- FindMarkers(c5, ident.1 = "c5_WT", ident.2 = "c5_KO",max.cells.per.ident = 270,logfc.threshold = 0.3)


# Dotplots
my_levels<-c("c4_WT","c4_KO")
c4@active.ident <- factor(x = c4@active.ident, levels = my_levels)
c4_dot <- c("Shox2","Nrxn1","Fgf10","Sema6d","Lhx2","Tcf21","Grid2","Isl1","Ebf1","Foxl2","Cacna2d3","Nav3","Wnt5a")
DotPlot(c4, features = c4_dot) + coord_flip()


my_levels<-c("c7_WT","c7_KO")
c7@active.ident <- factor(x = c7@active.ident, levels = my_levels)
c7_dot <- c("Shox2","Nrxn1","Fgf10","Sema6d","Lhx2","Tcf21","Grid2","Isl1","Ebf1","Foxl2","Cacna2d3","Nav3","Wnt5a")
DotPlot(c7, features = c7_dot) + coord_flip()


my_levels<-c("c3_WT","c3_KO")
c3@active.ident <- factor(x = c3@active.ident, levels = my_levels)
c3_dot <- c("Wnt5a","Sema3c","Fgf5","Rspo3","Sulf1","Chsy3","Nrp1","Adgrv1","Pcdh9")
DotPlot(c3, features = c3_dot) + coord_flip()

my_levels<-c("c2_WT","c2_KO")
c2@active.ident <- factor(x = c2@active.ident, levels = my_levels)
c2_dot <- c("Foxl2","Msx1","Rspo1","Itga4","Tenm2","Eya4","Negr1","Cadm2","Ntrk2","Efnb2")
DotPlot(c2, features = c2_dot) + coord_flip()
