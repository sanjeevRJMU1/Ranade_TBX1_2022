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
Tbx1_E105_sct <- readRDS("Tbx1_E105_sct.RDS")
DimPlot(Tbx1_E105_sct, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()



# Endoderm Subset
Endoderm_subset <- SubsetData(Tbx1_E105_sct, ident.use = c(30,27,22,7,12,29))
Endoderm_subset <- RunUMAP(Endoderm_subset, dims = 1:30, reduction = "mnn", verbose = TRUE)
Endoderm_subset <- FindNeighbors(Endoderm_subset,reduction = "mnn", dims = 1:30, verbose = TRUE, force.recalc = T)
Endoderm_subset <- FindClusters(Endoderm_subset, verbose = TRUE, resolution = 0.8)
#p3 <- DimPlot(Endoderm_subset, reduction = "umap", label = T, pt.size = 0.2)
#p4 <- DimPlot(Endoderm_subset, reduction = "umap", group.by = "gem.group", pt.size = 0.2)
#CombinePlots(plots = list(p3, p4))
DimPlot(Endoderm_subset, reduction = "umap", label = T, pt.size = 0.2)
VlnPlot(Endoderm_subset, c("nFeature_RNA","nCount_RNA"), pt.size = 0.1, ncol = 1)

Endoderm_subset <- RunUMAP(Endoderm_subset, dims = 1:30, reduction = "mnn", verbose = TRUE)
Endoderm_subset <- FindNeighbors(Endoderm_subset,reduction = "mnn", dims = 1:30, verbose = TRUE, force.recalc = T)
Endoderm_subset <- FindClusters(Endoderm_subset, verbose = TRUE, resolution = 0.5)
DimPlot(Endoderm_subset, reduction = "umap", label = T, pt.size = 0.3)
table(Endoderm_subset@active.ident, Endoderm_subset$gem.group)

Endoderm_subset.markers <- FindAllMarkers(Endoderm_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Endoderm_subset.markers, file = "Endoderm_subset.markers.csv")
saveRDS(Endoderm_subset, file = "Endoderm_subset.RDS")


### DGE

### DGE
# subset each cluster
table(Endoderm_subset@active.ident, Endoderm_subset$gem.group)
c0 <- SubsetData(Endoderm_subset, ident.use = 0)
c1 <- SubsetData(Endoderm_subset, ident.use = 1)
c2 <- SubsetData(Endoderm_subset, ident.use = 2)
c3 <- SubsetData(Endoderm_subset, ident.use = 3)
c4 <- SubsetData(Endoderm_subset, ident.use = 4)
c5 <- SubsetData(Endoderm_subset, ident.use = 5)
c6 <- SubsetData(Endoderm_subset, ident.use = 6)
c7 <- SubsetData(Endoderm_subset, ident.use = 7)
c8 <- SubsetData(Endoderm_subset, ident.use = 8)

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


# 3. Run Findmarkers -- normally you set the lower cell # but since its n=1, just run it
c0_findmark <- FindMarkers(c0, ident.1 = "c0_WT", ident.2 = "c0_KO",max.cells.per.ident = 150,logfc.threshold = 0.3)
c1_findmark <- FindMarkers(c1, ident.1 = "c1_WT", ident.2 = "c1_KO",max.cells.per.ident = 135,logfc.threshold = 0.3)
c2_findmark <- FindMarkers(c2, ident.1 = "c2_WT", ident.2 = "c2_KO",max.cells.per.ident = 108,logfc.threshold = 0.3)
c3_findmark <- FindMarkers(c3, ident.1 = "c3_WT", ident.2 = "c3_KO",max.cells.per.ident = 127,logfc.threshold = 0.3)
c4_findmark <- FindMarkers(c4, ident.1 = "c4_WT", ident.2 = "c4_KO",max.cells.per.ident = 126,logfc.threshold = 0.3)
c5_findmark <- FindMarkers(c5, ident.1 = "c5_WT", ident.2 = "c5_KO",max.cells.per.ident = 79,logfc.threshold = 0.3)
c6_findmark <- FindMarkers(c6, ident.1 = "c6_WT", ident.2 = "c6_KO", max.cells.per.ident = 69,logfc.threshold = 0.3)
c7_findmark <- FindMarkers(c7, ident.1 = "c7_WT", ident.2 = "c7_KO",max.cells.per.ident = 43,logfc.threshold = 0.3)
c8_findmark <- FindMarkers(c8, ident.1 = "c8_WT", ident.2 = "c8_KO",max.cells.per.ident = 32,logfc.threshold = 0.3)



# save
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/E105_SR40/DGE/Endoderm")
write.csv(c0_findmark, file = "c0_findmark.csv")
write.csv(c1_findmark, file = "c1_findmark.csv")
write.csv(c2_findmark, file = "c2_findmark.csv")
write.csv(c3_findmark, file = "c3_findmark.csv")
write.csv(c4_findmark, file = "c4_findmark.csv")
write.csv(c5_findmark, file = "c5_findmark.csv")
write.csv(c6_findmark, file = "c6_findmark.csv")
write.csv(c7_findmark, file = "c7_findmark.csv")
write.csv(c8_findmark, file = "c8_findmark.csv")


## Exporting figures
Endoderm_subset <- readRDS("~/scRNA-seq/2020_Tbx1/E105_SR40/Endoderm_subset.RDS")
DimPlot(Endoderm_subset, reduction = "umap", label = T, pt.size = 0.3)
DimPlot(Endoderm_subset, reduction = "umap", label = F, pt.size = 0.3, group.by = "gem.group")

table(Endoderm_subset@active.ident, Endoderm_subset$gem.group)

## heatmap marker genes
features <- c("Nkx2-3","Pitx2",
              "Foxi2","Wnt6",
              "Rfx6","Foxa2",
              "Hapln1","Wnt3",
              "Pax9","Pax1",
              "Foxe1","Sox10",
              "Foxo2","Hand1","Rgs5",
              "Fgb","Fbb")

DoHeatmap(subset(Endoderm_subset, downsample=100), features = features,slot = "scale.data",assay = "SCT", label = F, raster = F) + scale_fill_gradientn(colors = c("cornflowerblue", "white", "red")) 


## Dotplots
# Dotplots
my_levels<-c("c5_WT","c5_KO")
c5@active.ident <- factor(x = c5@active.ident, levels = my_levels)
c5_dot <- c("Pdlim3","Hoxb4","Edn1","Foxl2","Cxcl12","Bmp6","Fgf8","Pax9","Ripply3","Nkx2-6")
DotPlot(c5, features = c5_dot) + coord_flip()











