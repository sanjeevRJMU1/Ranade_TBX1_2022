
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

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
##########
Tbx1_CNC <-readRDS(file = "Tbx1_CNC_annotated.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = T)

# cluster re-assignment occurs, which re-assigns clustering in my_levels (assuming you have 12 clusters in total)
my_levels <- c("WT_E925","KO_E925","WT_E105","KO_E105","WT_E115","KO_E115")
Tbx1_CNC@active.ident <- factor(x = Tbx1_CNC@active.ident, levels = my_levels)
DimPlot(Tbx1_CNC, reduction = "umap", label = F, group.by = "gem.group")



#############################
features <- c("Emx2","Smoc1","Ebf1","Dlx1","Dlx2","Dlx3","Dlx5","Dlx6",
              "Hand2","Hand1","Dlk1","Msx2","Foxf1","Gata3","Rgs5","Isl1",
              "Pou3f3","Barx1",
              "Sox2","Sox10","Cdh19",
              "Hoxb3","Hoxd4","Six2","Hoxa3","Tbx2","Meox1","Hoxa4","Mef2c","Hoxb4")

features2 <- c("Hoxa3","Hoxb3","Hoxd4","Six2","Erbb3","Sox10","Hand2","Dlk1","Dlx2","Dlx3","Dlx5")

markerGenes <- c("Erbb3","Sox10","Sox2","Pitx1","Dlx2","Dlx3","Dlx5","Pou3f3","Smoc1","Smoc2","Hand2","Dlk1","Isl1","Foxf1","Rgs5","Hoxa3","Hoxb3","Hoxd4","Hoxa4","Hoxb4","Hoxc4","Six2")
my_order <-c("Migratory_NeuralCrest","Craniofacial_NeuralCrest","Cardiac_NeuralCrest","PA3_Cardiac_NeuralCrest")
Tbx1_CNC@active.ident <- factor(x = Tbx1_CNC@active.ident, levels = my_order)

DotPlot(Tbx1_CNC, features = markerGenes, cols = c("grey","red")) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev)


DoHeatmap(subset(Tbx1_CNC, downsample=100), features = features,slot = "scale.data",assay = "SCT", label = F) + scale_fill_gradientn(colors = c("cornflowerblue", "white", "red")) 
#############################
FeaturePlot(Tbx1_CNC, c("Hoxa3","Hoxb4","Hoxa4","Hoxd4"))
FeaturePlot(Tbx1_CNC, c("Hand1","Rgs5","Isl1","Foxf1"))
FeaturePlot(Tbx1_CNC, c("Emx2","Dlx5","Dlx3","Smoc1"))
FeaturePlot(Tbx1_CNC, c("Sox2","Sox10","Cdh19","Erbb3"))

############
new.cluster.ids <- c("Craniofacial_NeuralCrest",
                     "Cardiac_NeuralCrest",
                     "Craniofacial_NeuralCrest",
                     "Migratory_NeuralCrest",
                     "PA3_Cardiac_NeuralCrest",
                     "Cardiac_NeuralCrest")
names(new.cluster.ids) <- levels(Tbx1_CNC)
Tbx1_CNC <- RenameIdents(Tbx1_CNC, new.cluster.ids)
DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()
saveRDS(Tbx1_CNC, file = "Tbx1_CNC_annotated.RDS")
Tbx1_CNC_annotated.markers <- FindAllMarkers(Tbx1_CNC_annotated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_CNC_annotated.markers, file = "Tbx1_CNC_annotated.markers.csv")
###
Tbx1_CNC <-readRDS(file = "Tbx1_CNC_annotated.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()


###########



## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/E925")
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

Tbx1_CNC_E925 <- Tbx1_CNC

#0
clust0_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 0)
clust0_E925_WT <- OldWhichCells(clust0_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
clust0_E925_KO <- OldWhichCells(clust0_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust0_E925, cells= clust0_E925_WT) <- "WT_E925"
Idents(object = clust0_E925, cells= clust0_E925_KO) <- "KO_E925"
clust0_E925_findmark <- FindMarkers(clust0_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust0_E925_findmark, file="clust0_E925_findmark.csv")

# #1
# clust1_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 1)
# clust1_E925_WT <- OldWhichCells(clust1_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
# clust1_E925_KO <- OldWhichCells(clust1_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
# Idents(object = clust1_E925, cells= clust1_E925_WT) <- "WT_E925"
# Idents(object = clust1_E925, cells= clust1_E925_KO) <- "KO_E925"
# clust1_E925_findmark <- FindMarkers(clust1_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
# write.csv(clust1_E925_findmark,file="clust1_E925_findmark.csv")

#2
clust2_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 2)
clust2_E925_WT <- OldWhichCells(clust2_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
clust2_E925_KO <- OldWhichCells(clust2_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust2_E925, cells= clust2_E925_WT) <- "WT_E925"
Idents(object = clust2_E925, cells= clust2_E925_KO) <- "KO_E925"
clust2_E925_findmark <- FindMarkers(clust2_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust2_E925_findmark,file="clust2_E925_findmark.csv")

#3
clust3_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 3)
clust3_E925_WT <- OldWhichCells(clust3_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
clust3_E925_KO <- OldWhichCells(clust3_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust3_E925, cells= clust3_E925_WT) <- "WT_E925"
Idents(object = clust3_E925, cells= clust3_E925_KO) <- "KO_E925"
clust3_E925_findmark <- FindMarkers(clust3_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust3_E925_findmark,file="clust3_E925_findmark.csv")

# #4
# clust4_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 4)
# clust4_E925_WT <- OldWhichCells(clust4_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
# clust4_E925_KO <- OldWhichCells(clust4_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
# Idents(object = clust4_E925, cells= clust4_E925_WT) <- "WT_E925"
# Idents(object = clust4_E925, cells= clust4_E925_KO) <- "KO_E925"
# clust4_E925_findmark <- FindMarkers(clust4_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
# write.csv(clust4_E925_findmark,file="clust4_E925_findmark.csv")
# 
# #5
# clust5_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 5)
# clust5_E925_WT <- OldWhichCells(clust5_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
# clust5_E925_KO <- OldWhichCells(clust5_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
# Idents(object = clust5_E925, cells= clust5_E925_WT) <- "WT_E925"
# Idents(object = clust5_E925, cells= clust5_E925_KO) <- "KO_E925"
# clust5_E925_findmark <- FindMarkers(clust5_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
# write.csv(clust5_E925_findmark,file="clust5_E925_findmark.csv")



## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/E105")
Tbx1_CNC_E105 <- Tbx1_CNC

#0
clust0_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 0)
clust0_E105_WT <- OldWhichCells(clust0_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust0_E105_KO <- OldWhichCells(clust0_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust0_E105, cells= clust0_E105_WT) <- "WT_E105"
Idents(object = clust0_E105, cells= clust0_E105_KO) <- "KO_E105"
clust0_E105_findmark <- FindMarkers(clust0_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust0_E105_findmark, file="clust0_E105_findmark.csv")

#1
clust1_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 1)
clust1_E105_WT <- OldWhichCells(clust1_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust1_E105_KO <- OldWhichCells(clust1_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust1_E105, cells= clust1_E105_WT) <- "WT_E105"
Idents(object = clust1_E105, cells= clust1_E105_KO) <- "KO_E105"
clust1_E105_findmark <- FindMarkers(clust1_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust1_E105_findmark,file="clust1_E105_findmark.csv")

#2
clust2_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 2)
clust2_E105_WT <- OldWhichCells(clust2_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust2_E105_KO <- OldWhichCells(clust2_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust2_E105, cells= clust2_E105_WT) <- "WT_E105"
Idents(object = clust2_E105, cells= clust2_E105_KO) <- "KO_E105"
clust2_E105_findmark <- FindMarkers(clust2_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust2_E105_findmark,file="clust2_E105_findmark.csv")

#3
clust3_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 3)
clust3_E105_WT <- OldWhichCells(clust3_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust3_E105_KO <- OldWhichCells(clust3_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust3_E105, cells= clust3_E105_WT) <- "WT_E105"
Idents(object = clust3_E105, cells= clust3_E105_KO) <- "KO_E105"
clust3_E105_findmark <- FindMarkers(clust3_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust3_E105_findmark,file="clust3_E105_findmark.csv")

#4
clust4_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 4)
clust4_E105_WT <- OldWhichCells(clust4_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust4_E105_KO <- OldWhichCells(clust4_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust4_E105, cells= clust4_E105_WT) <- "WT_E105"
Idents(object = clust4_E105, cells= clust4_E105_KO) <- "KO_E105"
clust4_E105_findmark <- FindMarkers(clust4_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust4_E105_findmark,file="clust4_E105_findmark.csv")

#5
clust5_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 5)
clust5_E105_WT <- OldWhichCells(clust5_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust5_E105_KO <- OldWhichCells(clust5_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust5_E105, cells= clust5_E105_WT) <- "WT_E105"
Idents(object = clust5_E105, cells= clust5_E105_KO) <- "KO_E105"
clust5_E105_findmark <- FindMarkers(clust5_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust5_E105_findmark,file="clust5_E105_findmark.csv")



## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/E115")
Tbx1_CNC_E115 <- Tbx1_CNC

#0
clust0_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 0)
clust0_E115_WT <- OldWhichCells(clust0_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust0_E115_KO <- OldWhichCells(clust0_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust0_E115, cells= clust0_E115_WT) <- "WT_E115"
Idents(object = clust0_E115, cells= clust0_E115_KO) <- "KO_E115"
clust0_E115_findmark <- FindMarkers(clust0_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust0_E115_findmark, file="clust0_E115_findmark.csv")

#1
clust1_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 1)
clust1_E115_WT <- OldWhichCells(clust1_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust1_E115_KO <- OldWhichCells(clust1_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust1_E115, cells= clust1_E115_WT) <- "WT_E115"
Idents(object = clust1_E115, cells= clust1_E115_KO) <- "KO_E115"
clust1_E115_findmark <- FindMarkers(clust1_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust1_E115_findmark,file="clust1_E115_findmark.csv")

#2
clust2_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 2)
clust2_E115_WT <- OldWhichCells(clust2_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust2_E115_KO <- OldWhichCells(clust2_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust2_E115, cells= clust2_E115_WT) <- "WT_E115"
Idents(object = clust2_E115, cells= clust2_E115_KO) <- "KO_E115"
clust2_E115_findmark <- FindMarkers(clust2_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust2_E115_findmark,file="clust2_E115_findmark.csv")

#3
clust3_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 3)
clust3_E115_WT <- OldWhichCells(clust3_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust3_E115_KO <- OldWhichCells(clust3_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust3_E115, cells= clust3_E115_WT) <- "WT_E115"
Idents(object = clust3_E115, cells= clust3_E115_KO) <- "KO_E115"
clust3_E115_findmark <- FindMarkers(clust3_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust3_E115_findmark,file="clust3_E115_findmark.csv")

#4
clust4_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 4)
clust4_E115_WT <- OldWhichCells(clust4_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust4_E115_KO <- OldWhichCells(clust4_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust4_E115, cells= clust4_E115_WT) <- "WT_E115"
Idents(object = clust4_E115, cells= clust4_E115_KO) <- "KO_E115"
clust4_E115_findmark <- FindMarkers(clust4_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust4_E115_findmark,file="clust4_E115_findmark.csv")

#5
clust5_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 5)
clust5_E115_WT <- OldWhichCells(clust5_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust5_E115_KO <- OldWhichCells(clust5_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust5_E115, cells= clust5_E115_WT) <- "WT_E115"
Idents(object = clust5_E115, cells= clust5_E115_KO) <- "KO_E115"
clust5_E115_findmark <- FindMarkers(clust5_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust5_E115_findmark,file="clust5_E115_findmark.csv")

##########
clust0v2 <- FindMarkers(Tbx1_CNC, ident.1 = 0, ident.2 = 2)

######################