## re-load object and annotate

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

##
Tbx1_all_aggr_sct <-readRDS(file = "Tbx1_all_aggr_sct.RDS")
DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = T)

Tbx1_all_aggr_sct.markers <- FindAllMarkers(Tbx1_all_aggr_sct, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)

## dotplot/featureplots
FeaturePlot(Tbx1_all_aggr_sct, c("Tnnt2","Epcam","Sox2","Pecam1"))

## cluster annotation
new.cluster.ids <- c("MesodermProgenitor",
                     "Ectoderm",
                     "Endoderm",
                     "MesodermProgenitor",
                     "MesodermProgenitor",
                     "MesodermProgenitor",
                     "MesodermProgenitor",
                     "NeuralCrest",
                     "NeuralCrest",
                     "MesodermProgenitor",
                     "NeuralCrest",
                     "MesodermProgenitor",
                     "Cardiomyocyte",
                     "NeuralCrest",
                     "Endothelium",
                     "NeuralCrest",
                     "MesodermProgenitor",
                     "Ectoderm",
                     "MesodermProgenitor",
                     "Endoderm",
                     "Blood",
                     "Endothelium",
                     "Ectoderm",
                     "MesodermProgenitor",
                     "Endoderm",
                     "Ectoderm",
                     "Blood",
                     "Endoderm",
                     "Ectoderm")
names(new.cluster.ids) <- levels(Tbx1_all_aggr_sct)
Tbx1_all_aggr_sct <- RenameIdents(Tbx1_all_aggr_sct, new.cluster.ids)
DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()
Tbx1_all_aggr_sct.markers <- FindAllMarkers(Tbx1_all_aggr_sct, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
write.csv(Tbx1_all_aggr_sct.markers, file = "Tbx1_all_aggr_sct_all-markers_annotated_210917.csv")
saveRDS(Tbx1_all_aggr_sct, file = "Tbx1_all_aggr_sct_annotated_210917.RDS")


###
markerGenes <- c("Actc1","Myl7","Tnnt2",
                 "Tbx1","Twist1","Osr1",
                 "Dlx5","Dlx1","Dlx2",
                 "Epcam","Cldn7",
                 "Sox2","Hes5",
                 "Plvap","Cdh5",
                 "Hba-x","Hbb-bt")
my_order <-c("Cardiomyocyte",
             "MesodermProgenitor",
             "NeuralCrest",
             "Endoderm",
             "Ectoderm",
             "Endothelium",
             "Blood")
Tbx1_all_aggr_sct@active.ident <- factor(x = Tbx1_all_aggr_sct@active.ident, levels = my_order)

DotPlot(Tbx1_all_aggr_sct, features = markerGenes, cols = c("grey","red"),col.max = 1.5) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev)


#####################
# re-starting on 220103

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

##
Tbx1_all_aggr_sct <-readRDS(file = "Tbx1_all_aggr_sct_annotated_210917.RDS")
DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = T)

Tbx1_all_aggr_sct$active.ident <- Tbx1_all_aggr_sct@active.ident

table(Tbx1_all_aggr_sct$gem.group,Tbx1_all_aggr_sct@active.ident)
