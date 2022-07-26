setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

### 
# re-load
Tbx1_all_aggr_sct <-readRDS(file = "Tbx1_all_aggr_sct.RDS")
DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = T)
Tbx1_all_aggr_sct
# An object of class Seurat 
# 54328 features across 53604 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

# plots
VlnPlot(Tbx1_all_aggr_sct, c("nFeature_RNA"), pt.size = 0.1, ncol = 1)
# all cells between 4k and 9k genes/cell
VlnPlot(Tbx1_all_aggr_sct, c("nCount_RNA"), pt.size = 0.1, ncol = 1)
# all cells < 80k

table(Tbx1_all_aggr_sct@active.ident,Tbx1_all_aggr_sct$gem.group)

#
FeaturePlot(Tbx1_all_aggr_sct, c("Epcam","Sox2","Dlx5","Tnnt2"))
FeaturePlot(Tbx1_all_aggr_sct, c("Dlx2","Sox10","Hand2","Hoxb4"))
FeaturePlot(Tbx1_all_aggr_sct, c("Phox2a","Phox2b","Pecam1",'Hba-x'))
FeaturePlot(Tbx1_all_aggr_sct, c("Tbx1","Meox1","Alx1","Mab21l2"))
FeaturePlot(Tbx1_all_aggr_sct, c("Isl1","Acta2","Tnnt2","Osr1"))
FeaturePlot(Tbx1_all_aggr_sct, c("Fgf8","Fgf10","Wnt5a","Sema3c"))

#
DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = T)

### subset
Tbx1_mesoderm_cnc<-SubsetData(Tbx1_all_aggr_sct, ident.use = c(15,8,10,7,10,13,6,5,3,4,11,16,12,23,9,0,18))
Tbx1_mesoderm_cnc
# An object of class Seurat 
# 54328 features across 38168 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)

table(Tbx1_mesoderm_cnc@active.ident,Tbx1_mesoderm_cnc$gem.group)

##
Tbx1_mesoderm_cnc <- RunUMAP(Tbx1_mesoderm_cnc, dims = 1:20, reduction = "mnn", verbose = TRUE)
Tbx1_mesoderm_cnc <- FindNeighbors(Tbx1_mesoderm_cnc,reduction = "mnn", dims = 1:20, verbose = TRUE, force.recalc = T)
Tbx1_mesoderm_cnc <- FindClusters(Tbx1_mesoderm_cnc, verbose = TRUE, resolution = 0.5)
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = T)
Tbx1_mesoderm_cnc.markers <- FindAllMarkers(Tbx1_mesoderm_cnc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_mesoderm_cnc.markers, file = "Tbx1_mesoderm_cnc-markers.csv")
saveRDS(Tbx1_mesoderm_cnc, file = "Tbx1_mesoderm_cnc.RDS")


###
FeaturePlot(Tbx1_mesoderm_cnc, c("Tbx1","Meox1","Osr1","Fgf8"))
FeaturePlot(Tbx1_mesoderm_cnc, c("Foxd1","Foxc1","Foxc2","Pax9"))
FeaturePlot(Tbx1_mesoderm_cnc, c("Alx1","Wt1","Mab21l2","Twist1"))


########################################################
new.cluster.ids <- c("Mesoderm Progenitors",
                     "Ectoderm",
                     "Endoderm",
                     "Mesoderm Progenitors",
                     "Mesoderm Progenitors",
                     "Mesoderm Progenitors",
                     "Mesoderm Progenitors",
                     "Neural Crest",
                     "Neural Crest",
                     "Mesoderm Progenitors",
                     "Neural Crest",
                     "Mesoderm Progenitors",
                     "Cardiomyocyte",
                     "Neural Crest",
                     "Endothelium",
                     "Neural Crest",
                     "Mesoderm Progenitors",
                     "Ectoderm",
                     "Mesoderm Progenitors",
                     "Endoderm",
                     "Blood",
                     "Endothelium",
                     "Ectoderm",
                     "Mesoderm Progenitors",
                     "Endoderm",
                     "Ectoderm",
                     "Blood",
                     "Endoderm",
                     "Ectoderm")
names(new.cluster.ids) <- levels(Tbx1_all_aggr_sct)
Tbx1_all_aggr_sct <- RenameIdents(Tbx1_all_aggr_sct, new.cluster.ids)
DimPlot(Tbx1_all_aggr_sct, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
saveRDS(Tbx1_mesoderm_cnc, file = "Tbx1_mesoderm_cnc_annotated.RDS")

Tbx1_mesoderm_cnc.markers <- FindAllMarkers(Tbx1_mesoderm_cnc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_mesoderm_cnc.markers, file = "Tbx1_mesoderm_cnc_annotated-markers.csv")

#####
## re-start 08/29
Tbx1_mesoderm_cnc <- readRDS(file = "Tbx1_mesoderm_cnc_annotated.RDS")
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
Tbx1_mesoderm_cnc.markers <- FindAllMarkers(Tbx1_mesoderm_cnc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VlnPlot(Tbx1_mesoderm_cnc, c("nCount_RNA"), pt.size = 0.1, ncol = 1)

Tbx1_mesoderm_cnc <- RunUMAP(Tbx1_mesoderm_cnc, dims = 1:20, reduction = "mnn", verbose = TRUE)
Tbx1_mesoderm_cnc <- FindNeighbors(Tbx1_mesoderm_cnc,reduction = "mnn", dims = 1:20, verbose = TRUE, force.recalc = T)
Tbx1_mesoderm_cnc <- FindClusters(Tbx1_mesoderm_cnc, verbose = TRUE, resolution = 0.5)
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = T)
saveRDS(Tbx1_mesoderm_cnc, file = "Tbx1_mesoderm_cnc_unannotated_210829.RDS")

new.cluster.ids <- c("Posterior_Second_Heart_Field",
                     "Paraxial_Mesoderm",
                     "Cardiopharyngeal_Mesoderm",
                     "Cardiopharyngeal_Mesenchyme",
                     "Neural_Crest",
                     "Cranial_Mesenchyme",
                     "Neural_Crest",
                     "Smooth_Muscle",
                     "Neural_Crest",
                     "Epicardium",
                     "Cardiopharyngeal_Mesenchyme",
                     "Anterior_Second_Heart_Field",
                     "Cardiomyocyte",
                     "Neural_Crest",
                     "Posterior_Second_Heart_Field",
                     "Epicardium",
                     "Epicardium")
names(new.cluster.ids) <- levels(Tbx1_mesoderm_cnc)
Tbx1_mesoderm_cnc <- RenameIdents(Tbx1_mesoderm_cnc, new.cluster.ids)
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()


features <- c("Actc1","Tnnt2","Wt1","Tbx18","Rgs5","Fgf8","Isl1","Dlx2","Dlx5","Lix1","Alx1","Tbx1","Foxd1","Pax9","Meox1","Cxcl13","Osr1")
DotPlot(Tbx1_mesoderm_cnc, features = features, cols = c("grey","red")) + RotatedAxis()
DoHeatmap(subset(Tbx1_mesoderm_cnc, downsample=100), features = features,slot = "scale.data",assay = "SCT", label = F) + scale_fill_gradientn(colors = c("cornflowerblue", "white", "red")) 

