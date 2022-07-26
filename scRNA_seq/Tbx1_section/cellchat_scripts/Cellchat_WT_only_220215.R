###### Cellchat of WT cells only to define which signaling pathways are most relevant for the 4 neural crest populations
#### here take all WT cells from mesoderm/cnc analysis


library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)
library(CellChat)
library(patchwork)

## load mesoderm_cnc subset
Tbx1_mesoderm_cnc <- readRDS(file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_mesoderm_cnc_annotated_210830.RDS")
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
table(Tbx1_mesoderm_cnc@active.ident, Tbx1_mesoderm_cnc$gem.group)

#                                 KO_E105 KO_E115 KO_E925 WT_E105 WT_E115 WT_E925
# Posterior_Second_Heart_Field     884    2863     299     355    1032     175
# Paraxial_Mesoderm                728    1549     284     563     459     300
# Cardiopharyngeal_Mesoderm        536    1561     169     409     445     149
# Cardiopharyngeal_Mesenchyme      525    2430     236     397    1012     185
# Neural_Crest                    2022    1084     482    1937    3152     262
# Cranial_Mesenchyme               408    1298      83     300     872      56
# Smooth_Muscle                    279     799     101     192     701      87
# Epicardium                       859    1147     406     591     286     265
# Anterior_Second_Heart_Field      126     346     146     240     543     109
# Cardiomyocyte                    164     562     294      64      90     270

# neural crest subset
Tbx1_CNC <- readRDS("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_CNC_annotated.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)
#                           KO_E105 KO_E115 KO_E925 WT_E105 WT_E115 WT_E925
# Craniofacial_NeuralCrest     845     172     279     997    1802     131
# Cardiac_NeuralCrest          935     527       1     780    1162       0
# Migratory_NeuralCrest        304     406     226     197     178     149
# PA3_Cardiac_NeuralCrest      124     313       1     135     561       6


# 1. subset wt cells
Tbx1_mesoderm_cnc_wt_only <- Tbx1_mesoderm_cnc
# wt = 15,498 cells
WT_e115 <- OldWhichCells(Tbx1_mesoderm_cnc_wt_only, subset.name = "gem.group", accept.value = c("WT_E115","WT_E925","WT_E105"))
Tbx1_mesoderm_cnc_wt_only <- SubsetData(Tbx1_mesoderm_cnc_wt_only, cells= c(WT_e115))
Tbx1_mesoderm_cnc_wt_only
# An object of class Seurat 
# 54328 features across 15498 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

###############
Tbx1_CNC_wt_only <- Tbx1_CNC
# wt = 6098 cells
WT_e115 <- OldWhichCells(Tbx1_CNC_wt_only, subset.name = "gem.group", accept.value = c("WT_E115","WT_E105","WT_E925"))
Tbx1_CNC_wt_only <- SubsetData(Tbx1_CNC_wt_only, cells= c(WT_e115))
Tbx1_CNC_wt_only
# An object of class Seurat 
# 50565 features across 6098 samples within 2 assays 
# Active assay: SCT (22938 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

##### add neural crest cells using cell barcodes from subsetted object
wt_Craniofacial_NeuralCrest <- OldWhichCells(Tbx1_CNC_wt_only, ident.keep = "Craniofacial_NeuralCrest", subset.name = "gem.group",accept.value = c("WT_E115","WT_E105","WT_E925"))
# 2930
wt_Cardiac_NeuralCrest  <- OldWhichCells(Tbx1_CNC_wt_only, ident.keep = "Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = c("WT_E115","WT_E105"))
# 1942
wt_Migratory_NeuralCrest <- OldWhichCells(Tbx1_CNC_wt_only, ident.keep = "Migratory_NeuralCrest", subset.name = "gem.group",accept.value = c("WT_E115","WT_E105","WT_E925"))
#524
wt_PA3_Cardiac_NeuralCrest <- OldWhichCells(Tbx1_CNC_wt_only, ident.keep = "PA3_Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = c("WT_E115","WT_E105","WT_E925"))
#702

Idents(object = Tbx1_mesoderm_cnc_wt_only, cells= wt_Craniofacial_NeuralCrest) <- "wt_Craniofacial_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_wt_only, cells= wt_Cardiac_NeuralCrest) <- "wt_Cardiac_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_wt_only, cells= wt_Migratory_NeuralCrest) <- "wt_Migratory_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_wt_only, cells= wt_PA3_Cardiac_NeuralCrest) <- "wt_PA3_Cardiac_NeuralCrest"

Tbx1_mesoderm_cnc_wt_only <- SubsetData(Tbx1_mesoderm_cnc_wt_only,ident.remove = "Neural_Crest")
Tbx1_mesoderm_cnc_wt_only@active.ident<-droplevels(Tbx1_mesoderm_cnc_wt_only@active.ident)
table(Tbx1_mesoderm_cnc_wt_only@active.ident)
# wt_PA3_Cardiac_NeuralCrest     wt_Migratory_NeuralCrest       wt_Cardiac_NeuralCrest 
# 695                          504                         1930 
# wt_Craniofacial_NeuralCrest Posterior_Second_Heart_Field            Paraxial_Mesoderm 
# 2908                         1558                         1316 
# Cardiopharyngeal_Mesoderm  Cardiopharyngeal_Mesenchyme           Cranial_Mesenchyme 
# 999                         1589                         1006 
# Smooth_Muscle                   Epicardium  Anterior_Second_Heart_Field 
# 488                         1142                          890 
# Cardiomyocyte 
# 424 


# order idents 
my_levels_wt <- c("wt_Migratory_NeuralCrest","wt_Craniofacial_NeuralCrest",
                  "wt_PA3_Cardiac_NeuralCrest","wt_Cardiac_NeuralCrest",
                  "Anterior_Second_Heart_Field","Cardiopharyngeal_Mesoderm","Cardiopharyngeal_Mesenchyme",
                  "Cranial_Mesenchyme","Smooth_Muscle","Epicardium","Paraxial_Mesoderm",
                  "Posterior_Second_Heart_Field","Cardiomyocyte")

# Re-level object@ident
Tbx1_mesoderm_cnc_wt_only@active.ident <- factor(x = Tbx1_mesoderm_cnc_wt_only@active.ident, levels = my_levels_wt)
table(Tbx1_mesoderm_cnc_wt_only@active.ident)
# wt_Migratory_NeuralCrest  wt_Craniofacial_NeuralCrest   wt_PA3_Cardiac_NeuralCrest 
# 504                         2908                          695 
# wt_Cardiac_NeuralCrest  Anterior_Second_Heart_Field    Cardiopharyngeal_Mesoderm 
# 1930                          890                          999 
# Cardiopharyngeal_Mesenchyme           Cranial_Mesenchyme                Smooth_Muscle 
# 1589                         1006                          488 
# Epicardium            Paraxial_Mesoderm Posterior_Second_Heart_Field 
# 1142                         1316                         1558 
# Cardiomyocyte 
# 424



