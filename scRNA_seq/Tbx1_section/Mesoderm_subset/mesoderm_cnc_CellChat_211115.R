setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

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

# neural crest subset
Tbx1_CNC <- readRDS("~/scRNA-seq/2020_Tbx1/all_aggr_210429/Tbx1_CNC_annotated.RDS")
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

# 1. subset e11.5 time point
Tbx1_mesoderm_cnc_e115 <- Tbx1_mesoderm_cnc
# wt = 8592 cells
WT_e115 <- OldWhichCells(Tbx1_mesoderm_cnc_e115, subset.name = "gem.group", accept.value = c("WT_E115"))
# ko = 13639 cells
KO_e115 <- OldWhichCells(Tbx1_mesoderm_cnc_e115, subset.name = "gem.group", accept.value = c("KO_E115"))
Tbx1_mesoderm_cnc_e115 <- SubsetData(Tbx1_mesoderm_cnc_e115, cells= c(WT_e115,KO_e115))
Tbx1_mesoderm_cnc_e115
# An object of class Seurat 
# 54328 features across 22231 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(Tbx1_mesoderm_cnc_e115@active.ident, Tbx1_mesoderm_cnc_e115$gem.group)
#                               KO_E115 WT_E115
# Posterior_Second_Heart_Field    2863    1032
# Paraxial_Mesoderm               1549     459
# Cardiopharyngeal_Mesoderm       1561     445
# Cardiopharyngeal_Mesenchyme     2430    1012
# Neural_Crest                    1084    3152
# Cranial_Mesenchyme              1298     872
# Smooth_Muscle                    799     701
# Epicardium                      1147     286
# Anterior_Second_Heart_Field      346     543
# Cardiomyocyte                    562      90

###############
Tbx1_CNC_e115 <- Tbx1_CNC
# wt = 3703 cells
WT_e115 <- OldWhichCells(Tbx1_CNC_e115, subset.name = "gem.group", accept.value = c("WT_E115"))
# ko = 1418 cells
KO_e115 <- OldWhichCells(Tbx1_CNC_e115, subset.name = "gem.group", accept.value = c("KO_E115"))
Tbx1_CNC_e115 <- SubsetData(Tbx1_CNC_e115, cells= c(WT_e115,KO_e115))
Tbx1_CNC_e115
# An object of class Seurat 
# 50565 features across 5121 samples within 2 assays 
# Active assay: SCT (22938 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(Tbx1_CNC_e115@active.ident, Tbx1_CNC_e115$gem.group)
#                             KO_E115 WT_E115
# Craniofacial_NeuralCrest     172    1802
# Cardiac_NeuralCrest          527    1162
# Migratory_NeuralCrest        406     178
# PA3_Cardiac_NeuralCrest      313     561

####
## from scATAC-seq, the following clusters have DARs:
# AHF, Cardiopharyngeal_Mesenchyme, Cardiopharyngeal_Mesoderm, Cranial_Mesenchyme, Paraxial_Mesoderm
# Craniofacial_NeuralCrest, Cardiac_NeuralCrest, Migratory_NeuralCrest, PA3_Cardiac_NeuralCrest
AHF_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Anterior_Second_Heart_Field"),subset.name = "gem.group", accept.value = "WT_E115")
Cardiopharyngeal_Mesenchyme_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Cardiopharyngeal_Mesenchyme"),subset.name = "gem.group", accept.value = "WT_E115")
Cardiopharyngeal_Mesoderm_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Cardiopharyngeal_Mesoderm"),subset.name = "gem.group", accept.value = "WT_E115")
Cranial_Mesenchyme_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Cranial_Mesenchyme"),subset.name = "gem.group", accept.value = "WT_E115")
Paraxial_Mesoderm_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Paraxial_Mesoderm"),subset.name = "gem.group", accept.value = "WT_E115")

AHF_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Anterior_Second_Heart_Field"),subset.name = "gem.group", accept.value = "KO_E115")
Cardiopharyngeal_Mesenchyme_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Cardiopharyngeal_Mesenchyme"),subset.name = "gem.group", accept.value = "KO_E115")
Cardiopharyngeal_Mesoderm_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Cardiopharyngeal_Mesoderm"),subset.name = "gem.group", accept.value = "KO_E115")
Cranial_Mesenchyme_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Cranial_Mesenchyme"),subset.name = "gem.group", accept.value = "KO_E115")
Paraxial_Mesoderm_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e115, ident.keep = c("Paraxial_Mesoderm"),subset.name = "gem.group", accept.value = "KO_E115")

Idents(object = Tbx1_mesoderm_cnc_e115, cells= AHF_wt) <- "wt_AHF"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= Cardiopharyngeal_Mesenchyme_wt) <- "wt_Cardiopharyngeal_Mesenchyme"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= Cardiopharyngeal_Mesoderm_wt) <- "wt_Cardiopharyngeal_Mesoderm"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= Cranial_Mesenchyme_wt) <- "wt_Cranial_Mesenchyme"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= Paraxial_Mesoderm_wt) <- "wt_Paraxial_Mesoderm"

Idents(object = Tbx1_mesoderm_cnc_e115, cells= AHF_ko) <- "wt_AHF_ko"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= Cardiopharyngeal_Mesenchyme_ko) <- "ko_Cardiopharyngeal_Mesenchyme"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= Cardiopharyngeal_Mesoderm_ko) <- "ko_Cardiopharyngeal_Mesoderm"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= Cranial_Mesenchyme_ko) <- "ko_Cranial_Mesenchyme"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= Paraxial_Mesoderm_ko) <- "ko_Paraxial_Mesoderm"
table(Tbx1_mesoderm_cnc_e115@active.ident)

##### add neural crest cells using cell barcodes from subsetted object
wt_Craniofacial_NeuralCrest <- OldWhichCells(Tbx1_CNC_e115, ident.keep = "Craniofacial_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E115")
wt_Cardiac_NeuralCrest  <- OldWhichCells(Tbx1_CNC_e115, ident.keep = "Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E115")
wt_Migratory_NeuralCrest <- OldWhichCells(Tbx1_CNC_e115, ident.keep = "Migratory_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E115")
wt_PA3_Cardiac_NeuralCrest <- OldWhichCells(Tbx1_CNC_e115, ident.keep = "PA3_Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E115")

ko_Craniofacial_NeuralCrest <- OldWhichCells(Tbx1_CNC_e115, ident.keep = "Craniofacial_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E115")
ko_Cardiac_NeuralCrest  <- OldWhichCells(Tbx1_CNC_e115, ident.keep = "Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E115")
ko_Migratory_NeuralCrest <- OldWhichCells(Tbx1_CNC_e115, ident.keep = "Migratory_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E115")
ko_PA3_Cardiac_NeuralCrest <- OldWhichCells(Tbx1_CNC_e115, ident.keep = "PA3_Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E115")

Idents(object = Tbx1_mesoderm_cnc_e115, cells= wt_Craniofacial_NeuralCrest) <- "wt_Craniofacial_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= wt_Cardiac_NeuralCrest) <- "wt_Cardiac_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= wt_Migratory_NeuralCrest) <- "wt_Migratory_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= wt_PA3_Cardiac_NeuralCrest) <- "wt_PA3_Cardiac_NeuralCrest"

Idents(object = Tbx1_mesoderm_cnc_e115, cells= ko_Craniofacial_NeuralCrest) <- "ko_Craniofacial_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= ko_Cardiac_NeuralCrest) <- "ko_Cardiac_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= ko_Migratory_NeuralCrest) <- "ko_Migratory_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e115, cells= ko_PA3_Cardiac_NeuralCrest) <- "ko_PA3_Cardiac_NeuralCrest"

table(Tbx1_mesoderm_cnc_e115@active.ident)

# subset for cellchat
wt_cellchat <- SubsetData(Tbx1_mesoderm_cnc_e115, ident.use = c("wt_Craniofacial_NeuralCrest","wt_Cardiac_NeuralCrest",
                                                                "wt_Migratory_NeuralCrest","wt_PA3_Cardiac_NeuralCrest",
                                                                "wt_AHF","wt_Cardiopharyngeal_Mesenchyme",
                                                                "wt_Cardiopharyngeal_Mesoderm","wt_Cranial_Mesenchyme",
                                                                "wt_Paraxial_Mesoderm"))

wt_cellchat@active.ident <- droplevels(wt_cellchat@active.ident)
table(wt_cellchat@active.ident)
# wt_PA3_Cardiac_NeuralCrest       wt_Migratory_NeuralCrest         wt_Cardiac_NeuralCrest    wt_Craniofacial_NeuralCrest 
# 559                            176                           1158                           1799 
# wt_Paraxial_Mesoderm          wt_Cranial_Mesenchyme   wt_Cardiopharyngeal_Mesoderm wt_Cardiopharyngeal_Mesenchyme 
# 455                            691                            441                           1007 
# wt_AHF 
# 542 










