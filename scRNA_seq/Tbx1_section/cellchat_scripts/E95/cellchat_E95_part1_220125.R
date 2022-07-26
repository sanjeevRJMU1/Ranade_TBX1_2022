### Part 1
### filter cells and create 2 cellchat objects: 1 wt, 1 ko - both with all 4 NC clusters + mesoderm clusters with DARs

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

# # neural crest subset
# Tbx1_CNC <- readRDS("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_CNC_annotated.RDS")
# DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
# table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

# 1. subset E9.5 time point
Tbx1_mesoderm_cnc_e95 <- Tbx1_mesoderm_cnc
# wt = 8592 cells
WT_E925 <- OldWhichCells(Tbx1_mesoderm_cnc_e95, subset.name = "gem.group", accept.value = c("WT_E925"))
# ko = 13639 cells
KO_E925 <- OldWhichCells(Tbx1_mesoderm_cnc_e95, subset.name = "gem.group", accept.value = c("KO_E925"))
Tbx1_mesoderm_cnc_e95 <- SubsetData(Tbx1_mesoderm_cnc_e95, cells= c(WT_E925,KO_E925))
Tbx1_mesoderm_cnc_e95
# An object of class Seurat 
# 54328 features across 4358 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(Tbx1_mesoderm_cnc_e95@active.ident, Tbx1_mesoderm_cnc_e95$gem.group)
#                                 KO_E925 WT_E925
# Posterior_Second_Heart_Field     299     175
# Paraxial_Mesoderm                284     300
# Cardiopharyngeal_Mesoderm        169     149
# Cardiopharyngeal_Mesenchyme      236     185
# Neural_Crest                     482     262
# Cranial_Mesenchyme                83      56
# Smooth_Muscle                    101      87
# Epicardium                       406     265
# Anterior_Second_Heart_Field      146     109
# Cardiomyocyte                    294     270


####
## from scATAC-seq, the following clusters have DARs:
# AHF, Cardiopharyngeal_Mesenchyme, Cardiopharyngeal_Mesoderm, Cranial_Mesenchyme, Paraxial_Mesoderm
# Craniofacial_NeuralCrest, Cardiac_NeuralCrest, Migratory_NeuralCrest, PA3_Cardiac_NeuralCrest
pSHF_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Posterior_Second_Heart_Field"),subset.name = "gem.group", accept.value = "WT_E925")
Paraxial_Mesoderm_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Paraxial_Mesoderm"),subset.name = "gem.group", accept.value = "WT_E925")
Cardiopharyngeal_Mesoderm_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Cardiopharyngeal_Mesoderm"),subset.name = "gem.group", accept.value = "WT_E925")
Cardiopharyngeal_Mesenchyme_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Cardiopharyngeal_Mesenchyme"),subset.name = "gem.group", accept.value = "WT_E925")
Neural_Crest_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Neural_Crest"),subset.name = "gem.group", accept.value = "WT_E925")
Cranial_Mesenchyme_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Cranial_Mesenchyme"),subset.name = "gem.group", accept.value = "WT_E925")
Smooth_Muscle_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Smooth_Muscle"),subset.name = "gem.group", accept.value = "WT_E925")
Epicardium_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Epicardium"),subset.name = "gem.group", accept.value = "WT_E925")
Anterior_Second_Heart_Field_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Anterior_Second_Heart_Field"),subset.name = "gem.group", accept.value = "WT_E925")
Cardiomyocyte_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Cardiomyocyte"),subset.name = "gem.group", accept.value = "WT_E925")


pSHF_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Posterior_Second_Heart_Field"),subset.name = "gem.group", accept.value = "KO_E925")
Paraxial_Mesoderm_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Paraxial_Mesoderm"),subset.name = "gem.group", accept.value = "KO_E925")
Cardiopharyngeal_Mesoderm_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Cardiopharyngeal_Mesoderm"),subset.name = "gem.group", accept.value = "KO_E925")
Cardiopharyngeal_Mesenchyme_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Cardiopharyngeal_Mesenchyme"),subset.name = "gem.group", accept.value = "KO_E925")
Neural_Crest_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Neural_Crest"),subset.name = "gem.group", accept.value = "KO_E925")
Cranial_Mesenchyme_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Cranial_Mesenchyme"),subset.name = "gem.group", accept.value = "KO_E925")
Smooth_Muscle_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Smooth_Muscle"),subset.name = "gem.group", accept.value = "KO_E925")
Epicardium_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Epicardium"),subset.name = "gem.group", accept.value = "KO_E925")
Anterior_Second_Heart_Field_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Anterior_Second_Heart_Field"),subset.name = "gem.group", accept.value = "KO_E925")
Cardiomyocyte_KO <- OldWhichCells(Tbx1_mesoderm_cnc_e95, ident.keep = c("Cardiomyocyte"),subset.name = "gem.group", accept.value = "KO_E925")

##
Idents(object = Tbx1_mesoderm_cnc_e95, cells= pSHF_wt) <- "pSHF_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Paraxial_Mesoderm_wt) <- "Paraxial_Mesoderm_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Cardiopharyngeal_Mesoderm_wt) <- "Cardiopharyngeal_Mesoderm_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Cardiopharyngeal_Mesenchyme_wt) <- "Cardiopharyngeal_Mesenchyme_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Neural_Crest_wt) <- "Neural_Crest_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Cranial_Mesenchyme_wt) <- "Cranial_Mesenchyme_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Smooth_Muscle_wt) <- "Smooth_Muscle_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Epicardium_wt) <- "Epicardium_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Anterior_Second_Heart_Field_wt) <- "Anterior_Second_Heart_Field_wt"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Cardiomyocyte_wt) <- "Cardiomyocyte_wt"


Idents(object = Tbx1_mesoderm_cnc_e95, cells= pSHF_KO) <- "pSHF_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Paraxial_Mesoderm_KO) <- "Paraxial_Mesoderm_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Cardiopharyngeal_Mesoderm_KO) <- "Cardiopharyngeal_Mesoderm_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Cardiopharyngeal_Mesenchyme_KO) <- "Cardiopharyngeal_Mesenchyme_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Neural_Crest_KO) <- "Neural_Crest_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Cranial_Mesenchyme_KO) <- "Cranial_Mesenchyme_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Smooth_Muscle_KO) <- "Smooth_Muscle_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Epicardium_KO) <- "Epicardium_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Anterior_Second_Heart_Field_KO) <- "Anterior_Second_Heart_Field_KO"
Idents(object = Tbx1_mesoderm_cnc_e95, cells= Cardiomyocyte_KO) <- "Cardiomyocyte_KO"


table(Tbx1_mesoderm_cnc_e95@active.ident)

# subset for cellchat
wt_cellchat <- SubsetData(Tbx1_mesoderm_cnc_e95, ident.use = c("Anterior_Second_Heart_Field_wt",
                                                               "Cardiopharyngeal_Mesoderm_wt",
                                                               "Neural_Crest_wt",
                                                               "Cardiopharyngeal_Mesenchyme_wt",
                                                               "Paraxial_Mesoderm_wt",
                                                               "Cranial_Mesenchyme_wt",
                                                               "Cardiomyocyte_wt",
                                                               "Epicardium_wt",
                                                               "Smooth_Muscle_wt",
                                                               "pSHF_wt"))

wt_cellchat@active.ident <- droplevels(wt_cellchat@active.ident)
table(wt_cellchat@active.ident)
# Cardiomyocyte_wt Anterior_Second_Heart_Field_wt                  Epicardium_wt 
# 270                            109                            265 
# Smooth_Muscle_wt          Cranial_Mesenchyme_wt                Neural_Crest_wt 
# 87                             56                            262 
# Cardiopharyngeal_Mesenchyme_wt   Cardiopharyngeal_Mesoderm_wt           Paraxial_Mesoderm_wt 
# 185                            149                            300 
# pSHF_wt 
# 175
# order idents 
my_levels_wt <- c("Anterior_Second_Heart_Field_wt",
                  "Cardiopharyngeal_Mesoderm_wt",
                  "Neural_Crest_wt",
                  "Cardiopharyngeal_Mesenchyme_wt",
                  "Paraxial_Mesoderm_wt",
                  "Cranial_Mesenchyme_wt",
                  "Cardiomyocyte_wt",
                  "Epicardium_wt",
                  "Smooth_Muscle_wt",
                  "pSHF_wt")

# Re-level object@ident
wt_cellchat@active.ident <- factor(x = wt_cellchat@active.ident, levels = my_levels_wt)
table(wt_cellchat@active.ident)
# wt_Migratory_NeuralCrest    wt_Craniofacial_NeuralCrest     wt_PA3_Cardiac_NeuralCrest         wt_Cardiac_NeuralCrest 
# 176                           1799                            559                           1158 
# wt_AHF wt_Cardiopharyngeal_Mesenchyme   wt_Cardiopharyngeal_Mesoderm          wt_Cranial_Mesenchyme 
# 542                           1007                            441                            691 
# wt_Paraxial_Mesoderm 
# 455


# subset for cellchat
KO_cellchat <- SubsetData(Tbx1_mesoderm_cnc_e95, ident.use = c("Anterior_Second_Heart_Field_KO",
                                                               "Cardiopharyngeal_Mesoderm_KO",
                                                               "Neural_Crest_KO",
                                                               "Cardiopharyngeal_Mesenchyme_KO",
                                                               "Paraxial_Mesoderm_KO",
                                                               "Cranial_Mesenchyme_KO",
                                                               "Cardiomyocyte_KO",
                                                               "Epicardium_KO",
                                                               "Smooth_Muscle_KO",
                                                               "pSHF_KO"))

KO_cellchat@active.ident <- droplevels(KO_cellchat@active.ident)
table(KO_cellchat@active.ident)
# Cardiomyocyte_KO Anterior_Second_Heart_Field_KO                  Epicardium_KO 
# 270                            109                            265 
# Smooth_Muscle_KO          Cranial_Mesenchyme_KO                Neural_Crest_KO 
# 87                             56                            262 
# Cardiopharyngeal_Mesenchyme_KO   Cardiopharyngeal_Mesoderm_KO           Paraxial_Mesoderm_KO 
# 185                            149                            300 
# pSHF_KO 
# 175
# order idents 
my_levels_KO <- c("Anterior_Second_Heart_Field_KO",
                  "Cardiopharyngeal_Mesoderm_KO",
                  "Neural_Crest_KO",
                  "Cardiopharyngeal_Mesenchyme_KO",
                  "Paraxial_Mesoderm_KO",
                  "Cranial_Mesenchyme_KO",
                  "Cardiomyocyte_KO",
                  "Epicardium_KO",
                  "Smooth_Muscle_KO",
                  "pSHF_KO")

# Re-level object@ident
KO_cellchat@active.ident <- factor(x = KO_cellchat@active.ident, levels = my_levels_KO)
table(KO_cellchat@active.ident)
# Anterior_Second_Heart_Field_KO   Cardiopharyngeal_Mesoderm_KO                Neural_Crest_KO 
# 146                            169                            482 
# Cardiopharyngeal_Mesenchyme_KO           Paraxial_Mesoderm_KO          Cranial_Mesenchyme_KO 
# 236                            284                             83 
# Cardiomyocyte_KO                  Epicardium_KO               Smooth_Muscle_KO 
# 294                            406                            101 
# pSHF_KO 
# 299 







