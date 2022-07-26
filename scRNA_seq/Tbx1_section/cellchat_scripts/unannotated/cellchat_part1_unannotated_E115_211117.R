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

## load unannotatced mesoderm_cnc subset
Tbx1_mesoderm_cnc_unannotated_210829 <- readRDS("~/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_mesoderm_cnc_unannotated_210829.RDS")
DimPlot(Tbx1_mesoderm_cnc_unannotated_210829, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
table(Tbx1_mesoderm_cnc_unannotated_210829@active.ident, Tbx1_mesoderm_cnc_unannotated_210829$gem.group)

# neural crest subset
Tbx1_CNC <- readRDS("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_CNC_annotated.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

# 1. subset E115 time point
Tbx1_mesoderm_cnc_E115 <- Tbx1_mesoderm_cnc_unannotated_210829
# wt = 5048 cells
WT_E115 <- OldWhichCells(Tbx1_mesoderm_cnc_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
# ko = 6531 cells
KO_E115 <- OldWhichCells(Tbx1_mesoderm_cnc_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Tbx1_mesoderm_cnc_E115 <- SubsetData(Tbx1_mesoderm_cnc_E115, cells= c(WT_E115,KO_E115))
Tbx1_mesoderm_cnc_E115
# An object of class Seurat 
# 54328 features across 11579 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(Tbx1_mesoderm_cnc_E115@active.ident, Tbx1_mesoderm_cnc_E115$gem.group)
# KO_E115 WT_E115
# 0      735     283
# 1      728     563
# 2      536     409
# 3      385     291
# 4      692     661
# 5      408     300
# 6      507     600
# 7      279     192
# 8      525     478
# 9      433     310
# 10     140     106
# 11     126     240
# 12     164      64
# 13     298     198
# 14     149      72
# 15     278     171
# 16     148     110

###############
Tbx1_CNC_E115 <- Tbx1_CNC
# wt = 2109 cells
WT_E115 <- OldWhichCells(Tbx1_CNC_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
# ko = 2208 cells
KO_E115 <- OldWhichCells(Tbx1_CNC_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Tbx1_CNC_E115 <- SubsetData(Tbx1_CNC_E115, cells= c(WT_E115,KO_E115))
Tbx1_CNC_E115
# An object of class Seurat 
# 50565 features across 4317 samples within 2 assays 
# Active assay: SCT (22938 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(Tbx1_CNC_E115@active.ident, Tbx1_CNC_E115$gem.group)
#                             KO_E115 WT_E115
# Craniofacial_NeuralCrest     845     997
# Cardiac_NeuralCrest          935     780
# Migratory_NeuralCrest        304     197
# PA3_Cardiac_NeuralCrest      124     135

####
## for unannotated object, choose clusters: 1,2,3,5,10,11
cluster1_wt <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 1,subset.name = "gem.group", accept.value = "WT_E115")
cluster2_wt <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 2,subset.name = "gem.group", accept.value = "WT_E115")
cluster3_wt <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 3,subset.name = "gem.group", accept.value = "WT_E115")
cluster5_wt <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 5,subset.name = "gem.group", accept.value = "WT_E115")
cluster10_wt <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 10,subset.name = "gem.group", accept.value = "WT_E115")
cluster11_wt <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 11,subset.name = "gem.group", accept.value = "WT_E115")

cluster1_ko <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 1,subset.name = "gem.group", accept.value = "KO_E115")
cluster2_ko <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 2,subset.name = "gem.group", accept.value = "KO_E115")
cluster3_ko <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 3,subset.name = "gem.group", accept.value = "KO_E115")
cluster5_ko <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 5,subset.name = "gem.group", accept.value = "KO_E115")
cluster10_ko <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 10,subset.name = "gem.group", accept.value = "KO_E115")
cluster11_ko <- OldWhichCells(Tbx1_mesoderm_cnc_E115, ident.keep = 11,subset.name = "gem.group", accept.value = "KO_E115")

Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster1_wt) <- "cluster1_wt"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster2_wt) <- "cluster2_wt"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster3_wt) <- "cluster3_wt"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster5_wt) <- "cluster5_wt"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster10_wt) <- "cluster10_wt"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster11_wt) <- "cluster11_wt"

Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster1_ko) <- "cluster1_ko"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster2_ko) <- "cluster2_ko"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster3_ko) <- "cluster3_ko"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster5_ko) <- "cluster5_ko"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster10_ko) <- "cluster10_ko"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= cluster11_ko) <- "cluster11_ko"


##### add neural crest cells using cell barcodes from subsetted object
wt_Craniofacial_NeuralCrest <- OldWhichCells(Tbx1_CNC_E115, ident.keep = "Craniofacial_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E115")
wt_Cardiac_NeuralCrest  <- OldWhichCells(Tbx1_CNC_E115, ident.keep = "Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E115")
wt_Migratory_NeuralCrest <- OldWhichCells(Tbx1_CNC_E115, ident.keep = "Migratory_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E115")
wt_PA3_Cardiac_NeuralCrest <- OldWhichCells(Tbx1_CNC_E115, ident.keep = "PA3_Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E115")

ko_Craniofacial_NeuralCrest <- OldWhichCells(Tbx1_CNC_E115, ident.keep = "Craniofacial_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E115")
ko_Cardiac_NeuralCrest  <- OldWhichCells(Tbx1_CNC_E115, ident.keep = "Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E115")
ko_Migratory_NeuralCrest <- OldWhichCells(Tbx1_CNC_E115, ident.keep = "Migratory_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E115")
ko_PA3_Cardiac_NeuralCrest <- OldWhichCells(Tbx1_CNC_E115, ident.keep = "PA3_Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E115")

Idents(object = Tbx1_mesoderm_cnc_E115, cells= wt_Craniofacial_NeuralCrest) <- "wt_Craniofacial_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= wt_Cardiac_NeuralCrest) <- "wt_Cardiac_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= wt_Migratory_NeuralCrest) <- "wt_Migratory_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= wt_PA3_Cardiac_NeuralCrest) <- "wt_PA3_Cardiac_NeuralCrest"

Idents(object = Tbx1_mesoderm_cnc_E115, cells= ko_Craniofacial_NeuralCrest) <- "ko_Craniofacial_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= ko_Cardiac_NeuralCrest) <- "ko_Cardiac_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= ko_Migratory_NeuralCrest) <- "ko_Migratory_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_E115, cells= ko_PA3_Cardiac_NeuralCrest) <- "ko_PA3_Cardiac_NeuralCrest"

## subset out remaining clusters
table(Tbx1_mesoderm_cnc_E115@active.ident) # check to see that tables line up... (a few mismatches but oveall ok)
Tbx1_mesoderm_cnc_E115 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.remove = c(0,4,6,7,8,9,12,14,15,16))
Tbx1_mesoderm_cnc_E115@active.ident <-droplevels(Tbx1_mesoderm_cnc_E115@active.ident)
table(Tbx1_mesoderm_cnc_E115@active.ident)

# ko_PA3_Cardiac_NeuralCrest    ko_Migratory_NeuralCrest      ko_Cardiac_NeuralCrest ko_Craniofacial_NeuralCrest 
# 312                         405                         527                         172 
# wt_PA3_Cardiac_NeuralCrest    wt_Migratory_NeuralCrest      wt_Cardiac_NeuralCrest wt_Craniofacial_NeuralCrest 
# 559                         176                        1158                        1799 
# cluster11_ko                cluster10_ko                 cluster5_ko                 cluster3_ko 
# 344                         801                        1172                        1624 
# cluster2_ko                 cluster1_ko                cluster11_wt                cluster10_wt 
# 1560                        1546                         542                         401 
# cluster5_wt                 cluster3_wt                 cluster2_wt                 cluster1_wt 
# 691                         606                         441                         455 

# subset for cellchat
wt_cellchat <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = c("wt_Migratory_NeuralCrest","wt_Craniofacial_NeuralCrest",
                                                                "wt_PA3_Cardiac_NeuralCrest","wt_Cardiac_NeuralCrest",
                                                                "cluster1_wt","cluster2_wt","cluster3_wt",
                                                                "cluster5_wt","cluster10_wt","cluster11_wt"))

table(wt_cellchat@active.ident)
wt_cellchat@active.ident <- droplevels(wt_cellchat@active.ident)
table(wt_cellchat@active.ident)

# order idents 
my_levels_wt <- c("wt_Migratory_NeuralCrest","wt_Craniofacial_NeuralCrest",
               "wt_PA3_Cardiac_NeuralCrest","wt_Cardiac_NeuralCrest",
               "cluster1_wt","cluster2_wt","cluster3_wt",
               "cluster5_wt","cluster10_wt","cluster11_wt")

# Re-level object@ident
wt_cellchat@active.ident <- factor(x = wt_cellchat@active.ident, levels = my_levels_wt)
table(wt_cellchat@active.ident)
# wt_Migratory_NeuralCrest wt_Craniofacial_NeuralCrest  wt_PA3_Cardiac_NeuralCrest      wt_Cardiac_NeuralCrest 
# 176                        1799                         559                        1158 
# cluster1_wt                 cluster2_wt                 cluster3_wt                 cluster5_wt 
# 455                         441                         606                         691 
# cluster10_wt                cluster11_wt 
# 401                         542 

# subset for cellchat
ko_cellchat <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = c("ko_Migratory_NeuralCrest","ko_Craniofacial_NeuralCrest",
                                                                "ko_PA3_Cardiac_NeuralCrest","ko_Cardiac_NeuralCrest",
                                                                "cluster1_ko","cluster2_ko","cluster3_ko",
                                                                "cluster5_ko","cluster10_ko","cluster11_ko"))

table(ko_cellchat@active.ident)
ko_cellchat@active.ident <- droplevels(ko_cellchat@active.ident)
table(ko_cellchat@active.ident)

# order idents 
my_levels_ko <- c("ko_Migratory_NeuralCrest","ko_Craniofacial_NeuralCrest",
                  "ko_PA3_Cardiac_NeuralCrest","ko_Cardiac_NeuralCrest",
                  "cluster1_ko","cluster2_ko","cluster3_ko",
                  "cluster5_ko","cluster10_ko","cluster11_ko")

# Re-level object@ident
ko_cellchat@active.ident <- factor(x = ko_cellchat@active.ident, levels = my_levels_ko)
table(ko_cellchat@active.ident)
# ko_Migratory_NeuralCrest ko_Craniofacial_NeuralCrest  ko_PA3_Cardiac_NeuralCrest      ko_Cardiac_NeuralCrest 
# 405                         172                         312                         527 
# cluster1_ko                 cluster2_ko                 cluster3_ko                 cluster5_ko 
# 1546                        1560                        1624                        1172 
# cluster10_ko                cluster11_ko 
# 801                         344 









