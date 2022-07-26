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

# 1. subset e10.5 time point
Tbx1_mesoderm_cnc_e105 <- Tbx1_mesoderm_cnc_unannotated_210829
# wt = 5048 cells
WT_e105 <- OldWhichCells(Tbx1_mesoderm_cnc_e105, subset.name = "gem.group", accept.value = c("WT_E105"))
# ko = 6531 cells
KO_e105 <- OldWhichCells(Tbx1_mesoderm_cnc_e105, subset.name = "gem.group", accept.value = c("KO_E105"))
Tbx1_mesoderm_cnc_e105 <- SubsetData(Tbx1_mesoderm_cnc_e105, cells= c(WT_e105,KO_e105))
Tbx1_mesoderm_cnc_e105
# An object of class Seurat 
# 54328 features across 11579 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(Tbx1_mesoderm_cnc_e105@active.ident, Tbx1_mesoderm_cnc_e105$gem.group)
# KO_E105 WT_E105
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
Tbx1_CNC_e105 <- Tbx1_CNC
# wt = 2109 cells
WT_e105 <- OldWhichCells(Tbx1_CNC_e105, subset.name = "gem.group", accept.value = c("WT_E105"))
# ko = 2208 cells
KO_e105 <- OldWhichCells(Tbx1_CNC_e105, subset.name = "gem.group", accept.value = c("KO_E105"))
Tbx1_CNC_e105 <- SubsetData(Tbx1_CNC_e105, cells= c(WT_e105,KO_e105))
Tbx1_CNC_e105
# An object of class Seurat 
# 50565 features across 4317 samples within 2 assays 
# Active assay: SCT (22938 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(Tbx1_CNC_e105@active.ident, Tbx1_CNC_e105$gem.group)
#                             KO_E105 WT_E105
# Craniofacial_NeuralCrest     845     997
# Cardiac_NeuralCrest          935     780
# Migratory_NeuralCrest        304     197
# PA3_Cardiac_NeuralCrest      124     135

####
## for unannotated object, choose clusters: 1,2,3,5,10,11
cluster1_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 1,subset.name = "gem.group", accept.value = "WT_E105")
cluster2_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 2,subset.name = "gem.group", accept.value = "WT_E105")
cluster3_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 3,subset.name = "gem.group", accept.value = "WT_E105")
cluster5_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 5,subset.name = "gem.group", accept.value = "WT_E105")
cluster10_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 10,subset.name = "gem.group", accept.value = "WT_E105")
cluster11_wt <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 11,subset.name = "gem.group", accept.value = "WT_E105")

cluster1_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 1,subset.name = "gem.group", accept.value = "KO_E105")
cluster2_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 2,subset.name = "gem.group", accept.value = "KO_E105")
cluster3_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 3,subset.name = "gem.group", accept.value = "KO_E105")
cluster5_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 5,subset.name = "gem.group", accept.value = "KO_E105")
cluster10_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 10,subset.name = "gem.group", accept.value = "KO_E105")
cluster11_ko <- OldWhichCells(Tbx1_mesoderm_cnc_e105, ident.keep = 11,subset.name = "gem.group", accept.value = "KO_E105")

Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster1_wt) <- "cluster1_wt"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster2_wt) <- "cluster2_wt"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster3_wt) <- "cluster3_wt"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster5_wt) <- "cluster5_wt"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster10_wt) <- "cluster10_wt"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster11_wt) <- "cluster11_wt"

Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster1_ko) <- "cluster1_ko"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster2_ko) <- "cluster2_ko"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster3_ko) <- "cluster3_ko"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster5_ko) <- "cluster5_ko"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster10_ko) <- "cluster10_ko"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= cluster11_ko) <- "cluster11_ko"


##### add neural crest cells using cell barcodes from subsetted object
wt_Craniofacial_NeuralCrest <- OldWhichCells(Tbx1_CNC_e105, ident.keep = "Craniofacial_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E105")
wt_Cardiac_NeuralCrest  <- OldWhichCells(Tbx1_CNC_e105, ident.keep = "Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E105")
wt_Migratory_NeuralCrest <- OldWhichCells(Tbx1_CNC_e105, ident.keep = "Migratory_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E105")
wt_PA3_Cardiac_NeuralCrest <- OldWhichCells(Tbx1_CNC_e105, ident.keep = "PA3_Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "WT_E105")

ko_Craniofacial_NeuralCrest <- OldWhichCells(Tbx1_CNC_e105, ident.keep = "Craniofacial_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E105")
ko_Cardiac_NeuralCrest  <- OldWhichCells(Tbx1_CNC_e105, ident.keep = "Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E105")
ko_Migratory_NeuralCrest <- OldWhichCells(Tbx1_CNC_e105, ident.keep = "Migratory_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E105")
ko_PA3_Cardiac_NeuralCrest <- OldWhichCells(Tbx1_CNC_e105, ident.keep = "PA3_Cardiac_NeuralCrest", subset.name = "gem.group",accept.value = "KO_E105")

Idents(object = Tbx1_mesoderm_cnc_e105, cells= wt_Craniofacial_NeuralCrest) <- "wt_Craniofacial_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= wt_Cardiac_NeuralCrest) <- "wt_Cardiac_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= wt_Migratory_NeuralCrest) <- "wt_Migratory_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= wt_PA3_Cardiac_NeuralCrest) <- "wt_PA3_Cardiac_NeuralCrest"

Idents(object = Tbx1_mesoderm_cnc_e105, cells= ko_Craniofacial_NeuralCrest) <- "ko_Craniofacial_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= ko_Cardiac_NeuralCrest) <- "ko_Cardiac_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= ko_Migratory_NeuralCrest) <- "ko_Migratory_NeuralCrest"
Idents(object = Tbx1_mesoderm_cnc_e105, cells= ko_PA3_Cardiac_NeuralCrest) <- "ko_PA3_Cardiac_NeuralCrest"

## subset out remaining clusters
table(Tbx1_mesoderm_cnc_e105@active.ident) # check to see that tables line up... (a few mismatches but oveall ok)
Tbx1_mesoderm_cnc_e105 <- SubsetData(Tbx1_mesoderm_cnc_e105, ident.remove = c(0,4,6,7,8,9,12,14,15,16))
Tbx1_mesoderm_cnc_e105@active.ident <-droplevels(Tbx1_mesoderm_cnc_e105@active.ident)
table(Tbx1_mesoderm_cnc_e105@active.ident)

# ko_PA3_Cardiac_NeuralCrest    ko_Migratory_NeuralCrest      ko_Cardiac_NeuralCrest ko_Craniofacial_NeuralCrest 
# 123                         298                         927                         831 
# wt_PA3_Cardiac_NeuralCrest    wt_Migratory_NeuralCrest      wt_Cardiac_NeuralCrest wt_Craniofacial_NeuralCrest 
# 131                         196                         772                         986 
# cluster11_ko                cluster10_ko                 cluster5_ko                 cluster3_ko 
# 125                         140                         368                         384 
# cluster2_ko                 cluster1_ko                cluster11_wt                cluster10_wt 
# 534                         722                         239                         106 
# cluster5_wt                 cluster3_wt                 cluster2_wt                 cluster1_wt 
# 259                         291                         409                         561 


# subset for cellchat
wt_cellchat <- SubsetData(Tbx1_mesoderm_cnc_e105, ident.use = c("wt_Migratory_NeuralCrest","wt_Craniofacial_NeuralCrest",
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
# 196                         986                         131                         772 
# cluster1_wt                 cluster2_wt                 cluster3_wt                 cluster5_wt 
# 561                         409                         291                         259 
# cluster10_wt                cluster11_wt 
# 106                         239 

# subset for cellchat
ko_cellchat <- SubsetData(Tbx1_mesoderm_cnc_e105, ident.use = c("ko_Migratory_NeuralCrest","ko_Craniofacial_NeuralCrest",
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
# 298                         831                         123                         927 
# cluster1_ko                 cluster2_ko                 cluster3_ko                 cluster5_ko 
# 722                         534                         384                         368 
# cluster10_ko                cluster11_ko 
# 140                         125









