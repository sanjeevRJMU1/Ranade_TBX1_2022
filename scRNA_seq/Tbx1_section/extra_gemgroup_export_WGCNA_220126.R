
## Script 1: Clustering

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

### 1. Initialize Seurat object
Tbx1_all_aggr.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
Tbx1_all_aggr_sct <- CreateSeuratObject(counts = Tbx1_all_aggr.data, project = "Tbx1_all_aggr", min.cells = 3, min.features = 200)
Tbx1_all_aggr_sct
# An object of class Seurat 
# 27627 features across 104897 samples within 1 assay 
# Active assay: RNA (27627 features, 0 variable features)

### 2. Add mito percentage and gem group to meta data
# mito percentage 
Tbx1_all_aggr_sct[["percent.mt"]] <- PercentageFeatureSet(Tbx1_all_aggr_sct, pattern = "^mt-")

# gem group naming
classification.vec <- as.numeric(gsub(".*-","", (colnames(x = Tbx1_all_aggr_sct))))
names(classification.vec) <- colnames(x = Tbx1_all_aggr_sct)
classification.vec[classification.vec=="1"] <- "WT_E925_1"
classification.vec[classification.vec=="2"] <- "KO_E925_1"
classification.vec[classification.vec=="3"] <- "WT_E925_2"
classification.vec[classification.vec=="4"] <- "KO_E925_2"
classification.vec[classification.vec=="5"] <- "KO_E105_1"
classification.vec[classification.vec=="6"] <- "KO_E105_2"
classification.vec[classification.vec=="7"] <- "WT_E105_1"
classification.vec[classification.vec=="8"] <- "WT_E105_2"
classification.vec[classification.vec=="9"] <- "KO_E115_1"
classification.vec[classification.vec=="10"] <- "WT_E115_1"
classification.vec[classification.vec=="11"] <- "WT_E115_2"
classification.vec[classification.vec=="12"] <- "KO_E115_2"
classification.vec[classification.vec=="13"] <- "KO_E115_3"
Tbx1_all_aggr_sct$"gem.group" <- classification.vec
head(Tbx1_all_aggr_sct@meta.data)
tail(Tbx1_all_aggr_sct@meta.data)
write.csv(Tbx1_all_aggr_sct@meta.data, "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/exports/filename.csv")




