
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1

######## 
# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517")
mesoderm_subset
######## 

## label each time point as a sample and plot umap with 3 time points

# Break up each cluster into timepoint and add as a metadata column as Clusters3
table(mesoderm_subset$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 
# 2741        1777        2065        2742        3622        3717        3799        3321        4507 
# j_E115_WT_1 k_E115_WT_2 l_E115_KO_2 m_E115_KO_3 
# 4472        6678        2919        7387 

condition.vec.01 <- substr(mesoderm_subset$Sample,3,9)

unique(condition.vec.01)
#[1] "E105_WT" "E105_KO" "E115_KO" "E115_WT" "E925_WT" "E925_KO"
table(condition.vec.01)
# condition.vec.01
# E105_KO E105_WT E115_KO E115_WT E925_KO E925_WT 
# 7339    7120   14813   11150    4519    4806

mesoderm_subset$Clusters8 <- as.character(condition.vec.01)
# spot check
table(mesoderm_subset$Clusters8)

p17 <- plotEmbedding(ArchRProj = wt_atlas, colorBy = "cellColData", name = "Clusters3", embedding = "UMAPHarmony")


p17


#########
### plot separate UMAPs based on time point
table(mesoderm_subset$Sample,mesoderm_subset$Clusters6)
## E925_WT
# 1. clone each project 
E925_WT_subset <- mesoderm_subset
table(E925_WT_subset$Clusters8)
# 2. subset
idxSample_E925_WT <- BiocGenerics::which(mesoderm_subset$Clusters8 %in% "E925_WT")
cellsSample_E925_WT <- mesoderm_subset$cellNames[idxSample_E925_WT]
E925_WT_subset<-E925_WT_subset[cellsSample_E925_WT, ]
E925_WT_subset
# numberOfCells(1): 4806
# medianTSS(1): 10.415
# medianFrags(1): 50200.5

# confirm it looks good before saving
plotEmbedding(ArchRProj = E925_WT_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")


## E925_KO
# 1. clone each project 
E925_KO_subset <- mesoderm_subset
# 2. subset
idxSample_E925_KO <- BiocGenerics::which(mesoderm_subset$Clusters8 %in% "E925_KO")
cellsSample_E925_KO <- mesoderm_subset$cellNames[idxSample_E925_KO]
E925_KO_subset<-E925_KO_subset[cellsSample_E925_KO, ]
E925_KO_subset
# numberOfCells(1): 4519
# medianTSS(1): 10.402
# medianFrags(1): 51529

# confirm it looks good before saving
plotEmbedding(ArchRProj = E925_KO_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")


## E105_WT
# 1. clone each project 
E105_WT_subset <- mesoderm_subset
# 2. subset
idxSample_E105_WT <- BiocGenerics::which(mesoderm_subset$Clusters8 %in% "E105_WT")
cellsSample_E105_WT <- mesoderm_subset$cellNames[idxSample_E105_WT]
E105_WT_subset<-E105_WT_subset[cellsSample_E105_WT, ]
E105_WT_subset
# numberOfCells(1): 7120
# medianTSS(1): 13.486
# medianFrags(1): 63804.5

# confirm it looks good before saving
plotEmbedding(ArchRProj = E105_WT_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")


## E105_KO
# 1. clone each project 
E105_KO_subset <- mesoderm_subset
# 2. subset
idxSample_E105_KO <- BiocGenerics::which(mesoderm_subset$Clusters8 %in% "E105_KO")
cellsSample_E105_KO <- mesoderm_subset$cellNames[idxSample_E105_KO]
E105_KO_subset<-E105_KO_subset[cellsSample_E105_KO, ]
E105_KO_subset
# numberOfCells(1): 7339
# medianTSS(1): 13.658
# medianFrags(1): 61630

# confirm it looks good before saving
plotEmbedding(ArchRProj = E105_KO_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")


## E115_WT
# 1. clone each project 
E115_WT_subset <- mesoderm_subset
# 2. subset
idxSample_E115_WT <- BiocGenerics::which(mesoderm_subset$Clusters8 %in% "E115_WT")
cellsSample_E115_WT <- mesoderm_subset$cellNames[idxSample_E115_WT]
E115_WT_subset<-E115_WT_subset[cellsSample_E115_WT, ]
E115_WT_subset
# numberOfCells(1): 11150
# medianTSS(1): 11.9845
# medianFrags(1): 41089.5

# confirm it looks good before saving
plotEmbedding(ArchRProj = E115_WT_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")


## E115_KO
# 1. clone each project 
E115_KO_subset <- mesoderm_subset
# 2. subset
idxSample_E115_KO <- BiocGenerics::which(mesoderm_subset$Clusters8 %in% "E115_KO")
cellsSample_E115_KO <- mesoderm_subset$cellNames[idxSample_E115_KO]
E115_KO_subset<-E115_KO_subset[cellsSample_E115_KO, ]
E115_KO_subset
# numberOfCells(1): 14813
# medianTSS(1): 12.485
# medianFrags(1): 39221

# confirm it looks good before saving
plotEmbedding(ArchRProj = E115_KO_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")




# save
plotEmbedding(ArchRProj = E925_WT_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")
plotEmbedding(ArchRProj = E925_KO_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")
plotEmbedding(ArchRProj = E105_WT_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")
plotEmbedding(ArchRProj = E105_KO_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")
plotEmbedding(ArchRProj = E115_WT_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")
plotEmbedding(ArchRProj = E115_KO_subset, colorBy = "cellColData",name = "Clusters6", embedding = "UMAPHarmony")



