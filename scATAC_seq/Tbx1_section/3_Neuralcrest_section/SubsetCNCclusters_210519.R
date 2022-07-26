###### Tbx1 WT v KO all aggr
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
## here subset out only the CNC cells and then re-analyze downstream


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


#plots
plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
table(mesoderm_subset$Clusters,mesoderm_subset$Sample)

# Subset clusters of interest in the hackiest fucking way 

# assign Y/N values to only the clusters you want
library(plyr)
temp_cluster_names <- as.character(revalue(mesoderm_subset$Clusters, c(
  "C1"="False",
  "C2"="False",
  "C3"="False",
  "C4"="False",
  "C5"="True",
  "C6"="True",
  "C7"="True",
  "C8"="False",
  "C9"="False",
  "C10"="True",
  "C11"="False",
  "C12"="False",
  "C13"="False",
  "C14"="False",
  "C15"="False",
  "C16"="False"
)))
length(unique(temp_cluster_names))
mesoderm_subset$ClustersCNC <- temp_cluster_names
mesoderm_subset@cellColData@listData[["ClustersCNC"]] <- temp_cluster_names
table(mesoderm_subset$ClustersCNC,mesoderm_subset$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 j_E115_WT_1
# False        2419        1599        1989        2440        2512        2188        2277        1895        3991        2676
# True          322         178          76         302        1110        1529        1522        1426         516        1796
# 
# k_E115_WT_2 l_E115_KO_2 m_E115_KO_3
# False        4169        2620        6709
# True         2509         299         678

## subset
idxSample <- BiocGenerics::which(mesoderm_subset$ClustersCNC %in% "True")
cellsSample <- mesoderm_subset$cellNames[idxSample]
subsetArchRProject(
  ArchRProj = mesoderm_subset,
  cells = cellsSample,
  outputDirectory = "/Users/sranade/scATAC-seq/CNC_subset_210519", force = T
)
## Dropping ImputeWeights Since You Are Subsetting Cells! ImputeWeights is a cell-x-cell Matrix!
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/CNC_subset_210519 
# samples(13): m_E115_KO_3 k_E115_WT_2 ... b_E925_KO_1 l_E115_KO_2
# sampleColData names(1): ArrowFiles
# cellColData names(23): Sample TSSEnrichment ... Clusters5 ClustersCNC
# numberOfCells(1): 9091
# medianTSS(1): 11.854
# medianFrags(1): 56857












