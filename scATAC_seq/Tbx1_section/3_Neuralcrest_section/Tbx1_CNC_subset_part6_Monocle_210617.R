###### Tbx1 WT v KO CNC subset
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


######## 
# Re-load ArchR project!
CNC_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/CNC_subset_210607")
CNC_subset
######## 
# numberOfCells(1): 12263
# medianTSS(1): 12.054
# medianFrags(1): 54893
########
table(CNC_subset$Clusters7)

# assign Y/N values to only the clusters you want
library(plyr)
temp_cluster_names <- as.character(revalue(CNC_subset$Clusters7, c(
  "Cardiac_NeuralCrest_05_K"="False",
  "Cardiac_NeuralCrest_05_W"="False",
  "Cardiac_NeuralCrest_15_K"="False",
  "Cardiac_NeuralCrest_15_W"="False",
  "Cardiac_NeuralCrest_25_K"="False",
  "Cardiac_NeuralCrest_25_W"="False",
  "Craniofacial_NeuralCrest_05_K"="False",
  "Craniofacial_NeuralCrest_05_W"="False",
  "Craniofacial_NeuralCrest_15_K"="False",
  "Craniofacial_NeuralCrest_15_W"="False",
  "Migratory_NeuralCrest_05_K"="False",
  "Migratory_NeuralCrest_05_W"="False",
  "Migratory_NeuralCrest_15_K"="False",
  "Migratory_NeuralCrest_15_W"="False",
  "Migratory_NeuralCrest_25_K"="False",
  "Migratory_NeuralCrest_25_W"="False",
  "PA3_Cardiac_NeuralCrest_05_K"="True",
  "PA3_Cardiac_NeuralCrest_05_W"="True",
  "PA3_Cardiac_NeuralCrest_15_K"="True",
  "PA3_Cardiac_NeuralCrest_15_W"="True",
  "PA3_Cardiac_NeuralCrest_25_K"="False"
)))
length(unique(temp_cluster_names))
CNC_subset$MonocleSubset <- temp_cluster_names
CNC_subset@cellColData@listData[["MonocleSubset"]] <- temp_cluster_names
table(CNC_subset$MonocleSubset,CNC_subset$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 j_E115_WT_1
# False        2419        1599        1989        2440        2512        2188        2277        1895        3991        2676
# True          322         178          76         302        1110        1529        1522        1426         516        1796
# 
# k_E115_WT_2 l_E115_KO_2 m_E115_KO_3
# False        4169        2620        6709
# True         2509         299         678

## subset
idxSample <- BiocGenerics::which(CNC_subset$MonocleSubset %in% "True")
cellsSample <- CNC_subset$cellNames[idxSample]
subsetArchRProject(
  ArchRProj = CNC_subset,
  cells = cellsSample,
  outputDirectory = "/Users/sranade/scATAC-seq/CNC_subset_Monocle_210617", force = T
)
