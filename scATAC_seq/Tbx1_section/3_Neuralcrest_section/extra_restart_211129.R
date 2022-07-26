## scripts on 11/29

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
table(CNC_subset$Clusters6)
# Cardiac_NeuralCrest Craniofacial_NeuralCrest    Migratory_NeuralCrest  PA3_Cardiac_NeuralCrest 
# 3644                     5964                     1012                     1643


##### Dar plots, save 5x5 pdf
markerTest_Cardiac_NeuralCrest_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiac_NeuralCrest_05_W",
  bgdGroups = "Cardiac_NeuralCrest_05_K"
)
pma_Cardiac_NeuralCrest_E105 <- plotMarkers(seMarker = markerTest_Cardiac_NeuralCrest_E105, name = "Cardiac_NeuralCrest_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cardiac_NeuralCrest_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiac_NeuralCrest_15_W",
  bgdGroups = "Cardiac_NeuralCrest_15_K"
)
pma_Cardiac_NeuralCrest_E115 <- plotMarkers(seMarker = markerTest_Cardiac_NeuralCrest_E115, name = "Cardiac_NeuralCrest_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


markerTest_PA3_Cardiac_NeuralCrest_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PA3_Cardiac_NeuralCrest_05_W",
  bgdGroups = "PA3_Cardiac_NeuralCrest_05_K"
)
pma_PA3_Cardiac_NeuralCrest_E105 <- plotMarkers(seMarker = markerTest_PA3_Cardiac_NeuralCrest_E105, name = "PA3_Cardiac_NeuralCrest_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")



pma_Cardiac_NeuralCrest_E105
pma_Cardiac_NeuralCrest_E115
pma_PA3_Cardiac_NeuralCrest_E105

####### make UMAP of WT v KO per time point (so there should be 6 total)
# Break up each cluster into timepoint and add as a metadata column as Clusters9
table(CNC_subset$Sample)

condition.vec.01 <- substr(CNC_subset$Sample,3,9)

unique(condition.vec.01)
#[1] "E105_WT" "E105_KO" "E115_KO" "E115_WT" "E925_WT" "E925_KO"

length(unique(condition.vec.01))
# 6

table(condition.vec.01)
# condition.vec.01
# E105_KO E105_WT E115_KO E115_WT E925_KO E925_WT 
# 2639    2948    1493    4305     480     398 

CNC_subset$Clusters9 <- as.character(condition.vec.01)
# spot check
table(CNC_subset$Clusters9)

p <- plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Clusters9", embedding = "UMAPHarmony")
p
plotPDF(p,
        name = "Plot-UMAP-Manual-Annotations-by-TimexGenotype.pdf",
        ArchRProj = CNC_subset,
        addDOC = FALSE, width = 5, height = 5)







##### for last figure
markerTest_Migratory_NeuralCrest_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Migratory_NeuralCrest_25_W",
  bgdGroups = "Migratory_NeuralCrest_25_K"
)
pma_Migratory_NeuralCrest_E925 <- plotMarkers(seMarker = markerTest_Migratory_NeuralCrest_E925, name = "Migratory_NeuralCrest_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Migratory_NeuralCrest_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Migratory_NeuralCrest_05_W",
  bgdGroups = "Migratory_NeuralCrest_05_K"
)
pma_Migratory_NeuralCrest_E105 <- plotMarkers(seMarker = markerTest_Migratory_NeuralCrest_E105, name = "Migratory_NeuralCrest_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Migratory_NeuralCrest_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Migratory_NeuralCrest_15_W",
  bgdGroups = "Migratory_NeuralCrest_15_K"
)
pma_Migratory_NeuralCrest_E115 <- plotMarkers(seMarker = markerTest_Migratory_NeuralCrest_E115, name = "Migratory_NeuralCrest_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


pma_Migratory_NeuralCrest_E925
pma_Migratory_NeuralCrest_E105
pma_Migratory_NeuralCrest_E115
