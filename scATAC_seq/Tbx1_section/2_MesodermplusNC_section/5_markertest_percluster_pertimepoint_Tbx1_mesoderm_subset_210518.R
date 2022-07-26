###### Tbx1 WT v KO all aggr
###### E11.5 only, SR43
### Chapter 11: Differential Accessibility Testing
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


######## 
# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517")
mesoderm_subset
# numberOfCells(1): 49747
# medianTSS(1): 12.221
# medianFrags(1): 47380
######## 



###############################################
# II. Script part 1, running DARs on each cluster in Clusters5

table(mesoderm_subset$Clusters5)

### II. Script part 2, running DARs on each cluster in manually annotated Clusters3
## For all clusters, use FDR <= 0.05 ad log2FC > 1. This is prob the safe for FDR and still allows you to catch as many 
## peaks as possible since its not entirely clear what "1.8x more accessible actually means"

#
markerTest_C1_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1_25_W",
  bgdGroups = "C1_25_K"
)
pma_C1_E925 <- plotMarkers(seMarker = markerTest_C1_E925, name = "C1_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C1_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1_05_W",
  bgdGroups = "C1_05_K"
)
pma_C1_E105 <- plotMarkers(seMarker = markerTest_C1_E105, name = "C1_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C1_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1_15_W",
  bgdGroups = "C1_15_K"
)
pma_C1_E115 <- plotMarkers(seMarker = markerTest_C1_E115, name = "C1_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C2_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2_25_W",
  bgdGroups = "C2_25_K"
)
pma_C2_E925 <- plotMarkers(seMarker = markerTest_C2_E925, name = "C2_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C2_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2_05_W",
  bgdGroups = "C2_05_K"
)
pma_C2_E105 <- plotMarkers(seMarker = markerTest_C2_E105, name = "C2_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C2_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2_15_W",
  bgdGroups = "C2_15_K"
)
pma_C2_E115 <- plotMarkers(seMarker = markerTest_C2_E115, name = "C2_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C3_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C3_25_W",
  bgdGroups = "C3_25_K"
)
pma_C3_E925 <- plotMarkers(seMarker = markerTest_C3_E925, name = "C3_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C3_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C3_05_W",
  bgdGroups = "C3_05_K"
)
pma_C3_E105 <- plotMarkers(seMarker = markerTest_C3_E105, name = "C3_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C3_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C3_15_W",
  bgdGroups = "C3_15_K"
)
pma_C3_E115 <- plotMarkers(seMarker = markerTest_C3_E115, name = "C3_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C4_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4_25_W",
  bgdGroups = "C4_25_K"
)
pma_C4_E925 <- plotMarkers(seMarker = markerTest_C4_E925, name = "C4_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C4_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4_05_W",
  bgdGroups = "C4_05_K"
)
pma_C4_E105 <- plotMarkers(seMarker = markerTest_C4_E105, name = "C4_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C4_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4_15_W",
  bgdGroups = "C4_15_K"
)
pma_C4_E115 <- plotMarkers(seMarker = markerTest_C4_E115, name = "C4_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C5_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_25_W",
  bgdGroups = "C5_25_K"
)
pma_C5_E925 <- plotMarkers(seMarker = markerTest_C5_E925, name = "C5_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C5_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_05_W",
  bgdGroups = "C5_05_K"
)
pma_C5_E105 <- plotMarkers(seMarker = markerTest_C5_E105, name = "C5_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C5_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_15_W",
  bgdGroups = "C5_15_K"
)
pma_C5_E115 <- plotMarkers(seMarker = markerTest_C5_E115, name = "C5_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C6_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_25_W",
  bgdGroups = "C6_25_K"
)
pma_C6_E925 <- plotMarkers(seMarker = markerTest_C6_E925, name = "C6_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C6_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_05_W",
  bgdGroups = "C6_05_K"
)
pma_C6_E105 <- plotMarkers(seMarker = markerTest_C6_E105, name = "C6_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C6_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_15_W",
  bgdGroups = "C6_15_K"
)
pma_C6_E115 <- plotMarkers(seMarker = markerTest_C6_E115, name = "C6_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C7_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C7_25_W",
  bgdGroups = "C7_25_K"
)
pma_C7_E925 <- plotMarkers(seMarker = markerTest_C7_E925, name = "C7_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C7_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C7_05_W",
  bgdGroups = "C7_05_K"
)
pma_C7_E105 <- plotMarkers(seMarker = markerTest_C7_E105, name = "C7_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C7_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C7_15_W",
  bgdGroups = "C7_15_K"
)
pma_C7_E115 <- plotMarkers(seMarker = markerTest_C7_E115, name = "C7_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C8_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C8_25_W",
  bgdGroups = "C8_25_K"
)
pma_C8_E925 <- plotMarkers(seMarker = markerTest_C8_E925, name = "C8_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C8_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C8_05_W",
  bgdGroups = "C8_05_K"
)
pma_C8_E105 <- plotMarkers(seMarker = markerTest_C8_E105, name = "C8_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C8_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C8_15_W",
  bgdGroups = "C8_15_K"
)
pma_C8_E115 <- plotMarkers(seMarker = markerTest_C8_E115, name = "C8_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C9_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C9_25_W",
  bgdGroups = "C9_25_K"
)
pma_C9_E925 <- plotMarkers(seMarker = markerTest_C9_E925, name = "C9_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C9_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C9_05_W",
  bgdGroups = "C9_05_K"
)
pma_C9_E105 <- plotMarkers(seMarker = markerTest_C9_E105, name = "C9_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C9_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C9_15_W",
  bgdGroups = "C9_15_K"
)
pma_C9_E115 <- plotMarkers(seMarker = markerTest_C9_E115, name = "C9_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C10_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C10_25_W",
  bgdGroups = "C10_25_K"
)
pma_C10_E925 <- plotMarkers(seMarker = markerTest_C10_E925, name = "C10_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C10_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C10_05_W",
  bgdGroups = "C10_05_K"
)
pma_C10_E105 <- plotMarkers(seMarker = markerTest_C10_E105, name = "C10_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C10_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C10_15_W",
  bgdGroups = "C10_15_K"
)
pma_C10_E115 <- plotMarkers(seMarker = markerTest_C10_E115, name = "C10_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C11_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C11_25_W",
  bgdGroups = "C11_25_K"
)
pma_C11_E925 <- plotMarkers(seMarker = markerTest_C11_E925, name = "C11_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C11_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C11_05_W",
  bgdGroups = "C11_05_K"
)
pma_C11_E105 <- plotMarkers(seMarker = markerTest_C11_E105, name = "C11_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C11_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C11_15_W",
  bgdGroups = "C11_15_K"
)
pma_C11_E115 <- plotMarkers(seMarker = markerTest_C11_E115, name = "C11_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C12_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C12_25_W",
  bgdGroups = "C12_25_K"
)
pma_C12_E925 <- plotMarkers(seMarker = markerTest_C12_E925, name = "C12_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C12_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C12_05_W",
  bgdGroups = "C12_05_K"
)
pma_C12_E105 <- plotMarkers(seMarker = markerTest_C12_E105, name = "C12_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C12_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C12_15_W",
  bgdGroups = "C12_15_K"
)
pma_C12_E115 <- plotMarkers(seMarker = markerTest_C12_E115, name = "C12_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C13_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C13_25_W",
  bgdGroups = "C13_25_K"
)
pma_C13_E925 <- plotMarkers(seMarker = markerTest_C13_E925, name = "C13_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C13_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C13_05_W",
  bgdGroups = "C13_05_K"
)
pma_C13_E105 <- plotMarkers(seMarker = markerTest_C13_E105, name = "C13_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C13_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C13_15_W",
  bgdGroups = "C13_15_K"
)
pma_C13_E115 <- plotMarkers(seMarker = markerTest_C13_E115, name = "C13_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C14_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C14_25_W",
  bgdGroups = "C14_25_K"
)
pma_C14_E925 <- plotMarkers(seMarker = markerTest_C14_E925, name = "C14_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C14_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C14_05_W",
  bgdGroups = "C14_05_K"
)
pma_C14_E105 <- plotMarkers(seMarker = markerTest_C14_E105, name = "C14_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C14_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C14_15_W",
  bgdGroups = "C14_15_K"
)
pma_C14_E115 <- plotMarkers(seMarker = markerTest_C14_E115, name = "C14_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C15_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C15_25_W",
  bgdGroups = "C15_25_K"
)
pma_C15_E925 <- plotMarkers(seMarker = markerTest_C15_E925, name = "C15_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C15_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C15_05_W",
  bgdGroups = "C15_05_K"
)
pma_C15_E105 <- plotMarkers(seMarker = markerTest_C15_E105, name = "C15_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C15_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C15_15_W",
  bgdGroups = "C15_15_K"
)
pma_C15_E115 <- plotMarkers(seMarker = markerTest_C15_E115, name = "C15_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C16_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C16_25_W",
  bgdGroups = "C16_25_K"
)
pma_C16_E925 <- plotMarkers(seMarker = markerTest_C16_E925, name = "C16_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C16_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C16_05_W",
  bgdGroups = "C16_05_K"
)
pma_C16_E105 <- plotMarkers(seMarker = markerTest_C16_E105, name = "C16_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C16_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters5",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C16_15_W",
  bgdGroups = "C16_15_K"
)
pma_C16_E115 <- plotMarkers(seMarker = markerTest_C16_E115, name = "C16_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

































# export PDF into illustrator






