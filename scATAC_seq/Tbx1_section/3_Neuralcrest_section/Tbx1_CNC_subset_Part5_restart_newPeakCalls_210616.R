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
table(CNC_subset$Clusters6)
# Cardiac_NeuralCrest Craniofacial_NeuralCrest    Migratory_NeuralCrest  PA3_Cardiac_NeuralCrest 
# 3644                     5964                     1012                     1643

### Peak calling
## 9. Pseudo-bulk Replicates in ArchR
CNC_subset <- addGroupCoverages(ArchRProj = CNC_subset, groupBy = "Clusters6", minCells = 10, maxCells = 5000, force = T, maxFragments = 75 * 10^6)
#2021-06-16 16:48:36 : Finished Creation of Coverage Files!, 14.205 mins elapsed.

## 10. Peak calling
CNC_subset <- addReproduciblePeakSet(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters6", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 5000,
  maxPeaks = 150000,
  minCells = 10,
  force = T
)

# look at the peak set and then save it as a separate object
getPeakSet(CNC_subset)
# GRanges object with 303570 ranges and 13 metadata columns:

merged_peak_setGR <- getPeakSet(CNC_subset)

# add peakmatrix to object
CNC_subset <- addPeakMatrix(CNC_subset)
getAvailableMatrices(CNC_subset)

## Marker Peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters6",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks


heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.5",
  transpose = F,invert = F
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")



plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = CNC_subset, addDOC = FALSE)

#
getGroupBW(ArchRProj = CNC_subset, groupBy = "Clusters6", tileSize = 100, ceiling = 4, normMethod = "None")
getGroupBW(ArchRProj = CNC_subset, groupBy = "Clusters7", tileSize = 100, ceiling = 4, normMethod = "None")


###################################
#
markerTest_PA3_CNC <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PA3_Cardiac_NeuralCrest",
  bgdGroups = "Cardiac_NeuralCrest"
)
pma_PA3_CNC <- plotMarkers(seMarker = markerTest_PA3_CNC, name = "PA3_Cardiac_NeuralCrest", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
pma_PA3_CNC

###################################
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 2", returnGR = TRUE)
markerList


#### hack: export the markerlist genomic ranges object results for a specific cluster to a data frame
Craniofacial_NeuralCrest_df <- as.data.frame(markerList@listData[["Craniofacial_NeuralCrest"]])
Cardiac_NeuralCrest_df <- as.data.frame(markerList@listData[["Cardiac_NeuralCrest"]])
Migratory_NeuralCrest_df <- as.data.frame(markerList@listData[["Migratory_NeuralCrest"]])
PA3_Cardiac_NeuralCrest_df <- as.data.frame(markerList@listData[["PA3_Cardiac_NeuralCrest"]])

## 
# Craniofacial_NeuralCrest = 2734
# Cardiac_NeuralCrest = 5429
# Migratory_NeuralCrest = 3170
# PA3_Cardiac_NeuralCrest = 5050

## format
Craniofacial_NeuralCrest_df <- Craniofacial_NeuralCrest_df[order(Craniofacial_NeuralCrest_df$MeanDiff, decreasing = TRUE),]
Craniofacial_NeuralCrest_df$rank<-seq_len(nrow(Craniofacial_NeuralCrest_df))
Craniofacial_NeuralCrest_df$blank<-""
Craniofacial_NeuralCrest_df<- Craniofacial_NeuralCrest_df[,c(1,2,3,9,10,5)]

Cardiac_NeuralCrest_df <- Cardiac_NeuralCrest_df[order(Cardiac_NeuralCrest_df$MeanDiff, decreasing = TRUE),]
Cardiac_NeuralCrest_df$rank<-seq_len(nrow(Cardiac_NeuralCrest_df))
Cardiac_NeuralCrest_df$blank<-""
Cardiac_NeuralCrest_df<- Cardiac_NeuralCrest_df[,c(1,2,3,9,10,5)]

Migratory_NeuralCrest_df <- Migratory_NeuralCrest_df[order(Migratory_NeuralCrest_df$MeanDiff, decreasing = TRUE),]
Migratory_NeuralCrest_df$rank<-seq_len(nrow(Migratory_NeuralCrest_df))
Migratory_NeuralCrest_df$blank<-""
Migratory_NeuralCrest_df<- Migratory_NeuralCrest_df[,c(1,2,3,9,10,5)]

PA3_Cardiac_NeuralCrest_df <- PA3_Cardiac_NeuralCrest_df[order(PA3_Cardiac_NeuralCrest_df$MeanDiff, decreasing = TRUE),]
PA3_Cardiac_NeuralCrest_df$rank<-seq_len(nrow(PA3_Cardiac_NeuralCrest_df))
PA3_Cardiac_NeuralCrest_df$blank<-""
PA3_Cardiac_NeuralCrest_df<- PA3_Cardiac_NeuralCrest_df[,c(1,2,3,9,10,5)]

### export marker peaks
setwd("/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/bed_files/Cluster_Markers")
write.table(Craniofacial_NeuralCrest_df, file = "Craniofacial_NeuralCrest_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(Cardiac_NeuralCrest_df, file = "Cardiac_NeuralCrest_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(Migratory_NeuralCrest_df, file = "Migratory_NeuralCrest_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(PA3_Cardiac_NeuralCrest_df, file = "PA3_Cardiac_NeuralCrest_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)

###
# Browser Track
markerGenes  <- c("Pitx2","Rgs2","Cux2")

p <- plotBrowserTrack(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters6", 
  geneSymbol = markerGenes, 
  upstream = 100000,
  downstream = 1000
)
grid::grid.newpage()
grid::grid.draw(p$Cux2)
ArchRBrowser(CNC_subset)


##
# 
markerList_PA3vCushion <- getMarkers(markerTest_PA3_CNC, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_PA3vCushion_df <- as.data.frame(markerList_PA3vCushion@listData[["PA3_Cardiac_NeuralCrest"]])
markerList_PA3vCushion_df_homer <- markerList_PA3vCushion_df[order(markerList_PA3vCushion_df$MeanDiff, decreasing = TRUE),]
markerList_PA3vCushion_df_homer$rank <- seq_len(nrow(markerList_PA3vCushion_df_homer))
markerList_PA3vCushion_df_homer$blank <- ""
markerList_PA3vCushion_df_homer <- markerList_PA3vCushion_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_PA3vCushion_df_homer, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/bed_files/Cluster_Markers/markerList_PA3vCushion_df_homer.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


##############


