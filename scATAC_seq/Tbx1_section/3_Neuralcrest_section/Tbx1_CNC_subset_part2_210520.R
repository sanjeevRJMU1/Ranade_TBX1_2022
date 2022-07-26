## 9. Pseudo-bulk Replicates in ArchR
CNC_subset <- addGroupCoverages(ArchRProj = CNC_subset, groupBy = "Clusters6", minCells = 40, maxCells = 1500, force = T)
#020-09-28 19:35:16 : Finished Creation of Coverage Files!, 39.599 mins elapsed.

## 10. Peak calling
CNC_subset <- addReproduciblePeakSet(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters6", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 2000,
  maxPeaks = 150000,
  minCells = 25,
  force = T
)

#2021-05-19 14:19:55 : Finished Creating Union Peak Set (306228)!, 49.041 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(CNC_subset)
merged_peak_setGR <- getPeakSet(CNC_subset)

# add peakmatrix to object
CNC_subset <- addPeakMatrix(CNC_subset)
getAvailableMatrices(CNC_subset)

saveArchRProject(ArchRProj = CNC_subset, outputDirectory = "/Users/sranade/scATAC-seq/CNC_subset_210519", load = TRUE)
#####################
###############################################
# II. Script part 1, running DARs on each cluster in Clusters6

table(CNC_subset$Clusters6)

### II. Script part 2, running DARs on each cluster in manually annotated Clusters3
## For all clusters, use FDR <= 0.05 ad log2FC > 1. This is prob the safe for FDR and still allows you to catch as many 
## peaks as possible since its not entirely clear what "1.8x more accessible actually means"

#
markerTest_C1_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1_25_W",
  bgdGroups = "C1_25_K"
)
pma_C1_E925 <- plotMarkers(seMarker = markerTest_C1_E925, name = "C1_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C1_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1_05_W",
  bgdGroups = "C1_05_K"
)
pma_C1_E105 <- plotMarkers(seMarker = markerTest_C1_E105, name = "C1_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C1_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1_15_W",
  bgdGroups = "C1_15_K"
)
pma_C1_E115 <- plotMarkers(seMarker = markerTest_C1_E115, name = "C1_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C2_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2_25_W",
  bgdGroups = "C2_25_K"
)
pma_C2_E925 <- plotMarkers(seMarker = markerTest_C2_E925, name = "C2_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C2_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2_05_W",
  bgdGroups = "C2_05_K"
)
pma_C2_E105 <- plotMarkers(seMarker = markerTest_C2_E105, name = "C2_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C2_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2_15_W",
  bgdGroups = "C2_15_K"
)
pma_C2_E115 <- plotMarkers(seMarker = markerTest_C2_E115, name = "C2_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C3_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C3_25_W",
  bgdGroups = "C3_25_K"
)
pma_C3_E925 <- plotMarkers(seMarker = markerTest_C3_E925, name = "C3_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C3_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C3_05_W",
  bgdGroups = "C3_05_K"
)
pma_C3_E105 <- plotMarkers(seMarker = markerTest_C3_E105, name = "C3_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C3_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C3_15_W",
  bgdGroups = "C3_15_K"
)
pma_C3_E115 <- plotMarkers(seMarker = markerTest_C3_E115, name = "C3_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C4_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4_25_W",
  bgdGroups = "C4_25_K"
)
pma_C4_E925 <- plotMarkers(seMarker = markerTest_C4_E925, name = "C4_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C4_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4_05_W",
  bgdGroups = "C4_05_K"
)
pma_C4_E105 <- plotMarkers(seMarker = markerTest_C4_E105, name = "C4_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C4_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4_15_W",
  bgdGroups = "C4_15_K"
)
pma_C4_E115 <- plotMarkers(seMarker = markerTest_C4_E115, name = "C4_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C5_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_25_W",
  bgdGroups = "C5_25_K"
)
pma_C5_E925 <- plotMarkers(seMarker = markerTest_C5_E925, name = "C5_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C5_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_05_W",
  bgdGroups = "C5_05_K"
)
pma_C5_E105 <- plotMarkers(seMarker = markerTest_C5_E105, name = "C5_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C5_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_15_W",
  bgdGroups = "C5_15_K"
)
pma_C5_E115 <- plotMarkers(seMarker = markerTest_C5_E115, name = "C5_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C6_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_25_W",
  bgdGroups = "C6_25_K"
)
pma_C6_E925 <- plotMarkers(seMarker = markerTest_C6_E925, name = "C6_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C6_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_05_W",
  bgdGroups = "C6_05_K"
)
pma_C6_E105 <- plotMarkers(seMarker = markerTest_C6_E105, name = "C6_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C6_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_15_W",
  bgdGroups = "C6_15_K"
)
pma_C6_E115 <- plotMarkers(seMarker = markerTest_C6_E115, name = "C6_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C7_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C7_25_W",
  bgdGroups = "C7_25_K"
)
pma_C7_E925 <- plotMarkers(seMarker = markerTest_C7_E925, name = "C7_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C7_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C7_05_W",
  bgdGroups = "C7_05_K"
)
pma_C7_E105 <- plotMarkers(seMarker = markerTest_C7_E105, name = "C7_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_C7_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters6",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C7_15_W",
  bgdGroups = "C7_15_K"
)
pma_C7_E115 <- plotMarkers(seMarker = markerTest_C7_E115, name = "C7_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
##########################################################################################





##########################################################################################
###############
## 11. Marker Peaks (can use either PeakMatrix or GeneScoreMatrix!)
markersPeaks <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
# class: SummarizedExperiment 
# dim: 289804 7 
# metadata(2): MatchInfo Params
# assays(7): Log2FC Mean ... AUC MeanBGD
# rownames(289804): 1 2 ... 289803 289804
# rowData names(4): seqnames idx start end
# colnames(7): C1 C2 ... C6 C7
# colData names(0):


markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = CNC_subset, addDOC = FALSE)

#### hack: export the markerlist genomic ranges object results for a specific cluster to a data frame
c1_df <- as.data.frame(markerList@listData[["C1"]])
c2_df <- as.data.frame(markerList@listData[["C2"]])
c3_df <- as.data.frame(markerList@listData[["C3"]])
c4_df <- as.data.frame(markerList@listData[["C4"]])
c5_df <- as.data.frame(markerList@listData[["C5"]])
c6_df <- as.data.frame(markerList@listData[["C6"]])
c7_df <- as.data.frame(markerList@listData[["C7"]])
# c8_df <- as.data.frame(markerList@listData[["C8"]])
# c9_df <- as.data.frame(markerList@listData[["C9"]])
# c10_df <- as.data.frame(markerList@listData[["C10"]])
# c11_df <- as.data.frame(markerList@listData[["C11"]])
# c12_df <- as.data.frame(markerList@listData[["C12"]])
# c13_df <- as.data.frame(markerList@listData[["C13"]])
# c14_df <- as.data.frame(markerList@listData[["C14"]])
# c15_df <- as.data.frame(markerList@listData[["C15"]])
# c16_df <- as.data.frame(markerList@listData[["C16"]])
# c17_df <- as.data.frame(markerList@listData[["C17"]])
# c18_df <- as.data.frame(markerList@listData[["C18"]])
# c19_df <- as.data.frame(markerList@listData[["C19"]])

## format
c1_df <- c1_df[order(c1_df$Log2FC, decreasing = TRUE),]
c1_df$rank<-seq_len(nrow(c1_df))
c1_df$blank<-""
c1_df<- c1_df[,c(1,2,3,9,10,5)]

c2_df <- c2_df[order(c2_df$Log2FC, decreasing = TRUE),]
c2_df$rank<-seq_len(nrow(c2_df))
c2_df$blank<-""
c2_df<- c2_df[,c(1,2,3,9,10,5)]

c3_df <- c3_df[order(c3_df$Log2FC, decreasing = TRUE),]
c3_df$rank<-seq_len(nrow(c3_df))
c3_df$blank<-""
c3_df<- c3_df[,c(1,2,3,9,10,5)]

c4_df <- c4_df[order(c4_df$Log2FC, decreasing = TRUE),]
c4_df$rank<-seq_len(nrow(c4_df))
c4_df$blank<-""
c4_df<- c4_df[,c(1,2,3,9,10,5)]

c5_df <- c5_df[order(c5_df$Log2FC, decreasing = TRUE),]
c5_df$rank<-seq_len(nrow(c5_df))
c5_df$blank<-""
c5_df<- c5_df[,c(1,2,3,9,10,5)]

c6_df <- c6_df[order(c6_df$Log2FC, decreasing = TRUE),]
c6_df$rank<-seq_len(nrow(c6_df))
c6_df$blank<-""
c6_df<- c6_df[,c(1,2,3,9,10,5)]

c7_df <- c7_df[order(c7_df$Log2FC, decreasing = TRUE),]
c7_df$rank<-seq_len(nrow(c7_df))
c7_df$blank<-""
c7_df<- c7_df[,c(1,2,3,9,10,5)]

# c8_df <- c8_df[order(c8_df$Log2FC, decreasing = TRUE),]
# c8_df$rank<-seq_len(nrow(c8_df))
# c8_df$blank<-""
# c8_df<- c8_df[,c(1,2,3,9,10,5)]
# 
# c9_df <- c9_df[order(c9_df$Log2FC, decreasing = TRUE),]
# c9_df$rank<-seq_len(nrow(c9_df))
# c9_df$blank<-""
# c9_df<- c9_df[,c(1,2,3,9,10,5)]
# 
# c10_df <- c10_df[order(c10_df$Log2FC, decreasing = TRUE),]
# c10_df$rank<-seq_len(nrow(c10_df))
# c10_df$blank<-""
# c10_df<- c10_df[,c(1,2,3,9,10,5)]
# 
# c11_df <- c11_df[order(c11_df$Log2FC, decreasing = TRUE),]
# c11_df$rank<-seq_len(nrow(c11_df))
# c11_df$blank<-""
# c11_df<- c11_df[,c(1,2,3,9,10,5)]
# 
# c12_df <- c12_df[order(c12_df$Log2FC, decreasing = TRUE),]
# c12_df$rank<-seq_len(nrow(c12_df))
# c12_df$blank<-""
# c12_df<- c12_df[,c(1,2,3,9,10,5)]
# 
# c13_df <- c13_df[order(c13_df$Log2FC, decreasing = TRUE),]
# c13_df$rank<-seq_len(nrow(c13_df))
# c13_df$blank<-""
# c13_df<- c13_df[,c(1,2,3,9,10,5)]
# 
# c14_df <- c14_df[order(c14_df$Log2FC, decreasing = TRUE),]
# c14_df$rank<-seq_len(nrow(c14_df))
# c14_df$blank<-""
# c14_df<- c14_df[,c(1,2,3,9,10,5)]
# 
# c15_df <- c15_df[order(c15_df$Log2FC, decreasing = TRUE),]
# c15_df$rank<-seq_len(nrow(c15_df))
# c15_df$blank<-""
# c15_df<- c15_df[,c(1,2,3,9,10,5)]
# 
# c16_df <- c16_df[order(c16_df$Log2FC, decreasing = TRUE),]
# c16_df$rank<-seq_len(nrow(c16_df))
# c16_df$blank<-""
# c16_df<- c16_df[,c(1,2,3,9,10,5)]
# 
# c17_df <- c17_df[order(c17_df$Log2FC, decreasing = TRUE),]
# c17_df$rank<-seq_len(nrow(c17_df))
# c17_df$blank<-""
# c17_df<- c17_df[,c(1,2,3,9,10,5)]
# 
# c18_df <- c18_df[order(c18_df$Log2FC, decreasing = TRUE),]
# c18_df$rank<-seq_len(nrow(c18_df))
# c18_df$blank<-""
# c18_df<- c18_df[,c(1,2,3,9,10,5)]
# 
# c19_df <- c19_df[order(c19_df$Log2FC, decreasing = TRUE),]
# c19_df$rank<-seq_len(nrow(c19_df))
# c19_df$blank<-""
# c19_df<- c19_df[,c(1,2,3,9,10,5)]

### export marker peaks
setwd("/Users/sranade/scATAC-seq/CNC_subset_210519/Exports/CNC_Clusters_markerpeaks")
write.table(c1_df, file = "c1_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c2_df, file = "c2_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c3_df, file = "c3_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c4_df, file = "c4_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c5_df, file = "c5_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c6_df, file = "c6_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c7_df, file = "c7_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c8_df, file = "c8_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c9_df, file = "c9_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c10_df, file = "c10_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c11_df, file = "c11_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c12_df, file = "c12_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c13_df, file = "c13_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c14_df, file = "c14_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c15_df, file = "c15_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c16_df, file = "c16_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c17_df, file = "c17_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c18_df, file = "c18_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
# write.table(c19_df, file = "c19_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)

### what is going on with cluster 2? compare to cluster 1
markerTest_clust2v1 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2",
  bgdGroups = "C1"
)
pma_C2 <- plotMarkers(seMarker = markerTest_clust2v1, name = "C2", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#######################################
#### 7.3 Identifying Marker Genes
markersGS <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  binarize = F,
  maxCells = 1000
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")

#### Visualize as heatmap
markerGenes <- c("Tnnt2","Dlx5","Wt1","Isl1","Foxf1","Fgf8","Tbx1","Mab21l2","Meox1","Osr1","Fgf10","Wnt5a","Lix1")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
library(ComplexHeatmap)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = CNC_subset, addDOC = FALSE)

## spot check GAS
plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

CNC_subset <- addImputeWeights(CNC_subset, reducedDims = "Harmony")
markerGenes_spotmarkerGenes_spot <-c("Foxc1","Foxc2","Foxf1","Foxa1","Foxl2","Foxd2")
p_spotcheck <- plotEmbedding(
  ArchRProj = CNC_subset,
  colorBy = "GeneScoreMatrix",
  name = markerGenes_spotmarkerGenes_spot,
  continuousSet = "horizonExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(CNC_subset)
)
p_spotcheck$Foxd2

