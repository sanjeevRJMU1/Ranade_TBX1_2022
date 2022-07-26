###### Tbx1 WT v KO mesoderm subset1
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
#### Chapters 11 and 12: Marker Peaks and Motif Enrichment

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset1_210517")
mesoderm_subset
######## 

## 11. Marker Peaks (can use either PeakMatrix or GeneScoreMatrix!)
markersPeaks <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
# class: SummarizedExperiment 
# dim: 413016 19 
# metadata(2): MatchInfo Params
# assays(7): Log2FC Mean ... AUC MeanBGD
# rownames(413016): 1 2 ... 413015 413016
# rowData names(4): seqnames idx start end
# colnames(19): C1 C2 ... C18 C19
# colData names(0):


markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = mesoderm_subset, addDOC = FALSE)

#### hack: export the markerlist genomic ranges object results for a specific cluster to a data frame
c1_df <- as.data.frame(markerList@listData[["C1"]])
c2_df <- as.data.frame(markerList@listData[["C2"]])
c3_df <- as.data.frame(markerList@listData[["C3"]])
c4_df <- as.data.frame(markerList@listData[["C4"]])
c5_df <- as.data.frame(markerList@listData[["C5"]])
c6_df <- as.data.frame(markerList@listData[["C6"]])
c7_df <- as.data.frame(markerList@listData[["C7"]])
c8_df <- as.data.frame(markerList@listData[["C8"]])
c9_df <- as.data.frame(markerList@listData[["C9"]])
c10_df <- as.data.frame(markerList@listData[["C10"]])
c11_df <- as.data.frame(markerList@listData[["C11"]])
c12_df <- as.data.frame(markerList@listData[["C12"]])
c13_df <- as.data.frame(markerList@listData[["C13"]])
c14_df <- as.data.frame(markerList@listData[["C14"]])
c15_df <- as.data.frame(markerList@listData[["C15"]])
c16_df <- as.data.frame(markerList@listData[["C16"]])
c17_df <- as.data.frame(markerList@listData[["C17"]])
c18_df <- as.data.frame(markerList@listData[["C18"]])
c19_df <- as.data.frame(markerList@listData[["C19"]])

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

c8_df <- c8_df[order(c8_df$Log2FC, decreasing = TRUE),]
c8_df$rank<-seq_len(nrow(c8_df))
c8_df$blank<-""
c8_df<- c8_df[,c(1,2,3,9,10,5)]

c9_df <- c9_df[order(c9_df$Log2FC, decreasing = TRUE),]
c9_df$rank<-seq_len(nrow(c9_df))
c9_df$blank<-""
c9_df<- c9_df[,c(1,2,3,9,10,5)]

c10_df <- c10_df[order(c10_df$Log2FC, decreasing = TRUE),]
c10_df$rank<-seq_len(nrow(c10_df))
c10_df$blank<-""
c10_df<- c10_df[,c(1,2,3,9,10,5)]

c11_df <- c11_df[order(c11_df$Log2FC, decreasing = TRUE),]
c11_df$rank<-seq_len(nrow(c11_df))
c11_df$blank<-""
c11_df<- c11_df[,c(1,2,3,9,10,5)]

c12_df <- c12_df[order(c12_df$Log2FC, decreasing = TRUE),]
c12_df$rank<-seq_len(nrow(c12_df))
c12_df$blank<-""
c12_df<- c12_df[,c(1,2,3,9,10,5)]

c13_df <- c13_df[order(c13_df$Log2FC, decreasing = TRUE),]
c13_df$rank<-seq_len(nrow(c13_df))
c13_df$blank<-""
c13_df<- c13_df[,c(1,2,3,9,10,5)]

c14_df <- c14_df[order(c14_df$Log2FC, decreasing = TRUE),]
c14_df$rank<-seq_len(nrow(c14_df))
c14_df$blank<-""
c14_df<- c14_df[,c(1,2,3,9,10,5)]

c15_df <- c15_df[order(c15_df$Log2FC, decreasing = TRUE),]
c15_df$rank<-seq_len(nrow(c15_df))
c15_df$blank<-""
c15_df<- c15_df[,c(1,2,3,9,10,5)]

c16_df <- c16_df[order(c16_df$Log2FC, decreasing = TRUE),]
c16_df$rank<-seq_len(nrow(c16_df))
c16_df$blank<-""
c16_df<- c16_df[,c(1,2,3,9,10,5)]

c17_df <- c17_df[order(c17_df$Log2FC, decreasing = TRUE),]
c17_df$rank<-seq_len(nrow(c17_df))
c17_df$blank<-""
c17_df<- c17_df[,c(1,2,3,9,10,5)]

c18_df <- c18_df[order(c18_df$Log2FC, decreasing = TRUE),]
c18_df$rank<-seq_len(nrow(c18_df))
c18_df$blank<-""
c18_df<- c18_df[,c(1,2,3,9,10,5)]

c19_df <- c19_df[order(c19_df$Log2FC, decreasing = TRUE),]
c19_df$rank<-seq_len(nrow(c19_df))
c19_df$blank<-""
c19_df<- c19_df[,c(1,2,3,9,10,5)]

### export marker peaks
setwd("/Users/sranade/scATAC-seq/Mesoderm_subset1_210517/Exports/Clusters1_markerpeaks")
write.table(c1_df, file = "c1_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c2_df, file = "c2_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c3_df, file = "c3_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c4_df, file = "c4_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c5_df, file = "c5_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c6_df, file = "c6_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c7_df, file = "c7_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c8_df, file = "c8_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c9_df, file = "c9_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c10_df, file = "c10_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c11_df, file = "c11_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c12_df, file = "c12_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c13_df, file = "c13_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c14_df, file = "c14_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c15_df, file = "c15_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c16_df, file = "c16_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c17_df, file = "c17_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c18_df, file = "c18_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)
write.table(c19_df, file = "c19_df.bed",quote = F, sep = '\t',col.names = F, row.names = F)

### what is going on with cluster 2? compare to cluster 1
#pSHF
markerTest_clust2v1 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2",
  bgdGroups = "C1"
)
pma_C2 <- plotMarkers(seMarker = markerTest_clust2v1, name = "C2", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


## spot check GAS
markerGenes_spotmarkerGenes_spot <-c("Scx","Sema3d")
p_spotcheck <- plotEmbedding(
  ArchRProj = mesoderm_subset,
  colorBy = "GeneScoreMatrix",
  name = markerGenes_spotmarkerGenes_spot,
  continuousSet = "horizonExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(mesoderm_subset)
)
p_spotcheck$Sema3d

















saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/mesoderm_subset/mesoderm_subset_Archproject", load = TRUE)











