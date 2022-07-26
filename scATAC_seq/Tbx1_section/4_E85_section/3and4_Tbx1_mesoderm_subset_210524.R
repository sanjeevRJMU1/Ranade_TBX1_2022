mesoderm_subset <- addImputeWeights(mesoderm_subset, reducedDims = "Harmony", dimsToUse = 2:25)

markerGenes <- c("Tnnt2","Isl1","Osr1","Fgf8","Tbx1","Wt1","Dlx2","Mab21l2","Meox1",
                 "Tdgf1","Rgs5","Fgf10","Irx4","Vsnl1","Cited1","Nr2f1","Foxd1","Foxc2","Alx1","Tlx1","Aplnr","Nrg1","Hand1","Tbx18","Foxf1","Bmp2","Rspo3")

p13 <- plotEmbedding(
  ArchRProj = mesoderm_subset, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(mesoderm_subset)
)
plotPDF(plotList = p13, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = mesoderm_subset, 
        addDOC = FALSE, width = 5, height = 5)

markerGenes_2 <- c("Nkx2-5","Tbx5","Id2")

p13 <- plotEmbedding(
  ArchRProj = mesoderm_subset, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes_2, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(mesoderm_subset)
)
plotPDF(plotList = p13, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation_2.pdf", 
        ArchRProj = mesoderm_subset, 
        addDOC = FALSE, width = 5, height = 5)



# call peaks on clusters1, wt and ko

table(mesoderm_subset$Clusters,mesoderm_subset$Sample)

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters4"
condition.vec.01 <- substr(mesoderm_subset$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(mesoderm_subset$Clusters)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
# condition.vec
# C1_K  C1_W C10_K C10_W C11_K C11_W C12_K C12_W C13_K C13_W  C2_K  C2_W  C3_K  C3_W  C4_K  C4_W  C5_K  C5_W  C6_K  C6_W  C7_K 
# 213   228   413   257   167   172   587   573   205   282    69    49    89    89   187   257   498   565   395   334   269 
# C7_W  C8_K  C8_W  C9_K  C9_W 
# 258   391   698   166   318  

mesoderm_subset$Clusters4 <- as.character(condition.vec)
table(mesoderm_subset$Clusters4)


## 9. Pseudo-bulk Replicates in ArchR
mesoderm_subset <- addGroupCoverages(ArchRProj = mesoderm_subset, groupBy = "Clusters4", minCells = 40, maxCells = 1500, force = T)
#2021-05-24 18:57:51 : Finished Creation of Coverage Files!, 45.376 mins elapsed.

## 10. Peak calling
mesoderm_subset <- addReproduciblePeakSet(
  ArchRProj = mesoderm_subset, 
  groupBy = "Clusters4", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 1500,
  maxPeaks = 125000,
  minCells = 25,
  force = T
)


# look at the peak set and then save it as a separate object
getPeakSet(mesoderm_subset)
merged_peak_setGR <- getPeakSet(mesoderm_subset)

# add peakmatrix to object
mesoderm_subset <- addPeakMatrix(mesoderm_subset)
getAvailableMatrices(mesoderm_subset)

########################################################
# DAR per cluster
###############################################
# II. Script part 1, running DARs on each cluster

table(mesoderm_subset$Clusters4)

### II. Script part 2, running DARs on each cluster in manually annotated Clusters3
## For all clusters, use FDR <= 0.05 ad log2FC > 1. This is prob the safe for FDR and still allows you to catch as many 
## peaks as possible since its not entirely clear what "1.8x more accessible actually means"

#
markerTest_C1 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1_W",
  bgdGroups = "C1_K"
)
pma_C1 <- plotMarkers(seMarker = markerTest_C1, name = "C1_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C2 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2_W",
  bgdGroups = "C2_K"
)
pma_C2 <- plotMarkers(seMarker = markerTest_C2, name = "C2_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C3 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C3_W",
  bgdGroups = "C3_K"
)
pma_C3 <- plotMarkers(seMarker = markerTest_C3, name = "C3_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C4 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4_W",
  bgdGroups = "C4_K"
)
pma_C4 <- plotMarkers(seMarker = markerTest_C4, name = "C4_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C5 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_W",
  bgdGroups = "C5_K"
)
pma_C5 <- plotMarkers(seMarker = markerTest_C5, name = "C5_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C6 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_W",
  bgdGroups = "C6_K"
)
pma_C6 <- plotMarkers(seMarker = markerTest_C6, name = "C6_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C7 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C7_W",
  bgdGroups = "C7_K"
)
pma_C7 <- plotMarkers(seMarker = markerTest_C7, name = "C7_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C8 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C8_W",
  bgdGroups = "C8_K"
)
pma_C8 <- plotMarkers(seMarker = markerTest_C8, name = "C8_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C9 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C9_W",
  bgdGroups = "C9_K"
)
pma_C9 <- plotMarkers(seMarker = markerTest_C9, name = "C9_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C10 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C10_W",
  bgdGroups = "C10_K"
)
pma_C10 <- plotMarkers(seMarker = markerTest_C10, name = "C10_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C11 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C11_W",
  bgdGroups = "C11_K"
)
pma_C11 <- plotMarkers(seMarker = markerTest_C11, name = "C11_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C12 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C12_W",
  bgdGroups = "C12_K"
)
pma_C12 <- plotMarkers(seMarker = markerTest_C12, name = "C12_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C13 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C13_W",
  bgdGroups = "C13_K"
)
pma_C13 <- plotMarkers(seMarker = markerTest_C13, name = "C13_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210524", load = TRUE)


