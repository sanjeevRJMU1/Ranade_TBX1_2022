###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
### Chapter 11: Differential Accessibility Testing



###############################################
# II. Script part 1, running DARs on each cluster in Clusters3

table(Tbx1_E825$Clusters4)

### II. Script part 2, running DARs on each cluster in manually annotated Clusters3
## For all clusters, use FDR <= 0.05 ad log2FC > 1. This is prob the safe for FDR and still allows you to catch as many 
## peaks as possible since its not entirely clear what "1.8x more accessible actually means"

#
markerTest_C1 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
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
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C13_W",
  bgdGroups = "C13_K"
)
pma_C13 <- plotMarkers(seMarker = markerTest_C13, name = "C13_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C14 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C14_W",
  bgdGroups = "C14_K"
)
pma_C14 <- plotMarkers(seMarker = markerTest_C14, name = "C14_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C15 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C15_W",
  bgdGroups = "C15_K"
)
pma_C15 <- plotMarkers(seMarker = markerTest_C15, name = "C15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C16 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C16_W",
  bgdGroups = "C16_K"
)
pma_C16 <- plotMarkers(seMarker = markerTest_C16, name = "C16_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C17 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C17_W",
  bgdGroups = "C17_K"
)
pma_C17 <- plotMarkers(seMarker = markerTest_C17, name = "C17_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C18 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C18_W",
  bgdGroups = "C18_K"
)
pma_C18 <- plotMarkers(seMarker = markerTest_C18, name = "C18_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C19 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C19_W",
  bgdGroups = "C19_K"
)
pma_C19 <- plotMarkers(seMarker = markerTest_C19, name = "C19_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C20 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C20_W",
  bgdGroups = "C20_K"
)
pma_C20 <- plotMarkers(seMarker = markerTest_C20, name = "C20_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C21 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C21_W",
  bgdGroups = "C21_K"
)
pma_C21 <- plotMarkers(seMarker = markerTest_C21, name = "C21_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C22 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C22_W",
  bgdGroups = "C22_K"
)
pma_C22 <- plotMarkers(seMarker = markerTest_C22, name = "C22_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C23 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C23_W",
  bgdGroups = "C23_K"
)
pma_C23 <- plotMarkers(seMarker = markerTest_C23, name = "C23_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C24 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C24_W",
  bgdGroups = "C24_K"
)
pma_C24 <- plotMarkers(seMarker = markerTest_C24, name = "C24_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C25 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C25_W",
  bgdGroups = "C25_K"
)
pma_C25 <- plotMarkers(seMarker = markerTest_C25, name = "C25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C26 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C26_W",
  bgdGroups = "C26_K"
)
pma_C26 <- plotMarkers(seMarker = markerTest_C26, name = "C26_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C27 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C27_W",
  bgdGroups = "C27_K"
)
pma_C27 <- plotMarkers(seMarker = markerTest_C27, name = "C27_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


# export PDF into illustrator


## II. Script part 3: Extract the differential peaks for the clusters you want using the same sig cut off for plotmarkers
## choose FDR 0.05 and logFC 1 

# AHF
markerList_AHF <- getMarkers(markerTest_AHF, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_AHF_df <- as.data.frame(markerList_AHF@listData[["AHF_W"]])
write.csv(markerList_AHF_df, file = "/Users/sranade/scATAC-seq/Tbx1_E825/DAR_homer/DAR_csv/AHF.csv")

# Pharyngeal Mesoderm
markerList_PharyngealMesoderm <- getMarkers(markerTest_PharyngealMesoderm, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_PharyngealMesoderm_df <- as.data.frame(markerList_PharyngealMesoderm@listData[["PharyngealMesoderm_W"]])
write.csv(markerList_PharyngealMesoderm_df, file = "/Users/sranade/scATAC-seq/Tbx1_E825/DAR_homer/DAR_csv/PharyngealMesoderm.csv")

# Endoderm
markerList_Endoderm <- getMarkers(markerTest_Endoderm, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Endoderm_df <- as.data.frame(markerList_Endoderm@listData[["Endoderm_W"]])
write.csv(markerList_Endoderm_df, file = "/Users/sranade/scATAC-seq/Tbx1_E825/DAR_homer/DAR_csv/Endoderm.csv")

# Ectoderm
markerList_Ectoderm <- getMarkers(markerTest_Ectoderm, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Ectoderm_df <- as.data.frame(markerList_Ectoderm@listData[["Ectoderm_W"]])
write.csv(markerList_Ectoderm_df, file = "/Users/sranade/scATAC-seq/Tbx1_E825/DAR_homer/DAR_csv/Ectoderm.csv")

#################
## Extract the differential peaks for the clusters for peaks that are MORE accessible in KO

# AHF
markerList_AHF <- getMarkers(markerTest_AHF, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_AHF_df <- as.data.frame(markerList_AHF@listData[["AHF_W"]])
write.csv(markerList_AHF_df, file = "/Users/sranade/scATAC-seq/Tbx1_E825/DAR_homer/DAR_csv/upKO/AHF_upKO.csv")

# Pharyngeal Mesoderm
markerList_PharyngealMesoderm <- getMarkers(markerTest_PharyngealMesoderm, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_PharyngealMesoderm_df <- as.data.frame(markerList_PharyngealMesoderm@listData[["PharyngealMesoderm_W"]])
write.csv(markerList_PharyngealMesoderm_df, file = "/Users/sranade/scATAC-seq/Tbx1_E825/DAR_homer/DAR_csv/upKO/PharyngealMesoderm_upKO.csv")

# Endoderm
markerList_Endoderm <- getMarkers(markerTest_Endoderm, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Endoderm_df <- as.data.frame(markerList_Endoderm@listData[["Endoderm_W"]])
write.csv(markerList_Endoderm_df, file = "/Users/sranade/scATAC-seq/Tbx1_E825/DAR_homer/DAR_csv/upKO/Endoderm_upKO.csv")

# Ectoderm
markerList_Ectoderm <- getMarkers(markerTest_Ectoderm, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Ectoderm_df <- as.data.frame(markerList_Ectoderm@listData[["Ectoderm_W"]])
write.csv(markerList_Ectoderm_df, file = "/Users/sranade/scATAC-seq/Tbx1_E825/DAR_homer/DAR_csv/upKO/Ectoderm_upKO.csv")

#### Proceed to script 5e for formating into Homer



