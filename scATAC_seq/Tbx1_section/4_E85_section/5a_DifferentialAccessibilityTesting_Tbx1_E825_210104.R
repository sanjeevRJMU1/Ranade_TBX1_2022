###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
### Chapter 11: Differential Accessibility Testing



###############################################
# II. Script part 1, running DARs on each cluster in Clusters3

table(Tbx1_E825$Clusters3)
# AHF_K                AHF_W          CM_Atrium_K          CM_Atrium_W             CM_OFT_K             CM_OFT_W 
# 182                  297                  368                  402                  127                  355 
# CM_Ventricle_K       CM_Ventricle_W                CNC_K                CNC_W           Ectoderm_K           Ectoderm_W 
# 397                  758                  259                   40                 2289                 1366 
# Endocardium_K        Endocardium_W           Endoderm_K           Endoderm_W         Epicardium_K         Epicardium_W 
# 300                  339                 1493                 1485                  514                  612 
# LPM_K                LPM_W PharyngealMesoderm_K PharyngealMesoderm_W                 PM_K                 PM_W 
# 216                  237                  496                  275                  630                  595 
# pSHF_K               pSHF_W 
# 642                  504 

### II. Script part 2, running DARs on each cluster in manually annotated Clusters3
## For all clusters, use FDR <= 0.05 ad log2FC > 1. This is prob the safe for FDR and still allows you to catch as many 
## peaks as possible since its not entirely clear what "1.8x more accessible actually means"

#AHF
markerTest_AHF <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 182,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "AHF_W",
  bgdGroups = "AHF_K"
)
pma_AHF <- plotMarkers(seMarker = markerTest_AHF, name = "AHF_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


#CM_Atrium
markerTest_CM_Atrium <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 368,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CM_Atrium_W",
  bgdGroups = "CM_Atrium_K"
)
pma_CM_Atrium <- plotMarkers(seMarker = markerTest_CM_Atrium, name = "CM_Atrium_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


#CM_OFT
markerTest_CM_OFT <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 127,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CM_OFT_W",
  bgdGroups = "CM_OFT_K"
)
pma_CM_OFT <- plotMarkers(seMarker = markerTest_CM_OFT, name = "CM_OFT_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


#CM_Ventricle
markerTest_CM_Ventricle <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 397,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CM_Ventricle_W",
  bgdGroups = "CM_Ventricle_K"
)
pma_CM_Ventricle <- plotMarkers(seMarker = markerTest_CM_Ventricle, name = "CM_Ventricle_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


#CNC
markerTest_CNC <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 40,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CNC_W",
  bgdGroups = "CNC_K"
)
pma_CNC <- plotMarkers(seMarker = markerTest_CNC, name = "CNC_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#Ectoderm
markerTest_Ectoderm <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 1366,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Ectoderm_W",
  bgdGroups = "Ectoderm_K"
)
pma_Ectoderm <- plotMarkers(seMarker = markerTest_Ectoderm, name = "Ectoderm_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#Endocardium
markerTest_Endocardium <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 300,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Endocardium_W",
  bgdGroups = "Endocardium_K"
)
pma_Endocardium <- plotMarkers(seMarker = markerTest_Endocardium, name = "Endocardium_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#Endoderm
markerTest_Endoderm <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 1485,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Endoderm_W",
  bgdGroups = "Endoderm_K"
)
pma_Endoderm <- plotMarkers(seMarker = markerTest_Endoderm, name = "Endoderm_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#Epicardium
markerTest_Epicardium <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 514,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Epicardium_W",
  bgdGroups = "Epicardium_K"
)
pma_Epicardium <- plotMarkers(seMarker = markerTest_Epicardium, name = "Epicardium_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#LPM
markerTest_LPM <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 216,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "LPM_W",
  bgdGroups = "LPM_K"
)
pma_LPM <- plotMarkers(seMarker = markerTest_LPM, name = "LPM_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#PM
markerTest_PM <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 595,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PM_W",
  bgdGroups = "PM_K"
)
pma_PM <- plotMarkers(seMarker = markerTest_PM, name = "PM_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#PharynMeso
markerTest_PharyngealMesoderm <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 275,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PharyngealMesoderm_W",
  bgdGroups = "PharyngealMesoderm_K"
)
pma_PharyngealMesoderm <- plotMarkers(seMarker = markerTest_PharyngealMesoderm, name = "PharyngealMesoderm_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


#pSHF
markerTest_pSHF <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 504,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "pSHF_W",
  bgdGroups = "pSHF_K"
)
pma_pSHF <- plotMarkers(seMarker = markerTest_pSHF, name = "pSHF_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


## Plots
pma_AHF
pma_CM_Atrium
pma_CM_Ventricle
pma_CM_OFT
pma_CNC
pma_Ectoderm
pma_Endocardium
pma_Endoderm
pma_Epicardium
pma_LPM
pma_PM
pma_PharyngealMesoderm
pma_pSHF
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



