##### Alternate peak call --> DAR workflow, starting from the point of psuedo- bulk replicates

# create a new object so you don't overwrite the other workflow
Tbx1_E925_new <- Tbx1_E925
table(Tbx1_E925_new$Clusters2,Tbx1_E925_new$Sample)

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters3
condition.vec.01 <- substr(Tbx1_E925_new$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_E925_new$Clusters2)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
# condition.vec
# AHF_K                AHF_W              Blood_K              Blood_W      Cardiomyocyte_K      Cardiomyocyte_W 
# 1047                 1327                  260                  308                  599                  722 
# CNC_K                CNC_W           Ectoderm_K           Ectoderm_W        Endocardium_K        Endocardium_W 
# 1003                  899                 2308                 2344                  398                  382 
# Endoderm_K           Endoderm_W         Epicardium_K         Epicardium_W   ParaxialMesoderm_K   ParaxialMesoderm_W 
# 1829                 2078                  310                  280                  439                  814 
# PharyngealMesoderm_K PharyngealMesoderm_W               pSHF_K               pSHF_W 
# 739                  565                 1029                  773 

Tbx1_E925_new$Clusters3 <- as.character(condition.vec)
table(Tbx1_E925_new$Clusters3)

### IMPORTANT, save as a new project. Then re-load the new project
saveArchRProject(ArchRProj = Tbx1_E925_new, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E925/Tbx1_E925_new_Archproject", load = TRUE)
Tbx1_E925_new <- loadArchRProject(path ="/Users/sranade/scATAC-seq/Tbx1_E925/Tbx1_E925_new_Archproject")
Tbx1_E925_new

## 9. Pseudo-bulk Replicates in ArchR
Tbx1_E925_new <- addGroupCoverages(ArchRProj = Tbx1_E925_new, groupBy = "Clusters3", minCells = 40, maxCells = 1500, force = T)
#020-09-28 19:35:16 : Finished Creation of Coverage Files!, 39.599 mins elapsed.

## 10. Peak calling
Tbx1_E925_new <- addReproduciblePeakSet(
  ArchRProj = Tbx1_E925_new, 
  groupBy = "Clusters3", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 25,
  force = T
)

#                                     Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
# AHF_K                               AHF_K   1047       1047           2  418  629   150000
# AHF_W                               AHF_W   1327       1327           2  594  733   150000
# Blood_K                           Blood_K    260        260           2  103  157   130000
# Blood_W                           Blood_W    308        308           2  104  204   150000
# Cardiomyocyte_K           Cardiomyocyte_K    599        599           2  288  311   150000
# Cardiomyocyte_W           Cardiomyocyte_W    722        722           2  214  508   150000
# CNC_K                               CNC_K   1003       1003           2  403  600   150000
# CNC_W                               CNC_W    899        899           2  246  653   150000
# Ectoderm_K                     Ectoderm_K   2308       2308           2  941 1367   150000
# Ectoderm_W                     Ectoderm_W   2344       2337           2  837 1500   150000
# Endocardium_K               Endocardium_K    398        398           2  188  210   150000
# Endocardium_W               Endocardium_W    382        382           2  150  232   150000
# Endoderm_K                     Endoderm_K   1829       1829           2  715 1114   150000
# Endoderm_W                     Endoderm_W   2078       2078           2  781 1297   150000
# Epicardium_K                 Epicardium_K    310        310           2  133  177   150000
# Epicardium_W                 Epicardium_W    280        280           2  115  165   140000
# ParaxialMesoderm_K     ParaxialMesoderm_K    439        439           2  184  255   150000
# ParaxialMesoderm_W     ParaxialMesoderm_W    814        814           2  337  477   150000
# PharyngealMesoderm_K PharyngealMesoderm_K    739        739           2  297  442   150000
# PharyngealMesoderm_W PharyngealMesoderm_W    565        565           2  159  406   150000
# pSHF_K                             pSHF_K   1029       1029           2  328  701   150000
# pSHF_W                             pSHF_W    773        773           2  333  440   150000




# 2020-09-28 21:59:40 : Finished Creating Union Peak Set (416062)!, 77.334 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(Tbx1_E925_new)
merged_peak_setGR <- getPeakSet(Tbx1_E925_new)

# add peakmatrix to object
Tbx1_E925_new <- addPeakMatrix(Tbx1_E925_new)
getAvailableMatrices(Tbx1_E925_new)
#[1] "GeneIntegrationMatrix" "GeneScoreMatrix"       "MotifMatrix"           "PeakMatrix"            "TileMatrix" 

table(Tbx1_E925_new$Clusters3)



# 1. AHF
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 1047,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "AHF_W",
  bgdGroups = "AHF_K"
)
pma <- markerPlot(seMarker = markerTest, name = "AHF_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 2. Blood
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 260,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Blood_W",
  bgdGroups = "Blood_K"
)
pma <- markerPlot(seMarker = markerTest, name = "Blood_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 3. Cardiomyocyte
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 599,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiomyocyte_W",
  bgdGroups = "Cardiomyocyte_K"
)
pma <- markerPlot(seMarker = markerTest, name = "Cardiomyocyte_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 4. CNC
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 899,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CNC_W",
  bgdGroups = "CNC_K"
)
pma <- markerPlot(seMarker = markerTest, name = "CNC_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 5. Ectoderm
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 2308,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Ectoderm_W",
  bgdGroups = "Ectoderm_K"
)
pma <- markerPlot(seMarker = markerTest, name = "Ectoderm_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 6. Endocardium
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 382,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Endocardium_W",
  bgdGroups = "Endocardium_K"
)
pma <- markerPlot(seMarker = markerTest, name = "Endocardium_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 7. Endoderm
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 1829,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Endoderm_W",
  bgdGroups = "Endoderm_K"
)
pma <- markerPlot(seMarker = markerTest, name = "Endoderm_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 8. Epicardium
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 280,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Epicardium_W",
  bgdGroups = "Epicardium_K"
)
pma <- markerPlot(seMarker = markerTest, name = "Epicardium_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma


# 9. ParaxialMesoderm
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 439,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ParaxialMesoderm_W",
  bgdGroups = "ParaxialMesoderm_K"
)
pma <- markerPlot(seMarker = markerTest, name = "ParaxialMesoderm_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 10. PharyngealMesoderm
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 565,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PharyngealMesoderm_W",
  bgdGroups = "PharyngealMesoderm_K"
)
pma <- markerPlot(seMarker = markerTest, name = "PharyngealMesoderm_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

# 11. pSHF
markerTest <- getMarkerFeatures(
  ArchRProj = Tbx1_E925_new, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters3",
  maxCells = 773,
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "pSHF_W",
  bgdGroups = "pSHF_K"
)
pma <- markerPlot(seMarker = markerTest, name = "pSHF_W", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma








