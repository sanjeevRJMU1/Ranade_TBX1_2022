###### Tbx1 E9.25 scATAC-seq. 2 WT (Embryos 1 and 4) and 2 KO (Embryos 3 and 5)
#### Chapters 9 and 10 Pseudo-bulk Replicates and Peak calling
#### This version = groupcoverage on Clusters3
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

######## 
# Re-load ArchR project!
Tbx1_E925 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_E925/Tbx1_E925_Archproject")
Tbx1_E925
######## 

# divide clusters by genotype
# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters3"
condition.vec.01 <- substr(Tbx1_E925$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_E925$Clusters2)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
# condition.vec
# AHF_K                AHF_W              Blood_K              Blood_W      Cardiomyocyte_K      Cardiomyocyte_W 
# 580                  861                  260                  308                  599                  722 
# CNC_K                CNC_W           Ectoderm_K           Ectoderm_W        Endocardium_K        Endocardium_W 
# 1003                  899                 2308                 2344                  398                  382 
# Endoderm_K           Endoderm_W         Epicardium_K         Epicardium_W   ParaxialMesoderm_K   ParaxialMesoderm_W 
# 1829                 2078                  310                  280                  439                  814 
# PharyngealMesoderm_K PharyngealMesoderm_W               pSHF_K               pSHF_W 
# 739                  565                 1496                 1239 

Tbx1_E925$Clusters3 <- as.character(condition.vec)
table(Tbx1_E925$Clusters3)

## 9. Pseudo-bulk Replicates in ArchR --> important!! Use Clusters3 to call peaks on WT and KO separately
Tbx1_E925 <- addGroupCoverages(ArchRProj = Tbx1_E925, groupBy = "Clusters3", minCells = 40, maxCells = 1500, force = T)
#2020-10-08 10:47:51 : Finished Creation of Coverage Files!, 47.066 mins elapsed.

## 10. Peak calling
Tbx1_E925 <- addReproduciblePeakSet(
  ArchRProj = Tbx1_E925, 
  groupBy = "Clusters3", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 25,
  reproducibility = "(n+1)/2",
  force = T
)

#2020-10-08 12:17:18 : Creating Union Peak Set!, 89.441 mins elapsed.
#2020-10-08 12:17:57 : Finished Creating Union Peak Set (468114)!, 90.096 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(Tbx1_E925)
merged_peak_setGR <- getPeakSet(Tbx1_E925)

# add peakmatrix to object
Tbx1_E925 <- addPeakMatrix(Tbx1_E925, force = T, binarize = F, ceiling = 4)
getAvailableMatrices(Tbx1_E925)


# get the groupBW files
Tbx1_E925 <- getGroupBW(ArchRProj = Tbx1_E925, groupBy = "Clusters3", tileSize = 100, ceiling = 4, normMethod = "None")

saveArchRProject(ArchRProj = Tbx1_E925, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E925/Tbx1_E925_Archproject", load = TRUE)

