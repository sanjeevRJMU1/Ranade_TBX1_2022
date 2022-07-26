##### Alternate peak call --> DAR workflow, starting from the point of psuedo- bulk replicates

# create a new object so you don't overwrite the other workflow
table(mesoderm_subset$Clusters,mesoderm_subset$Sample)

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters3
condition.vec.01 <- substr(mesoderm_subset$Sample,5,8)
unique(condition.vec.01)
condition.vec.02 <- as.vector(mesoderm_subset$Clusters)
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

mesoderm_subset$Clusters5 <- as.character(condition.vec)
table(mesoderm_subset$Clusters5)


## 9. Pseudo-bulk Replicates in ArchR
mesoderm_subset <- addGroupCoverages(ArchRProj = mesoderm_subset, groupBy = "Clusters5", minCells = 40, maxCells = 1500, force = T)
#020-09-28 19:35:16 : Finished Creation of Coverage Files!, 39.599 mins elapsed.

## 10. Peak calling
mesoderm_subset <- addReproduciblePeakSet(
  ArchRProj = mesoderm_subset, 
  groupBy = "Clusters5", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 1500,
  maxPeaks = 150000,
  minCells = 25,
  force = T
)

# 2020-09-28 21:59:40 : Finished Creating Union Peak Set (416062)!, 77.334 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(mesoderm_subset)
merged_peak_setGR <- getPeakSet(mesoderm_subset)

# add peakmatrix to object
mesoderm_subset <- addPeakMatrix(mesoderm_subset)
getAvailableMatrices(mesoderm_subset)
#[1] "GeneIntegrationMatrix" "GeneScoreMatrix"       "MotifMatrix"           "PeakMatrix"            "TileMatrix" 










