###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
#### Chapters 9 and 10 Pseudo-bulk Replicates and Peak calling
#### This version = groupcoverage on Clusters2, which calls peaks on WT and KO separately. This is nice bc it gives you 
#### browser tracks for both genotypes

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 0.9.5


######## 
# Re-load ArchR project!
Tbx1_E825 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825")
Tbx1_E825
######## 
# Confirm clusters are separated by genotype
table(Tbx1_E825$Clusters2)

## 9. Pseudo-bulk Replicates in ArchR --> important!! Use Clusters2 to call peaks on WT and KO separately
Tbx1_E825 <- addGroupCoverages(ArchRProj = Tbx1_E825, groupBy = "Clusters2", minCells = 30, maxCells = 1500, force = T)
#2021-01-04 14:59:58 : Finished Creation of Coverage Files!, 28.109 mins elapsed.

## 10. Peak calling
Tbx1_E825 <- addReproduciblePeakSet(
  ArchRProj = Tbx1_E825, 
  groupBy = "Clusters2", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 25,
  reproducibility = "2",
  force = T
)

#2020-10-08 12:17:18 : Creating Union Peak Set!, 89.441 mins elapsed.
#2021-01-04 16:09:50 : Finished Creating Union Peak Set (404481)!, 69.857 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(Tbx1_E825)
merged_peak_setGR <- getPeakSet(Tbx1_E825)

# add peakmatrix to object
Tbx1_E825 <- addPeakMatrix(Tbx1_E825, force = T, binarize = F, ceiling = 4)
getAvailableMatrices(Tbx1_E825)


# get the groupBW files
Tbx1_E825 <- getGroupBW(ArchRProj = Tbx1_E825, groupBy = "Clusters2", tileSize = 100, ceiling = 4, normMethod = "None")

saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/SR39_Tbx1_E825", load = TRUE)

