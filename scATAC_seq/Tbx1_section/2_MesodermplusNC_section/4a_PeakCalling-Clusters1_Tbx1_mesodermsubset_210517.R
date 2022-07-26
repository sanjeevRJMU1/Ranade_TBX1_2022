###### Tbx1 WT v KO mesoderm subset1
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
#### Chapters 9 and 10 Pseudo-bulk Replicates and Peak calling
#### This version = groupcoverage on Clusters4 = Clusters1 divided into WT and KO

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset1_210517")
mesoderm_subset
######## 


## 9. Pseudo-bulk Replicates in ArchR --> important!! Use Clusters3 to call peaks on WT and KO separately
mesoderm_subset <- addGroupCoverages(ArchRProj = mesoderm_subset, groupBy = "Clusters", minCells = 40, maxCells = 1000, force = T)
# 2021-01-24 14:24:00 : Further Sampled 21 Groups above the Max Fragments!, 0.076 mins elapsed.
#2020-10-08 10:47:51 : Finished Creation of Coverage Files!, 47.066 mins elapsed.

## 10. Peak calling
mesoderm_subset <- addReproduciblePeakSet(
  ArchRProj = mesoderm_subset, 
  groupBy = "Clusters", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 1000,
  maxPeaks = 150000,
  minCells = 25,
  reproducibility = "2",
  force = T
)

#2021-01-24 15:49:52 : Finished Creating Union Peak Set (363009)!, 59.769 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(mesoderm_subset)
merged_peak_setGR <- getPeakSet(mesoderm_subset)

# add peakmatrix to object
mesoderm_subset <- addPeakMatrix(mesoderm_subset, force = T, binarize = F, ceiling = 4)
getAvailableMatrices(mesoderm_subset)

saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/mesoderm_subset", load = TRUE)

