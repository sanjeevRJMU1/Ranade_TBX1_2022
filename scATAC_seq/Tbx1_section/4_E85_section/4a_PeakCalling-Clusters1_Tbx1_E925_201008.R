###### Tbx1 E9.25 scATAC-seq. 2 WT (Embryos 1 and 4) and 2 KO (Embryos 3 and 5)
#### Chapters 9 and 10 Pseudo-bulk Replicates and Peak calling
#### This version = groupcoverage on Clusters4 = Clusters1 divided into WT and KO
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
# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters4"
condition.vec.01 <- substr(Tbx1_E925$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_E925$Clusters)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
# condition.vec
# C1_K  C1_W C10_K C10_W C11_K C11_W C12_K C12_W C13_K C13_W C14_K C14_W C15_K C15_W C16_K C16_W C17_K C17_W C18_K C18_W C19_K C19_W 
# 260   308   498   396   323   489   182    14   310   285    88    97   310   280   599   722   288   411   423   331   292   450 
# C2_K  C2_W C20_K C20_W C21_K C21_W C22_K C22_W C23_K C23_W C24_K C24_W C25_K C25_W C26_K C26_W C27_K C27_W C28_K C28_W C29_K C29_W 
# 254   287   467   466   527   351   239   151   500   414   439   814    18    24    79    91   361   317    82    77   250   504 
# C3_K  C3_W C30_K C30_W C31_K C31_W C32_K C32_W  C4_K  C4_W  C5_K  C5_W  C6_K  C6_W  C7_K  C7_W  C8_K  C8_W  C9_K  C9_W 
# 461   827   150   146   274   386   460   296   520   617   350    33    70    71   397    87   256   422   234   328 

Tbx1_E925$Clusters4 <- as.character(condition.vec)
table(Tbx1_E925$Clusters4)

## 9. Pseudo-bulk Replicates in ArchR --> important!! Use Clusters3 to call peaks on WT and KO separately
Tbx1_E925 <- addGroupCoverages(ArchRProj = Tbx1_E925, groupBy = "Clusters4", minCells = 40, maxCells = 1000, force = T)
#2020-10-08 10:47:51 : Finished Creation of Coverage Files!, 47.066 mins elapsed.

## 10. Peak calling
Tbx1_E925 <- addReproduciblePeakSet(
  ArchRProj = Tbx1_E925, 
  groupBy = "Clusters4", 
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

saveArchRProject(ArchRProj = Tbx1_E925, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E925/Tbx1_E925_Archproject", load = TRUE)

