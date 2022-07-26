#### Chapters 9 and 10 Pseudo-bulk Replicates and Peak calling
#### 

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

######## 
# Re-load ArchR project!
wt_atlas <- loadArchRProject(path = "/Users/sranade/scATAC-seq/new_wt_section_211111")
wt_atlas
######## 

## 9. Pseudo-bulk Replicates in ArchR
table(wt_atlas$Clusters2,wt_atlas$Sample)

wt_atlas <- addGroupCoverages(ArchRProj = wt_atlas,minReplicates = 2,
                              maxReplicates = 18,
                              maxFragments = 75 * 10^6,
                              groupBy = "Clusters2",
                              minCells = 20,
                              maxCells = 3000,
                              sampleRatio = 0.8,
                              returnGroups = FALSE,
                              force = T)
# 2021-11-14 15:43:12 : Adding Kmer Bias to Coverage Files!, 74.784 mins elapsed.
# 2021-11-14 16:38:15 : Finished Creation of Coverage Files!, 129.833 mins elapsed.


## 10. Peak calling -- watch macs2 path!
wt_atlas <- addReproduciblePeakSet(
  ArchRProj = wt_atlas, 
  groupBy = "Clusters2", 
  pathToMacs2 = "/Users/sranade/opt/anaconda3/bin/macs2",
  peaksPerCell = 2500,
  maxPeaks = 200000,
  minCells = 25,
  reproducibility = "2",
  force = T
)
#                               Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
# Blood                         Blood   2004       1990          11   38  403    2e+05
# Cardiomyocyte         Cardiomyocyte   7602       7602          18   84 1050    2e+05
# Ectoderm                   Ectoderm   7270       7270          14   34 1452    2e+05
# EndMT                         EndMT    664        664           8   23  212    2e+05
# Endoderm                   Endoderm   9726       9722          17  145 1329    2e+05
# Endothelium             Endothelium   3482       3482          18   59  629    2e+05
# Epicardium               Epicardium   3861       3861          18   49  449    2e+05
# LPM                             LPM   1119       1102          12   56  153    2e+05
# Neural_Crest           Neural_Crest  10401      10391          17   26 2857    2e+05
# Paraxial_Mesoderm Paraxial_Mesoderm   6536       6533          17   35 1150    2e+05
# SHF_Progenitor       SHF_Progenitor  12291      12291          18  152 2172    2e+05

#

# look at the peak set and then save it as a separate object
getPeakSet(wt_atlas)
#GRanges object with 548312 ranges and 13 metadata columns:
merged_peak_setGR <- getPeakSet(wt_atlas)

# add peakmatrix to object
wt_atlas <- addPeakMatrix(wt_atlas)
getAvailableMatrices(wt_atlas)
#[1] "GeneScoreMatrix" "PeakMatrix"      "TileMatrix" 
# have not added the integration matrix yet!

getGroupBW(
  ArchRProj = wt_atlas,
  groupBy = "Clusters2",
  normMethod = "None",
  tileSize = 100,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)

