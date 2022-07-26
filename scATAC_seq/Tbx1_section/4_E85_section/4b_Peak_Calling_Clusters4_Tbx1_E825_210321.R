###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
#### Chapters 9 and 10 Pseudo-bulk Replicates and Peak calling
#### This version = groupcoverage on Clusters2, first separate WT v KO then call peaks

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

# divide clusters by genotype
# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters4"
table(Tbx1_E825$Sample)
condition.vec.01 <- substr(Tbx1_E825$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_E825$Clusters)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
Tbx1_E825$Clusters4 <- as.character(condition.vec)

# Confirm clusters are separated by genotype
table(Tbx1_E825$Clusters4)

## 9. Pseudo-bulk Replicates in ArchR --> important!! Use Clusters4 to call peaks on WT and KO separately
Tbx1_E825 <- addGroupCoverages(ArchRProj = Tbx1_E825, groupBy = "Clusters4", minCells = 30, maxCells = 1500, force = T)
#2021-01-04 14:59:58 : Finished Creation of Coverage Files!, 28.109 mins elapsed.

## 10. Peak calling
Tbx1_E825 <- addReproduciblePeakSet(
  ArchRProj = Tbx1_E825, 
  groupBy = "Clusters4", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 2000,
  maxPeaks = 150000,
  minCells = 25,
  reproducibility = "2",
  force = T
)

# Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
# C1_K   C1_K     95         83           2   39   44   150000
# C1_W   C1_W    182        180           3   43   79   150000
# C2_K   C2_K    325        325           3   91  135   150000
# C2_W   C2_W    458        458           4   84  173   150000
# C3_K   C3_K    238        238           3   42  121   150000
# C3_W   C3_W    353        353           4   44  113   150000
# C4_K   C4_K    124        124           2   48   76   150000
# C4_W   C4_W     83         72           2   35   37   144000
# C5_K   C5_K    677        677           3   49  419   150000
# C5_W   C5_W    381        372           3   79  165   150000
# C6_K   C6_K    191        191           3   47   92   150000
# C6_W   C6_W    112         89           2   44   45   150000
# C7_K   C7_K    107         94           2   43   51   150000
# C7_W   C7_W     27         24           2   14   17    48000
# C8_K   C8_K     48         45           2   30   30    90000
# C8_W   C8_W     79         79           2   30   49   150000
# C9_K   C9_K    632        632           3  124  281   150000
# C9_W   C9_W    457        457           3   49  237   150000
# C10_K C10_K    413        413           3   99  206   150000
# C10_W C10_W    200        194           2   53  141   150000
# C11_K C11_K    831        831           3  139  518   150000
# C11_W C11_W    442        428           3   57  217   150000
# C12_K C12_K     69         69           2   30   39   138000
# C12_W C12_W     47         45           2   30   30    90000
# C13_K C13_K    564        564           3  122  263   150000
# C13_W C13_W    560        560           3   88  307   150000
# C14_K C14_K    555        555           3  112  327   150000
# C14_W C14_W     17         17           2   11   14    34000
# C15_K C15_K      1          1           2    1    1     2000
# C15_W C15_W    367        367           4   34  161   150000
# C16_K C16_K     82         82           2   36   46   150000
# C16_W C16_W     82         82           2   41   41   150000
# C17_K C17_K    256        235           2   40  195   150000
# C17_W C17_W     40         40           2   30   30    80000
# C18_K C18_K     38         38           2   30   30    76000
# C18_W C18_W     37         34           2   21   21    68000
# C19_K C19_K    299        299           3   66  135   150000
# C19_W C19_W    334        334           4   61  116   150000
# C20_K C20_K    341        341           3   87  130   150000
# C20_W C20_W    367        367           4   38  158   150000
# C21_K C21_K    418        418           3  102  203   150000
# C21_W C21_W    780        780           4  107  272   150000
# C22_K C22_K     97         97           3   31   35   150000
# C22_W C22_W    212        212           4   40   65   150000
# C23_K C23_K    663        663           3  108  314   150000
# C23_W C23_W    385        377           3   40  214   150000
# C24_K C24_K    222        222           3   51  104   150000
# C24_W C24_W    195        193           3   39  106   150000
# C25_K C25_K     75         75           1   75   75   150000
# C25_W C25_W    458        458           4   81  175   150000
# C26_K C26_K    230        230           3   36  148   150000
# C26_W C26_W    249        248           3   56  113   150000
# C27_K C27_K    322        322           3   63  150   150000
# C27_W C27_W    361        353           3   89  145   150000

# 2021-03-21 16:19:38 : Batching Peak Calls!, 0.008 mins elapsed.
# 2021-03-21 16:19:38 : Batch Execution w/ safelapply!, 0 mins elapsed.
# 2021-03-21 17:57:54 : Identifying Reproducible Peaks!, 98.266 mins elapsed.
# 2021-03-21 18:00:43 : Creating Union Peak Set!, 101.093 mins elapsed.
# Converged after 14 iterations!
#   Plotting Ggplot!
#   2021-03-21 18:01:26 : Finished Creating Union Peak Set (481439)!, 101.802 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(Tbx1_E825)
merged_peak_setGR <- getPeakSet(Tbx1_E825)

# add peakmatrix to object
Tbx1_E825 <- addPeakMatrix(Tbx1_E825, force = T, binarize = F, ceiling = 4)
getAvailableMatrices(Tbx1_E825)


# get the groupBW files
#getGroupBW(ArchRProj = Tbx1_E825, groupBy = "Clusters4", tileSize = 100, ceiling = 4, normMethod = "None")

saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825", load = TRUE)

