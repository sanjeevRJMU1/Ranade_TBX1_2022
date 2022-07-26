
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1

######## 
# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517")
mesoderm_subset
######## 
getPeakSet(mesoderm_subset)
#GRanges object with 475450 ranges and 13 metadata columns:

p1 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Clusters6", embedding = "UMAPHarmony")
p1

#######
table(mesoderm_subset$Clusters8,mesoderm_subset$Sample)
