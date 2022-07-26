###### Tbx1 WT v KO mesoderm subset
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
##### All scripts for making figures

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

p1 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p9,p10, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = mesoderm_subset, addDOC = FALSE, width = 5, height = 5)


###### 
# Re-naming idents
## My manual labeling --> this goes to Clusters6 as rest are already taken
library(plyr)
temp_cluster_names <- as.character(revalue(mesoderm_subset$Clusters, c(
  "C1"="Cardiomyocyte",
  "C2"="Cardiomyocyte",
  "C3"="Cardiopharyngeal_Mesenchyme",
  "C4"="Cranial_Mesenchyme",
  "C5"="Neural_Crest",
  "C6"="Neural_Crest",
  "C7"="Neural_Crest",
  "C8"="Epicardium",
  "C9"="Smooth_Muscle",
  "C10"="Neural_Crest",
  "C11"="Cardiopharyngeal_Mesenchyme",
  "C12"="Paraxial_Mesoderm",
  "C13"="Cardiopharyngeal_Mesoderm",
  "C14"="Paraxial_Mesoderm",
  "C15"="Posterior_Second_Heart_Field",
  "C16"="Anterior_Second_Heart_Field"
)))
length(unique(temp_cluster_names))
mesoderm_subset$Clusters6 <- temp_cluster_names
mesoderm_subset@cellColData@listData[["Clusters6"]] <- temp_cluster_names
table(mesoderm_subset$Clusters6,mesoderm_subset$Sample)

p3 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Clusters6", embedding = "UMAPHarmony")


p3

plotPDF(p3,
        name = "Plot-UMAP-Manual-Annotations_2.pdf",
        ArchRProj = mesoderm_subset,
        addDOC = FALSE, width = , height = 5)


saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517", load = TRUE)


########################




