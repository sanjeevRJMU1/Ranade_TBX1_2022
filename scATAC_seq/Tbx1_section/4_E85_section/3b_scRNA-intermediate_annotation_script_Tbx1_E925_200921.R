#### This is a temp script that fits between #2 and #3 ONLY for re-clustering
## run scripts 2 and 3 until the spot checking part, then run this temp script and then continue with #4
## bc its a temp script, I do not include the standard loading and closing lines. This is not meant to be the part
## where you start/re-start. Also, all plot saving lines are moved to the end so that you only save it once.

### 4. Iterative LSI (TF-IDF and SVD)
Tbx1_E925 <- addIterativeLSI(
  ArchRProj = Tbx1_E925,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(1.0), 
    sampleCells = 3500, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:20,
  LSIMethod = 2,
  seed = 1,
  force = T
)

#### 5. Clustering
Tbx1_E925 <- addClusters(
  input = Tbx1_E925,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1.0, force = T
)


#### 6. Single-cell Embeddings
### UMAP with IterativeLSI
Tbx1_E925 <- addUMAP(
  ArchRProj = Tbx1_E925, 
  reducedDims = "IterativeLSI", 
  name = "UMAP",
  metric = "cosine",
  force = T
)
plotEmbedding(ArchRProj = Tbx1_E925, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")




p14 <- plotEmbedding(ArchRProj = Tbx1_E925, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p15 <- plotEmbedding(ArchRProj = Tbx1_E925, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p14, p15, type = "h")


table(Tbx1_E925$Clusters,Tbx1_E925$Sample)

########################
########################
########################


##8.3 GAS spot sheck and manual assignment
markerGenes_spotmarkerGenes_spot <-c("Isl1","Hand2")
p_spotcheck <- plotEmbedding(
  ArchRProj = Tbx1_E925,
  colorBy = "GeneScoreMatrix",
  name = markerGenes_spotmarkerGenes_spot,
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Tbx1_E925)
)
p_spotcheck$Hand2


# clusters
p10

## My manual labeling --> make sure this is a separate metadata column from existing clusters
library(plyr)
temp_cluster_names <- as.character(revalue(Tbx1_E925$Clusters, c(
  "C1"="Blood",
  "C2"="Ectoderm",
  "C3"="Ectoderm",
  "C4"="Ectoderm",
  "C5"="Ectoderm",
  "C6"="Ectoderm",
  "C7"="Ectoderm",
  "C8"="Ectoderm",
  "C9"="Endoderm",
  "C10"="CNC",
  "C11"="CNC",
  "C12"="CNC",
  "C13"="Endocardium",
  "C14"="Endocardium",
  "C15"="Epicardium",
  "C16"="Cardiomyocyte",
  "C17"="AHF",
  "C18"="pSHF",
  "C19"="AHF",
  "C20"="pSHF",
  "C21"="pSHF",
  "C22"="PharyngealMesoderm",
  "C23"="PharyngealMesoderm",
  "C24"="ParaxialMesoderm",
  "C25"="Endoderm",
  "C26"="pSHF",
  "C27"="Endoderm",
  "C28"="Endoderm",
  "C29"="Endoderm",
  "C30"="Endoderm",
  "C31"="Endoderm",
  "C32"="Endoderm"
)))
length(unique(temp_cluster_names))
Tbx1_E925$Clusters2 <- temp_cluster_names
Tbx1_E925@cellColData@listData[["Clusters2"]] <- temp_cluster_names
table(Tbx1_E925$Clusters2,Tbx1_E925$Sample)

#                     a_WT_1 b_WT_2 c_KO_1 d_KO_2
# AHF                   541    320    268    312
# Blood                 104    204    103    157
# Cardiomyocyte         508    214    288    311
# CNC                   653    246    403    600
# Ectoderm             1507    837    941   1367
# Endocardium           232    150    188    210
# Endoderm             1297    781    715   1114
# Epicardium            165    115    133    177
# ParaxialMesoderm      337    477    184    255
# PharyngealMesoderm    406    159    297    442
# pSHF                  525    714    478   1018

p16 <- plotEmbedding(
  Tbx1_E925,
  colorBy = "cellColData",
  name = "Clusters2",
  embedding = "UMAP",
  pal = pal)

p16

### save all intermediate plots
p14 <- plotEmbedding(ArchRProj = Tbx1_E925, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p15 <- plotEmbedding(ArchRProj = Tbx1_E925, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p14, p15, type = "h")
plotPDF(p14,p15, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = Tbx1_E925, addDOC = FALSE, width = 5, height = 5)



plotPDF(p16,
        name = "Plot-UMAP-Manual-Annotations.pdf",
        ArchRProj = Tbx1_E925,
        addDOC = FALSE, width = , height = 5)

########################
## In this workflow, separate WT and KO per cluster and add a new metadata column called Clusters3
# create a new object so you don't overwrite the other workflow

table(Tbx1_E925$Clusters2,Tbx1_E925$Sample)

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters3
condition.vec.01 <- substr(Tbx1_E925$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_E925$Clusters2)
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

Tbx1_E925$Clusters3 <- as.character(condition.vec)
table(Tbx1_E925$Clusters3)

### IMPORTANT, save as a new project. Then re-load the new project
saveArchRProject(ArchRProj = Tbx1_E925, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E925/Tbx1_E925_Archproject", load = TRUE)



########################
########################
# Proceed to step 4


