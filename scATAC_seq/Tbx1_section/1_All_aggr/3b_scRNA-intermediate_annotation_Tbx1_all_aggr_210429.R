#### This is a temp script that fits between #2 and #3 ONLY for re-clustering
## run scripts 2 and 3 until the spot checking part, then run this temp script and then continue with #4
## bc its a temp script, I do not include the standard loading and closing lines. This is not meant to be the part
## where you start/re-start. Also, all plot saving lines are moved to the end so that you only save it once.

### 4. Iterative LSI (TF-IDF and SVD)
Tbx1_all_aggr <- addIterativeLSI(
  ArchRProj = Tbx1_all_aggr,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 3, 
  clusterParams = list(
    resolution = c(0.1,0.2),
    sampleCells = 10000,
    maxClusters = 40,
    n.start = 10), 
  varFeatures = 25000, 
  sampleCellsPre = 10000,
  dimsToUse = 2:20,
  LSIMethod = 2,
  seed = 1,
  force = T
)

# # # # Batch correction with Harmony --> from now on, use "Harmony" for reducedDims and "UMAP" for embedding
Tbx1_all_aggr <- addHarmony(
  ArchRProj = Tbx1_all_aggr,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  dimsToUse = 2:20,
  groupBy = "Sample", force = T
)

#### 5. Clustering
Tbx1_all_aggr <- addClusters(
  input = Tbx1_all_aggr,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  maxClusters = 40,
  dimsToUse = 2:20,
  resolution = 0.6, force = T
)

#### 6. Single-cell Embeddings
### UMAP with Harmony
Tbx1_all_aggr <- addUMAP(
  ArchRProj = Tbx1_all_aggr, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony",
  metric = "cosine",
  force = T
)

plotEmbedding(ArchRProj = Tbx1_all_aggr, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")



table(Tbx1_all_aggr$Clusters,Tbx1_all_aggr$Sample)
cluster_nums <- table(Tbx1_all_aggr$Clusters,Tbx1_all_aggr$Sample)
write.csv(cluster_nums, file = "/Users/sranade/scATAC-seq/Tbx1_all_aggr/cluster_numbers.csv")

########################
########################
########################


##8.3 GAS spot sheck and manual assignment
markerGenes_spotmarkerGenes_spot <-c("Dlx1","Dlx2","Dlx3","Dlx4","Dlx5","Dlx6")
p_spotcheck <- plotEmbedding(
  ArchRProj = Tbx1_all_aggr,
  colorBy = "GeneScoreMatrix",
  name = markerGenes_spotmarkerGenes_spot,
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Tbx1_all_aggr)
)
p_spotcheck$Dlx2


# clusters
p10

## My manual labeling --> make sure this is a separate metadata column from existing clusters
library(plyr)
temp_cluster_names <- as.character(revalue(Tbx1_all_aggr$Clusters, c(
  "C1"="Blood",
  "C2"="NeuralCrest",
  "C3"="NeuralCrest",
  "C4"="NeuralCrest",
  "C5"="Sox2Ectoderm",
  "C6"="Sox2Ectoderm",
  "C7"="Sox2Ectoderm",
  "C8"="Sox2Ectoderm",
  "C9"="Sox2Ectoderm",
  "C10"="Sox2Ectoderm",
  "C11"="Sox2Ectoderm",
  "C12"="Sox2Ectoderm",
  "C13"="Sox2Ectoderm",
  "C14"="EpcamEndoderm",
  "C15"="EpcamEndoderm",
  "C16"="EpcamEndoderm",
  "C17"="EpcamEndoderm",
  "C18"="EpcamEndoderm",
  "C19"="EpcamEndoderm",
  "C20"="EpcamEndoderm",
  "C21"="MesodermalProgenitors",
  "C22"="MesodermalProgenitors",
  "C23"="MesodermalProgenitors",
  "C24"="MesodermalProgenitors",
  "C25"="Cardiomyocyte",
  "C26"="MesodermalProgenitors",
  "C27"="MesodermalProgenitors",
  "C28"="Sox2Ectoderm",
  "C29"="MesodermalProgenitors",
  "C30"="MesodermalProgenitors",
  "C31"="NeuralCrest",
  "C32"="NeuralCrest",
  "C33"="MesodermalProgenitors",
  "C34"="MesodermalProgenitors",
  "C35"="MesodermalProgenitors",
  "C36"="EpcamEndoderm",
  "C37"="EpcamEndoderm",
  "C38"="Endothelium",
  "C39"="Blood",
  "C40"="Blood"
)))
length(unique(temp_cluster_names))
Tbx1_all_aggr$Clusters2 <- temp_cluster_names
Tbx1_all_aggr@cellColData@listData[["Clusters2"]] <- temp_cluster_names
table(Tbx1_all_aggr$Clusters2,Tbx1_all_aggr$Sample)

# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E825_KO_1 f_E825_WT_1 g_E825_KO_2 h_E825_KO_3 i_E825_WT_2
# Blood                         106         106         207         159           0           2           7           0           1
# Cardiomyocyte                 522         291         220         318         115         306         204         127         206
# Endothelium                   234         195         149         216          67         118         121          97          60
# EpcamEndoderm                1288         704         770        1097         292         559         816         376         210
# MesodermalProgenitors        1965        1337        1772        2180         650        1422        1552         937         338
# NeuralCrest                   657         416         254         614          21          24         204          60           0
# Sox2Ectoderm                 1503         949         845        1379         445         756        1188         634          97
# 
# j_E825_WT_3 k_E825_WT_4 l_E105_KO_1 m_E105_KO_2 n_E105_WT_1 o_E105_WT_2
# Blood                           0           0         260         337         401         394
# Cardiomyocyte                 120         239          86          53          36          30
# Endothelium                    79          77         222         198         222         179
# EpcamEndoderm                 398         322         804         628         780         575
# MesodermalProgenitors         842         566        2421        2139        2262        1838
# NeuralCrest                     4          16        1262        1813        1840        1556
# Sox2Ectoderm                  445          58         504        1771        1131        1053

p16 <- plotEmbedding(
  Tbx1_all_aggr,
  colorBy = "cellColData",
  name = "Clusters2",
  embedding = "UMAP",
  pal = pal)

p16

### save all intermediate plots
p14 <- plotEmbedding(ArchRProj = Tbx1_all_aggr, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p15 <- plotEmbedding(ArchRProj = Tbx1_all_aggr, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p14, p15, type = "h")
plotPDF(p14,p15, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = Tbx1_all_aggr, addDOC = FALSE, width = 5, height = 5)



plotPDF(p16,
        name = "Plot-UMAP-Manual-Annotations.pdf",
        ArchRProj = Tbx1_all_aggr,
        addDOC = FALSE, width = , height = 5)

########################
## In this workflow, separate WT and KO per cluster and add a new metadata column called Clusters3
# create a new object so you don't overwrite the other workflow

table(Tbx1_all_aggr$Clusters2,Tbx1_all_aggr$Sample)

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters3
condition.vec.01 <- substr(Tbx1_all_aggr$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_all_aggr$Clusters2)
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

Tbx1_all_aggr$Clusters3 <- as.character(condition.vec)
table(Tbx1_all_aggr$Clusters3)

### IMPORTANT, save as a new project. Then re-load the new project
saveArchRProject(ArchRProj = Tbx1_all_aggr, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_all_aggr", load = TRUE)



########################
########################
# Proceed to step 4


