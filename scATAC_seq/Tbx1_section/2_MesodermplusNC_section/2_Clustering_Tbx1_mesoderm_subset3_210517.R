###### Tbx1 WT v KO mesoderm subset1
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
##### Chapters 4 - 6, Clustering with batch correction

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

### 4. Iterative LSI (TF-IDF and SVD)
mesoderm_subset <- addIterativeLSI(
  ArchRProj = mesoderm_subset,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2),
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
mesoderm_subset <- addHarmony(
  ArchRProj = mesoderm_subset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  dimsToUse = 2:20,
  groupBy = "Sample", force = T
)

#### 5. Clustering
mesoderm_subset <- addClusters(
  input = mesoderm_subset,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  maxClusters = 40,
  dimsToUse = 2:20,
  resolution = 0.5, force = T
)

# ## Confusion matrix of sample number vs cluster 
cM <- confusionMatrix(paste0(mesoderm_subset$Clusters), paste0(mesoderm_subset$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p
plotPDF(p, name = "Plot-ConfusionMatrix-Heatmap.pdf", ArchRProj = mesoderm_subset, addDOC = FALSE, width = 5, height = 5)

#### 6. Single-cell Embeddings
### UMAP with Harmony
mesoderm_subset <- addUMAP(
  ArchRProj = mesoderm_subset, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony",
  metric = "cosine",
  force = T
)

p9 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p10 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p9, p10, type = "h")
plotPDF(p9,p10, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = mesoderm_subset, addDOC = FALSE, width = 5, height = 5)

table(mesoderm_subset$Clusters,mesoderm_subset$Sample)
#         a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 j_E115_WT_1
# C1          582         351         276         406         161          79          73          65         317         124
# C10         254         129          57         240         711         942         930         837         189         928
# C11           0           0           0           0           0           0           0           0          14           1
# C12           0           0           0           2          92         239         199         117          67         131
# C13           1           0           0           1          48          14          57          63          83         262
# C14          35          23          39          32          37          32          43          31          49          50
# C15         376         239         212         326         347         530         431         320         346         277
# C16         270         199         318         245         214         340         337         286         701         148
# C17           2           1           3           2           4           7           7           8           7         143
# C18          39          23          95          59          30         145         112         117         117         125
# C19         116          74          39         136          69          54         128         122         100         192
# C2            4           1           1           0          26          27          23          23          17          45
# C3          230         160         238         284         226         209         280         220         568         328
# C4          106          79         134         140         127         233         145         161         391         262
# C5          148          81          98         115         311         385         381         460         354         578
# C6          384         262         336         410         589         326         291         236         877         455
# C7            1           5           5           9           6           0           4           4          14         444
# C8          157         112         176         270         487          69         257         139         164          49
# C9           38          39          41          67         141          93         108         120         153          74
# 
#         k_E115_WT_2 l_E115_KO_2 m_E115_KO_3
# C1          166         172         457
# C10        1227          85         163
# C11           0           5         134
# C12         282          26         102
# C13         308          57         131
# C14          72          29          76
# C15         560         267         613
# C16         420         484        1331
# C17          26           5          17
# C18         187         129         144
# C19         263          70         155
# C2           78           8          24
# C3          557         284        1007
# C4          489         209         669
# C5          933         258         414
# C6          915         690        1500
# C7           22           7          21
# C8           44          75         379
# C9          155          69         201

saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset1_210517", load = TRUE)


####################################################
## Proceed to 3. scRNA-seq integration and cluster ID



