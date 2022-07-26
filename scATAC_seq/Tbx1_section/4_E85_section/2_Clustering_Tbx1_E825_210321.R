###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
##### Chapters 4 - 6, Clustering without batch correction

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

### 4. Iterative LSI (TF-IDF and SVD)
Tbx1_E825 <- addIterativeLSI(
  ArchRProj = Tbx1_E825,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), maxClusters=NULL), 
  varFeatures = 25000, 
  dimsToUse = 2:30,
  LSIMethod = 2,
  seed = 1,
  force = T
)

# Batch correction with Harmony --> from now on, use "Harmony" for reducedDims and "UMAP" for embedding
Tbx1_E825 <- addHarmony(
  ArchRProj = Tbx1_E825,
  dimsToUse = 2:30,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample", force = T
)

#### 5. Clustering
Tbx1_E825 <- addClusters(
  input = Tbx1_E825,
  reducedDims = "Harmony",
  dimsToUse = 2:30,
  method = "Seurat",
  name = "Clusters",
  maxClusters = NULL,
  resolution = 0.6, force = T
)

# ## Confusion matrix of sample number vs cluster 
cM <- confusionMatrix(paste0(Tbx1_E825$Clusters), paste0(Tbx1_E825$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p
plotPDF(p, name = "Plot-ConfusionMatrix-Heatmap.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)

#### 6. Single-cell Embeddings
### UMAP with IterativeLSI
Tbx1_E825 <- addUMAP(
  ArchRProj = Tbx1_E825, 
  reducedDims = "Harmony", 
  dimsToUse = 2:30,
  name = "UMAPHarmony",
  metric = "cosine",
  force = T
)

p9 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p10 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p9, p10, type = "h")
plotPDF(p9,p10, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)

table(Tbx1_E825$Clusters,Tbx1_E825$Sample)
# a_KO_1 b_WT_1 c_KO_2 d_KO_3 e_WT_2 f_WT_3 g_WT_4
# C1      39     58     46     12      9     43     83
# C10     18     42     48     20     10     24     14
# C11      7     34     38     18      4      9      2
# C12    140    207    517    174     15    145     40
# C13     22     25    196     41      1     10     14
# C14    132    340    238    156    181    157    237
# C15    144    292    300    226     78    167    177
# C16     37     78    138     46      9     47    103
# C17    246    471    556    387     55    237     17
# C18    132    322    300    240     81    172    118
# C19     60    174    123     75     55    103     75
# C2      66    114    132    100     61     84     78
# C3      24     31     69     25      8     18     31
# C4      49    166    369    206     45    118     83
# C5     177    288    312    133    142    206    124
# C6      51     50     92     61     23     47      0
# C7       8     20     21      9      8     14      4
# C8     126    304    271    177     86    170     23
# C9     112    171    326    125     41    117     55

saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825", load = TRUE)


####################################################
## Proceed to 3. scRNA-seq integration and cluster ID



