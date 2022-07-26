###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
##### Chapters 4 - 6, Clustering without batch correction

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

######## 
# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210524")
mesoderm_subset
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210524 
# samples(7): c_KO_2 b_WT_1 ... a_KO_1 e_WT_2
# sampleColData names(1): ArrowFiles
# cellColData names(24): Sample TSSEnrichment ... FRIP Clusters4
# numberOfCells(1): 7729
# medianTSS(1): 11.436
# medianFrags(1): 53266
######## 

### 4. Iterative LSI (TF-IDF and SVD)
mesoderm_subset <- addIterativeLSI(
  ArchRProj = mesoderm_subset,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), maxClusters=NULL), 
  varFeatures = 25000, 
  dimsToUse = 2:25,
  LSIMethod = 2,
  seed = 1,
  force = T
)

# Batch correction with Harmony --> from now on, use "Harmony" for reducedDims and "UMAP" for embedding
mesoderm_subset <- addHarmony(
  ArchRProj = mesoderm_subset,
  dimsToUse = 2:25,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample", force = T
)

#### 5. Clustering
mesoderm_subset <- addClusters(
  input = mesoderm_subset,
  reducedDims = "Harmony",
  dimsToUse = 2:25,
  method = "Seurat",
  name = "Clusters",
  maxClusters = NULL,
  resolution = 0.8, force = T
)

#### 6. Single-cell Embeddings
### UMAP with IterativeLSI
mesoderm_subset <- addUMAP(
  ArchRProj = mesoderm_subset, 
  reducedDims = "Harmony", 
  dimsToUse = 2:25,
  name = "UMAPHarmony",
  metric = "cosine",
  force = T
)

p1 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = mesoderm_subset, addDOC = FALSE, width = 5, height = 5)

table(mesoderm_subset$Clusters,mesoderm_subset$Sample)
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

saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210524", load = TRUE)


####################################################
## Proceed to 3. scRNA-seq integration and cluster ID



