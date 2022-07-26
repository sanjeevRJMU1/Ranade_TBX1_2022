##### Chapters 4 - 6, Clustering with batch correction

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 0.9.5


######## 
# Re-load ArchR project!
wt_atlas <- loadArchRProject(path = "/Users/sranade/scATAC-seq/new_wt_section_211111")
wt_atlas
######## 

### 4. Iterative LSI (TF-IDF and SVD)
wt_atlas <- addIterativeLSI(
  ArchRProj = wt_atlas,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1, 0.2, 0.4), 
    n.start = 10, maxClusters = NULL), 
  varFeatures = 25000, 
  dimsToUse = 2:25,
  LSIMethod = 2,
  seed = 1,
  force = T
)

# # Batch correction with Harmony --> from now on, use "Harmony" for reducedDims and "UMAP" for embedding
wt_atlas <- addHarmony(
  ArchRProj = wt_atlas,
  reducedDims = "IterativeLSI",
  dimsToUse = 2:25,
  name = "Harmony",
  groupBy = "Sample", force = T
)

#### 5. Clustering
wt_atlas <- addClusters(
  input = wt_atlas,
  reducedDims = "Harmony",
  dimsToUse = 2:25,
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,maxClusters = NULL, force = T
)

# ## Confusion matrix of sample number vs cluster 
cM <- confusionMatrix(paste0(wt_atlas$Clusters), paste0(wt_atlas$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p
plotPDF(p, name = "Plot-ConfusionMatrix-Heatmap.pdf", ArchRProj = wt_atlas, addDOC = FALSE, width = 5, height = 5)

#### 6. Single-cell Embeddings
### UMAP with harmony --> replace iterativeLSI with Harmony and UMAP with UMAP
wt_atlas <- addUMAP(
  ArchRProj = wt_atlas, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony",
  metric = "cosine",
  force = T
)

p9 <- plotEmbedding(ArchRProj = wt_atlas, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p10 <- plotEmbedding(ArchRProj = wt_atlas, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p9, p10, type = "h")
plotPDF(p9,p10, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = wt_atlas, addDOC = FALSE, width = 5, height = 5)

table(wt_atlas$Clusters,wt_atlas$Sample)


saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)


####################################################
## Proceed to 3. scRNA-seq integration and cluster ID



