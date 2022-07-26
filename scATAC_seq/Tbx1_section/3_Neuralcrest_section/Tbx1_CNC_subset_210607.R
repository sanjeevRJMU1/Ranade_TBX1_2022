###### Tbx1 WT v KO CNC subset
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
CNC_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/CNC_subset_210607")
CNC_subset
######## 
# numberOfCells(1): 12263
# medianTSS(1): 12.054
# medianFrags(1): 54893
########

### 4. Iterative LSI (TF-IDF and SVD)
CNC_subset <- addIterativeLSI(
  ArchRProj = CNC_subset,
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
  dimsToUse = 1:30,
  LSIMethod = 2,
  seed = 1,
  force = T
)

# # # # # Batch correction with Harmony --> from now on, use "Harmony" for reducedDims and "UMAP" for embedding
 CNC_subset <- addHarmony(
   ArchRProj = CNC_subset,
   reducedDims = "IterativeLSI",
   name = "Harmony",
   dimsToUse = 1:30,
   groupBy = "Sample", force = T
 )

#### 5. Clustering
CNC_subset <- addClusters(
  input = CNC_subset,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  maxClusters = 40,
  dimsToUse = 1:30,
  resolution = 0.3, force = T
)

#### 6. Single-cell Embeddings
### UMAP with Harmony
CNC_subset <- addUMAP(
  ArchRProj = CNC_subset, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony",
  metric = "cosine",
  dimsToUse = 1:30,
  force = T
)

p9 <- plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p10 <- plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p9, p10, type = "h")
plotPDF(p9,p10, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = CNC_subset, addDOC = FALSE, width = 5, height = 5)

table(CNC_subset$Clusters,CNC_subset$Sample)
table(CNC_subset$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 j_E115_WT_1
# C1            0           0           0           0         181         277         198         170          40         113
# C10           1           1           1           0          10           6          37          17          36         158
# C11           0           3           0           0          39          12          23          44          62         138
# C2            0           0           0           0         162         136         138         156          87         168
# C3            0           0           0           0          84          64          58          82          13          53
# C4           82          36          25          83          73          70          74          85           0          90
# C5            0           0           0           0          91         106         201         161          21         212
# C6            0           0           0           0         130         218         211         186          12         222
# C7          170          91          28         155          33          92          49          20           7          21
# C8            0           0           0           0          70         201         187         101          47         135
# C9            1           0           0           1           2          14           8           7          25          16
# 
# k_E115_WT_2 l_E115_KO_2 m_E115_KO_3
# C1          136           5          35
# C10         197          30         105
# C11         135          50         127
# C2          192          11          42
# C3           85          13          17
# C4          109           0           4
# C5          293           0           2
# C6          358          16          12
# C7           58           9          32
# C8          228          11          44
# C9           14          32          51

# GAS
CNC_subset <- addImputeWeights(CNC_subset, reducedDims = "Harmony",dimsToUse = 2:20)
markerGenes<-c("Emx2","Smoc1","Ebf1","Dlx1","Dlx2","Dlx3","Dlx5","Dlx6",
               "Hand2","Hand1","Dlk1","Msx2","Foxf1","Gata3","Rgs5","Isl1",
               "Pou3f3","Barx1","Sox2","Sox10","Cdh19",
               "Hoxb3","Hoxd4","Six2","Hoxa3","Tbx2","Meox1","Hoxa4","Mef2c","Hoxb4")
p13 <- plotEmbedding(
  ArchRProj = CNC_subset, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(CNC_subset)
)
plotPDF(plotList = p13, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = CNC_subset, 
        addDOC = FALSE, width = 5, height = 5)
#################

markerGenes_spotmarkerGenes_spot <- c("Tnnt2","Isl1","Osr1","Fgf8","Tbx1","Wt1","Mab21l2","Meox1","Hand2",
                                      "Tdgf1","Rgs5","Fgf10","Irx4","Vsnl1","Cited1","Nr2f1","Foxd1","Foxc2",
                                      "Alx1","Tlx1","Aplnr","Nrg1","Hand1","Tbx18","Foxf1","Bmp2","Rspo3",
                                      "Dlx1","Dlx2","Dlx3","Dlx4","Dlx5","Dlx6",
                                      "Tdgf1","Irx3","Irx1","Tbx2","Tbx3","Nkx2-5","Gata4",
                                      "Lix1","Tcf21","Lhx1","Lhx2","Gbx1","Otx2","Sox10","Sox2","Myf5","Pax3","Pitx2",
                                      "Twist1","Prrx1","Prrx2","Col1a1","Pdgfra","Ebf1","Ebf2","Lum","Sox9","Msx1","Msx2",
                                      "Pax1","Pax9","Pax5","Wnt4","Hoxa4","Hoxb4","Hoxc4","Hoxd4","Bdnf","Asb4","Zic1","Foxg1","Foxp2",
                                      "Tagln","Acta2","Des","Vcl")

p14 <- plotEmbedding(
  ArchRProj = CNC_subset, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes_spotmarkerGenes_spot, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(CNC_subset)
)
plotPDF(plotList = p14, 
        name = "Plot-UMAP-Marker-Genes-SpotCheck.pdf", 
        ArchRProj = CNC_subset, 
        addDOC = FALSE, width = 5, height = 5)

####################################################
## This now uses the scRNAseq from the same tissue used for scATACseq
library(Seurat)
scRNA <- readRDS("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/Tbx1_CNC_annotated.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)
# Craniofacial_NeuralCrest      Cardiac_NeuralCrest    Migratory_NeuralCrest  PA3_Cardiac_NeuralCrest 
# 4226                     3405                     1460                     1140 

seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)
# spot check that its the same as line 93

#### 8.1.1 Unconstrained Integration, no imputation.
CNC_subset <- addGeneIntegrationMatrix(
  ArchRProj = CNC_subset,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupRNA = "clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  useImputation = F,
  force = T
)

pal <- paletteDiscrete(values = seRNA$clusters)
p15 <- plotEmbedding(
  CNC_subset,
  colorBy = "cellColData",
  name = "predictedGroup",
  embedding = "UMAPHarmony",
  pal = pal)
p15
plotPDF(p15, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = CNC_subset, addDOC = FALSE, width = 5, height = 5)
####################################################
## My manual labeling --> this goes to Clusters6 as rest are already taken
library(plyr)
temp_cluster_names <- as.character(revalue(CNC_subset$Clusters, c(
  "C1"="Cardiac_NeuralCrest",
  "C2"="PA3_Cardiac_NeuralCrest",
  "C3"="Cardiac_NeuralCrest",
  "C4"="Cardiac_NeuralCrest",
  "C5"="PA3_Cardiac_NeuralCrest",
  "C6"="Craniofacial_NeuralCrest",
  "C7"="Craniofacial_NeuralCrest",
  "C8"="Migratory_NeuralCrest",
  "C9"="Migratory_NeuralCrest",
  "C10"="Craniofacial_NeuralCrest",
  "C11"="Craniofacial_NeuralCrest"
)))
length(unique(temp_cluster_names))
CNC_subset$Clusters6 <- temp_cluster_names
CNC_subset@cellColData@listData[["Clusters6"]] <- temp_cluster_names
table(CNC_subset$Clusters6,CNC_subset$Sample)

################################################################

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called Clusters7
condition.vec.01 <- substr(CNC_subset$Sample,5,8)
unique(condition.vec.01)
condition.vec.02 <- as.vector(CNC_subset$Clusters6)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)

CNC_subset$Clusters7 <- as.character(condition.vec)
# spot check
table(CNC_subset$Clusters7)

p20 <- plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Clusters6", embedding = "UMAPHarmony")


p20

plotPDF(p20,
        name = "Plot-UMAP-Manual-Annotations.pdf",
        ArchRProj = CNC_subset,
        addDOC = FALSE, width = , height = 5)
saveArchRProject(ArchRProj = CNC_subset, outputDirectory = "/Users/sranade/scATAC-seq/CNC_subset_210607", load = TRUE)
####################################################
### Peak calling
## 9. Pseudo-bulk Replicates in ArchR
CNC_subset <- addGroupCoverages(ArchRProj = CNC_subset, groupBy = "Clusters7", minCells = 40, maxCells = 2000, force = T, maxFragments = 75 * 10^6)
#2021-06-08 06:42:10 : Finished Creation of Coverage Files!, 17.396 mins elapsed.

## 10. Peak calling
CNC_subset <- addReproduciblePeakSet(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters7", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 5000,
  maxPeaks = 150000,
  minCells = 10,
  force = T
)

#2021-06-08 07:16:36 : Creating Union Peak Set!, 34.418 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(CNC_subset)
merged_peak_setGR <- getPeakSet(CNC_subset)

# add peakmatrix to object
CNC_subset <- addPeakMatrix(CNC_subset)
getAvailableMatrices(CNC_subset)

#saveArchRProject(ArchRProj = CNC_subset, outputDirectory = "/Users/sranade/scATAC-seq/CNC_subset_210607", load = TRUE)
#####################
###############################################
# II. Script part 1, running DARs on each cluster in Clusters6

table(CNC_subset$Clusters7)

# #
# markerTest_Craniofacial_NeuralCrest_E925 <- getMarkerFeatures(
#   ArchRProj = CNC_subset, 
#   useMatrix = "PeakMatrix",
#   groupBy = "Clusters7",
#   binarize = F,
#   testMethod = "wilcoxon",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   useGroups = "Craniofacial_NeuralCrest_25_W",
#   bgdGroups = "Craniofacial_NeuralCrest_25_K"
# )
# pma_Craniofacial_NeuralCrest_E925 <- plotMarkers(seMarker = markerTest_Craniofacial_NeuralCrest_E925, name = "Craniofacial_NeuralCrest_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
# 
# markerTest_Craniofacial_NeuralCrest_E105 <- getMarkerFeatures(
#   ArchRProj = CNC_subset, 
#   useMatrix = "PeakMatrix",
#   groupBy = "Clusters7",
#   binarize = F,
#   testMethod = "wilcoxon",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   useGroups = "Craniofacial_NeuralCrest_05_W",
#   bgdGroups = "Craniofacial_NeuralCrest_05_K"
# )
# pma_Craniofacial_NeuralCrest_E105 <- plotMarkers(seMarker = markerTest_Craniofacial_NeuralCrest_E105, name = "Craniofacial_NeuralCrest_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Craniofacial_NeuralCrest_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Craniofacial_NeuralCrest_15_W",
  bgdGroups = "Craniofacial_NeuralCrest_15_K"
)
pma_Craniofacial_NeuralCrest_E115 <- plotMarkers(seMarker = markerTest_Craniofacial_NeuralCrest_E115, name = "Craniofacial_NeuralCrest_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
pma_Craniofacial_NeuralCrest_E115
#
# #
# markerTest_Cardiac_NeuralCrest_E925 <- getMarkerFeatures(
#   ArchRProj = CNC_subset, 
#   useMatrix = "PeakMatrix",
#   groupBy = "Clusters7",
#   binarize = F,
#   testMethod = "wilcoxon",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   useGroups = "Cardiac_NeuralCrest_25_W",
#   bgdGroups = "Cardiac_NeuralCrest_25_K"
# )
# pma_Cardiac_NeuralCrest_E925 <- plotMarkers(seMarker = markerTest_Cardiac_NeuralCrest_E925, name = "Cardiac_NeuralCrest_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
# 
markerTest_Cardiac_NeuralCrest_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiac_NeuralCrest_05_W",
  bgdGroups = "Cardiac_NeuralCrest_05_K"
)
pma_Cardiac_NeuralCrest_E105 <- plotMarkers(seMarker = markerTest_Cardiac_NeuralCrest_E105, name = "Cardiac_NeuralCrest_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cardiac_NeuralCrest_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiac_NeuralCrest_15_W",
  bgdGroups = "Cardiac_NeuralCrest_15_K"
)
pma_Cardiac_NeuralCrest_E115 <- plotMarkers(seMarker = markerTest_Cardiac_NeuralCrest_E115, name = "Cardiac_NeuralCrest_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
pma_Cardiac_NeuralCrest_E115

#
markerTest_Migratory_NeuralCrest_E925 <- getMarkerFeatures(
  ArchRProj = CNC_subset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Migratory_NeuralCrest_25_W",
  bgdGroups = "Migratory_NeuralCrest_25_K"
)
pma_Migratory_NeuralCrest_E925 <- plotMarkers(seMarker = markerTest_Migratory_NeuralCrest_E925, name = "Migratory_NeuralCrest_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Migratory_NeuralCrest_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Migratory_NeuralCrest_05_W",
  bgdGroups = "Migratory_NeuralCrest_05_K"
)
pma_Migratory_NeuralCrest_E105 <- plotMarkers(seMarker = markerTest_Migratory_NeuralCrest_E105, name = "Migratory_NeuralCrest_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Migratory_NeuralCrest_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Migratory_NeuralCrest_15_W",
  bgdGroups = "Migratory_NeuralCrest_15_K"
)
pma_Migratory_NeuralCrest_E115 <- plotMarkers(seMarker = markerTest_Migratory_NeuralCrest_E115, name = "Migratory_NeuralCrest_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

pma_Migratory_NeuralCrest_E925
pma_Migratory_NeuralCrest_E105
pma_Migratory_NeuralCrest_E115
# #
# markerTest_PA3_Cardiac_NeuralCrest_E925 <- getMarkerFeatures(
#   ArchRProj = CNC_subset, 
#   useMatrix = "PeakMatrix",
#   groupBy = "Clusters7",
#   binarize = F,
#   testMethod = "wilcoxon",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   useGroups = "PA3_Cardiac_NeuralCrest_25_W",
#   bgdGroups = "PA3_Cardiac_NeuralCrest_25_K"
# )
# pma_PA3_Cardiac_NeuralCrest_E925 <- plotMarkers(seMarker = markerTest_PA3_Cardiac_NeuralCrest_E925, name = "PA3_Cardiac_NeuralCrest_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_PA3_Cardiac_NeuralCrest_E105 <- getMarkerFeatures(
  ArchRProj = CNC_subset,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PA3_Cardiac_NeuralCrest_05_W",
  bgdGroups = "PA3_Cardiac_NeuralCrest_05_K"
)
pma_PA3_Cardiac_NeuralCrest_E105 <- plotMarkers(seMarker = markerTest_PA3_Cardiac_NeuralCrest_E105, name = "PA3_Cardiac_NeuralCrest_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_PA3_Cardiac_NeuralCrest_E115 <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "PA3_Cardiac_NeuralCrest_15_W",
  bgdGroups = "PA3_Cardiac_NeuralCrest_15_K"
)
pma_PA3_Cardiac_NeuralCrest_E115 <- plotMarkers(seMarker = markerTest_PA3_Cardiac_NeuralCrest_E115, name = "PA3_Cardiac_NeuralCrest_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
pma_PA3_Cardiac_NeuralCrest_E115




