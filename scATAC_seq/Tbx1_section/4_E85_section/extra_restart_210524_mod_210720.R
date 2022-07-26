###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
##### Re-start project on 05/24

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


######## 
# Re-load ArchR project!
Tbx1_E825 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_E825_only/Tbx1_E825")
Tbx1_E825
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825 
# samples(7): c_KO_2 b_WT_1 ... f_WT_3 e_WT_2
# sampleColData names(1): ArrowFiles
# cellColData names(23): Sample TSSEnrichment ... ReadsInPeaks FRIP
# numberOfCells(1): 15178
# medianTSS(1): 11.38
# medianFrags(1): 51252
######## 

## Table of Clusters
# Clusters = simple clustering all cells
table(Tbx1_E825$Clusters)
# C1  C10  C11  C12  C13  C14  C15  C16  C17  C18  C19   C2  C20  C21  C22  C23  C24  C25   C3   C4   C5   C6   C7   C8   C9 
# 625 1333  131  296  103  295  234 1309  444  314  150  289 1038  660 1220  153  699  840 1162  714  205 1103  735  559  567

# Clusters2 = Clusters1 separated by WT and KO (this got over-written at some point, it used to be just the 
## manually annotated clusters)
table(Tbx1_E825$Clusters2)
# C1_K  C1_W C10_K C10_W C11_K C11_W C12_K C12_W C13_K C13_W C14_K C14_W C15_K C15_W C16_K C16_W C17_K C17_W C18_K C18_W C19_K 
# 292   333   471   862    56    75   157   139    47    56   253    42   111   123   654   655   216   228   202   112   112 
# C19_W  C2_K  C2_W C20_K C20_W C21_K C21_W C22_K C22_W C23_K C23_W C24_K C24_W C25_K C25_W  C3_K  C3_W  C4_K  C4_W  C5_K  C5_W 
# 38    96   193   611   427   421   239   823   397    50   103   274   425   393   447   581   581   453   261   116    89 
# C6_K  C6_W  C7_K  C7_W  C8_K  C8_W  C9_K  C9_W 
# 678   425   323   412   240   319   283   284 

# Clusters3 = manually annotated clusters, divided by WT v KO
table(Tbx1_E825$Clusters3)
# AHF_K                AHF_W          CM_Atrium_K          CM_Atrium_W             CM_OFT_K 
# 182                  297                  368                  402                  127 
# CM_OFT_W       CM_Ventricle_K       CM_Ventricle_W                CNC_K                CNC_W 
# 355                  397                  758                  259                   40 
# Ectoderm_K           Ectoderm_W        Endocardium_K        Endocardium_W           Endoderm_K 
# 2289                 1366                  300                  339                 1493 
# Endoderm_W         Epicardium_K         Epicardium_W                LPM_K                LPM_W 
# 1485                  514                  612                  216                  237 
# PharyngealMesoderm_K PharyngealMesoderm_W                 PM_K                 PM_W               pSHF_K 
# 496                  275                  630                  595                  642 
# pSHF_W 
# 504


########################################################################################################################
########################################################################################################################
# What I did was re-run the entire workflow from script #2 up to #5 (DAR). I tried different versions of 
# clustering and liked the harmony + dims 2:30. The peak calls right now are for harmony etc for Clusters 4
# overall, the conclusions are that the mesoderm progenitors are fucked, the myocytes are spared and one 
# population of endoderm cells (that have GAS for Tbx1) do show slight defects...
# also, i am really liking the whole aggr object + subset'ing. In the end, that may be the winner here!
# save project and move to hub on 03/22.
########################################################################################################################
########################################################################################################################
# re-start this on 05/16 and re-visit conclusions from E8.25 stage. Key question, include in paper? We have no 
# scRNA-seq for this time point. but maybe its ok. maybe the key conclusion is that mesoderm subsets are fucked
# but that does not mean any given cell type does not form...(unlike hand2 e.g.)

### as mentioned before, key may be to zoom in on mesoderm subset and be solid there. Then bring in neural crest at the very end
########################################################################################################################
########################################################################################################################
## re-re-vist on 05/24 after re-annotated E8.5 clusters
########################################################################################################################


# 1. re-visit cluster annotation.
p1 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p1, p2, type = "h")
# 19 clusters, looks nicely resolved, but re-visit manual annotation!
#################
#### Visualize GAS
Tbx1_E825 <- addImputeWeights(Tbx1_E825, reducedDims = "Harmony", dimsToUse = 2:20)

markerGenes <- c("Tnnt2","Epcam","Sox2","Isl1","Osr1","Fgf8","Tbx1","Pecam1","Wt1","Dlx2","Hba-x","Mab21l2","Meox1",
                 "Tdgf1","Rgs5","Fgf10","Irx4","Vsnl1","Cited1","Nr2f1","Foxd1","Foxc2","Alx1","Tlx1","Aplnr","Nrg1","Hand1","Tbx18","Foxf1","Bmp2","Rspo3")

p13 <- plotEmbedding(
  ArchRProj = Tbx1_E825, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(Tbx1_E825)
)
plotPDF(plotList = p13, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = Tbx1_E825, 
        addDOC = FALSE, width = 5, height = 5)

######## this is under-clustered...
########
# recluster
########

### 4. Iterative LSI (TF-IDF and SVD)
Tbx1_E825 <- addIterativeLSI(
  ArchRProj = Tbx1_E825,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.3), maxClusters=NULL), 
  varFeatures = 25000, 
  dimsToUse = 2:25,
  LSIMethod = 2,
  seed = 1,
  force = T
)

# Batch correction with Harmony --> from now on, use "Harmony" for reducedDims and "UMAP" for embedding
Tbx1_E825 <- addHarmony(
  ArchRProj = Tbx1_E825,
  dimsToUse = 2:25,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample", force = T
)

#### 5. Clustering
Tbx1_E825 <- addClusters(
  input = Tbx1_E825,
  reducedDims = "Harmony",
  dimsToUse = 2:25,
  method = "Seurat",
  name = "Clusters",
  maxClusters = NULL,
  resolution = 1.4, force = T
)

### 
Tbx1_E825 <- addUMAP(
  ArchRProj = Tbx1_E825, 
  reducedDims = "Harmony", 
  dimsToUse = 2:25,
  name = "UMAPHarmony",
  metric = "cosine",
  force = T
)

plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

### re-do GAS
###### Browser Tracks
p_browser <- plotBrowserTrack(
  ArchRProj = Tbx1_E825, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 10000,
  downstream = 10000
)
dev.off()
#grid::grid.newpage()
#grid::grid.draw(p_browser$Tbx20)
plotPDF(plotList = p_browser, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = Tbx1_E825, 
        addDOC = FALSE, width = 5, height = 5)
### re-integrate with scRNAseq






p9 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p10 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
#ggAlignPlots(p9, p10, type = "h")
plotPDF(p9,p10, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)





#################
library(Seurat)
scRNA <- readRDS("/Users/sranade/scRNA-seq/2019_scRNA-seq/WT_B6/E85_aggr_201022/E85_aggr_mnn_annotated_210524.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)
# Ectoderm        VentricleCM   ParaxialMesoderm           AtrialCM           Endoderm PharyngealMesoderm 
# 1792                700               1041                513               1403                435 
# OFT_CM                JCF        Endocardium               pSHF                AHF                CNC 
# 422                385                379                343                258                188 
# AVC_CM 
# 160  

seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)
# spot check that its the same as line 93

#### 8.1.1 Unconstrained Integration, no imputation. 
Tbx1_E825 <- addGeneIntegrationMatrix(
  ArchRProj = Tbx1_E825, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupRNA = "clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  nGenes = 3000,
  useImputation = F,
  force = T
)

pal <- paletteDiscrete(values = seRNA$clusters)
p15 <- plotEmbedding(
  Tbx1_E825, 
  colorBy = "cellColData",
  name = "predictedGroup",
  embedding = "UMAPHarmony",
  pal = pal)
p15
plotPDF(p15, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)

##########################################
# call peaks on clusters1, wt and ko

table(Tbx1_E825$Clusters,Tbx1_E825$Sample)

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters3
condition.vec.01 <- substr(Tbx1_E825$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_E825$Clusters)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
# condition.vec
# C1_K  C1_W C10_K C10_W C11_K C11_W C12_K C12_W C13_K C13_W C14_K C14_W C15_K C15_W C16_K C16_W C17_K C17_W C18_K C18_W C19_K 
# 292   333   471   862    56    75   157   139    47    56   253    42   111   123   654   655   216   228   202   112   112 
# C19_W  C2_K  C2_W C20_K C20_W C21_K C21_W C22_K C22_W C23_K C23_W C24_K C24_W C25_K C25_W  C3_K  C3_W  C4_K  C4_W  C5_K  C5_W 
# 38    96   193   611   427   421   239   823   397    50   103   274   425   393   447   581   581   453   261   116    89 
# C6_K  C6_W  C7_K  C7_W  C8_K  C8_W  C9_K  C9_W 
# 678   425   323   412   240   319   283   284  

Tbx1_E825$Clusters2 <- as.character(condition.vec)
table(Tbx1_E825$Clusters2)


## 9. Pseudo-bulk Replicates in ArchR
Tbx1_E825 <- addGroupCoverages(ArchRProj = Tbx1_E825, groupBy = "Clusters2", minCells = 40, maxCells = 1500, force = T)
#2021-05-24 18:57:51 : Finished Creation of Coverage Files!, 45.376 mins elapsed.

## 10. Peak calling
Tbx1_E825 <- addReproduciblePeakSet(
  ArchRProj = Tbx1_E825, 
  groupBy = "Clusters2", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 1500,
  maxPeaks = 125000,
  minCells = 25,
  force = T
)

# 2020-09-28 21:59:40 : Finished Creating Union Peak Set (416062)!, 77.334 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(Tbx1_E825)
merged_peak_setGR <- getPeakSet(Tbx1_E825)

# add peakmatrix to object
Tbx1_E825 <- addPeakMatrix(Tbx1_E825)
getAvailableMatrices(Tbx1_E825)

########################################################
# DAR per cluster
###############################################
# II. Script part 1, running DARs on each cluster in Clusters2

table(mesoderm_subset$Clusters2)

### II. Script part 2, running DARs on each cluster in manually annotated Clusters3
## For all clusters, use FDR <= 0.05 ad log2FC > 1. This is prob the safe for FDR and still allows you to catch as many 
## peaks as possible since its not entirely clear what "1.8x more accessible actually means"

#
markerTest_C1 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1_W",
  bgdGroups = "C1_K"
)
pma_C1 <- plotMarkers(seMarker = markerTest_C1, name = "C1_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C2 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C2_W",
  bgdGroups = "C2_K"
)
pma_C2 <- plotMarkers(seMarker = markerTest_C2, name = "C2_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C3 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C3_W",
  bgdGroups = "C3_K"
)
pma_C3 <- plotMarkers(seMarker = markerTest_C3, name = "C3_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C4 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C4_W",
  bgdGroups = "C4_K"
)
pma_C4 <- plotMarkers(seMarker = markerTest_C4, name = "C4_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C5 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C5_W",
  bgdGroups = "C5_K"
)
pma_C5 <- plotMarkers(seMarker = markerTest_C5, name = "C5_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C6 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C6_W",
  bgdGroups = "C6_K"
)
pma_C6 <- plotMarkers(seMarker = markerTest_C6, name = "C6_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C7 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C7_W",
  bgdGroups = "C7_K"
)
pma_C7 <- plotMarkers(seMarker = markerTest_C7, name = "C7_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C8 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C8_W",
  bgdGroups = "C8_K"
)
pma_C8 <- plotMarkers(seMarker = markerTest_C8, name = "C8_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C9 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C9_W",
  bgdGroups = "C9_K"
)
pma_C9 <- plotMarkers(seMarker = markerTest_C9, name = "C9_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C10 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C10_W",
  bgdGroups = "C10_K"
)
pma_C10 <- plotMarkers(seMarker = markerTest_C10, name = "C10_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C11 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C11_W",
  bgdGroups = "C11_K"
)
pma_C11 <- plotMarkers(seMarker = markerTest_C11, name = "C11_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C12 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C12_W",
  bgdGroups = "C12_K"
)
pma_C12 <- plotMarkers(seMarker = markerTest_C12, name = "C12_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C13 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C13_W",
  bgdGroups = "C13_K"
)
pma_C13 <- plotMarkers(seMarker = markerTest_C13, name = "C13_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C14 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C14_W",
  bgdGroups = "C14_K"
)
pma_C14 <- plotMarkers(seMarker = markerTest_C14, name = "C14_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C15 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C15_W",
  bgdGroups = "C15_K"
)
pma_C15 <- plotMarkers(seMarker = markerTest_C15, name = "C15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C16 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C16_W",
  bgdGroups = "C16_K"
)
pma_C16 <- plotMarkers(seMarker = markerTest_C16, name = "C16_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C17 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C17_W",
  bgdGroups = "C17_K"
)
pma_C17 <- plotMarkers(seMarker = markerTest_C17, name = "C17_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C18 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C18_W",
  bgdGroups = "C18_K"
)
pma_C18 <- plotMarkers(seMarker = markerTest_C18, name = "C18_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C19 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C19_W",
  bgdGroups = "C19_K"
)
pma_C19 <- plotMarkers(seMarker = markerTest_C19, name = "C19_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C20 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C20_W",
  bgdGroups = "C20_K"
)
pma_C20 <- plotMarkers(seMarker = markerTest_C20, name = "C20_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C21 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C21_W",
  bgdGroups = "C21_K"
)
pma_C21 <- plotMarkers(seMarker = markerTest_C21, name = "C21_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C22 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C22_W",
  bgdGroups = "C22_K"
)
pma_C22 <- plotMarkers(seMarker = markerTest_C22, name = "C22_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C23 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C23_W",
  bgdGroups = "C23_K"
)
pma_C23 <- plotMarkers(seMarker = markerTest_C23, name = "C23_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C24 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C24_W",
  bgdGroups = "C24_K"
)
pma_C24 <- plotMarkers(seMarker = markerTest_C24, name = "C24_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_C25 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C25_W",
  bgdGroups = "C25_K"
)
pma_C25 <- plotMarkers(seMarker = markerTest_C25, name = "C25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#####
setwd("/Users/sranade/scATAC-seq/Tbx1_E825/Exports")
# 
markerList_C6 <- getMarkers(markerTest_C6, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_C6_df <- as.data.frame(markerList_C6@listData[["C6_W"]])
markerList_C6_df_homer <- markerList_C6_df[order(markerList_C6_df$Log2FC, decreasing = TRUE),]
markerList_C6_df_homer$rank <- seq_len(nrow(markerList_C6_df_homer))
markerList_C6_df_homer$blank <- ""
markerList_C6_df_homer_trim <- markerList_C6_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C6_df_homer_trim, file = "markerList_C6_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

##############
# save object
##############
saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E825", load = TRUE)

############################
############################
# subset scATACseq object to include only mesoderm population

# Subset clusters of interest in the hackiest fucking way 
# clone object
mesoderm_subset <- Tbx1_E825

# assign Y/N values to only the clusters you want
library(plyr)
temp_cluster_names <- as.character(revalue(mesoderm_subset$Clusters, c(
  "C1"="False",
  "C2"="False",
  "C3"="True",
  "C4"="True",
  "C5"="False",
  "C6"="False",
  "C7"="False",
  "C8"="False",
  "C9"="True",
  "C10"="True",
  "C11"="True",
  "C12"="True",
  "C13"="False",
  "C14"="False",
  "C15"="True",
  "C16"="True",
  "C17"="True",
  "C18"="False",
  "C19"="False",
  "C20"="False",
  "C21"="False",
  "C22"="False",
  "C23"="False",
  "C24"="True",
  "C25"="True"
)))
length(unique(temp_cluster_names))
# Clusters3 = t/f. Clusters2 = WT/KO
mesoderm_subset$Clusters3 <- temp_cluster_names
mesoderm_subset@cellColData@listData[["Clusters3"]] <- temp_cluster_names
table(mesoderm_subset$Clusters3,mesoderm_subset$Sample)
# a_KO_1 b_WT_1 c_KO_2 d_KO_3 e_WT_2 f_WT_3 g_WT_4
# False    819   1428   2301   1144    365    920    472
# True     771   1759   1791   1087    547    968    806

## subset
idxSample <- BiocGenerics::which(mesoderm_subset$Clusters3 %in% "True")
cellsSample <- mesoderm_subset$cellNames[idxSample]
subsetArchRProject(
  ArchRProj = mesoderm_subset,
  cells = cellsSample,
  outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210524", force = T
)
## Dropping ImputeWeights Since You Are Subsetting Cells! ImputeWeights is a cell-x-cell Matrix!
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210516 
# samples(7): c_KO_2 b_WT_1 ... a_KO_1 e_WT_2
# sampleColData names(1): ArrowFiles
# cellColData names(25): Sample TSSEnrichment ... Clusters4 Clusters5
# numberOfCells(1): 7417
# medianTSS(1): 11.359
# medianFrags(1): 54455


##############
# 4. Make new project, new scripts and go back essentially to step 2
##############






