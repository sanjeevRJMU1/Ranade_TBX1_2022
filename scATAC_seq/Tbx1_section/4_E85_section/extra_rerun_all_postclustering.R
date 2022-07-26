

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 0.9.5


######## 
# Re-load ArchR project!
Tbx1_E825 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825")
Tbx1_E825




#### 7.3 Identifying Marker Genes
markersGS <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  binarize = F,
  maxCells = 1000
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1.25")

#### Visualize as heatmap
markerGenes <- c("Tnnt2","Epcam","Sox2","Isl1","Foxf1","Fgf8","Tbx1","Pecam1","Wt1","Shh","Dlx2","Hba-x")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
library(ComplexHeatmap)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = Tbx1_E825, addDOC = FALSE)

#### Visualize as featureplots
p11 <- plotEmbedding(
  ArchRProj = Tbx1_E825, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)


p12 <- lapply(p11, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p12))

plotPDF(plotList = p11, 
        name = "Plot-UMAPHarmony-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = Tbx1_E825, 
        addDOC = FALSE, width = 5, height = 5)

######## 7.5 Marker Genes Imputation with MAGIC (only for GAS not GIS)
Tbx1_E825 <- addImputeWeights(Tbx1_E825, reducedDims = "Harmony")

p13 <- plotEmbedding(
  ArchRProj = Tbx1_E825, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(Tbx1_E825)
)

#Rearrange for grid plotting
p14 <- lapply(p13, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p14))

plotPDF(plotList = p13, 
        name = "Plot-UMAPHarmony-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = Tbx1_E825, 
        addDOC = FALSE, width = 5, height = 5)


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

## Launch browser
#ArchRBrowser(Tbx1_E825)

###################################################################################################

################### Chapter 8 Defining Cluster Identity with scRNA-seq
## This now uses the scRNAseq from the same tissue used for scATACseq 
library(Seurat)
scRNA <- readRDS("/Users/sranade/scRNA-seq/2019_scRNA-seq/WT_B6/E85_aggr_201022/E85_aggr_mnn_annotated_201022.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)
# Ectoderm Pharyngeal_Mesoderm                pSHF                 CNC   Paraxial_Mesoderm            Endoderm 
# 2104                 580                1177                 901                 414                1732 
# Epicardium         Endocardium       Cardiomyocyte                 AHF     FHF_Progenitors               Blood 
# 375                 288                 529                 544                 209                 254 

seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)
# spot check that its the same as line 93

#### 8.1.1 Unconstrained Integration, no imputation. 
Tbx1_E825 <- addGeneIntegrationMatrix(
  ArchRProj = Tbx1_E825, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
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
plotPDF(p15, name = "Plot-UMAPHarmony-RNA-Integration.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)



# divide clusters by genotype
# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters4"
table(Tbx1_E825$Sample)
condition.vec.01 <- substr(Tbx1_E825$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_E825$Clusters)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
Tbx1_E825$Clusters4 <- as.character(condition.vec)

# Confirm clusters are separated by genotype
table(Tbx1_E825$Clusters4)

## 9. Pseudo-bulk Replicates in ArchR --> important!! Use Clusters4 to call peaks on WT and KO separately
Tbx1_E825 <- addGroupCoverages(ArchRProj = Tbx1_E825, groupBy = "Clusters4", minCells = 30, maxCells = 1500, force = T)
#2021-01-04 14:59:58 : Finished Creation of Coverage Files!, 28.109 mins elapsed.

## 10. Peak calling
Tbx1_E825 <- addReproduciblePeakSet(
  ArchRProj = Tbx1_E825, 
  groupBy = "Clusters4", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 2000,
  maxPeaks = 150000,
  minCells = 25,
  reproducibility = "2",
  force = T
)

# Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
# C1_K   C1_K     95         83           2   39   44   150000
# C1_W   C1_W    182        180           3   43   79   150000
# C2_K   C2_K    325        325           3   91  135   150000
# C2_W   C2_W    458        458           4   84  173   150000
# C3_K   C3_K    238        238           3   42  121   150000
# C3_W   C3_W    353        353           4   44  113   150000
# C4_K   C4_K    124        124           2   48   76   150000
# C4_W   C4_W     83         72           2   35   37   144000
# C5_K   C5_K    677        677           3   49  419   150000
# C5_W   C5_W    381        372           3   79  165   150000
# C6_K   C6_K    191        191           3   47   92   150000
# C6_W   C6_W    112         89           2   44   45   150000
# C7_K   C7_K    107         94           2   43   51   150000
# C7_W   C7_W     27         24           2   14   17    48000
# C8_K   C8_K     48         45           2   30   30    90000
# C8_W   C8_W     79         79           2   30   49   150000
# C9_K   C9_K    632        632           3  124  281   150000
# C9_W   C9_W    457        457           3   49  237   150000
# C10_K C10_K    413        413           3   99  206   150000
# C10_W C10_W    200        194           2   53  141   150000
# C11_K C11_K    831        831           3  139  518   150000
# C11_W C11_W    442        428           3   57  217   150000
# C12_K C12_K     69         69           2   30   39   138000
# C12_W C12_W     47         45           2   30   30    90000
# C13_K C13_K    564        564           3  122  263   150000
# C13_W C13_W    560        560           3   88  307   150000
# C14_K C14_K    555        555           3  112  327   150000
# C14_W C14_W     17         17           2   11   14    34000
# C15_K C15_K      1          1           2    1    1     2000
# C15_W C15_W    367        367           4   34  161   150000
# C16_K C16_K     82         82           2   36   46   150000
# C16_W C16_W     82         82           2   41   41   150000
# C17_K C17_K    256        235           2   40  195   150000
# C17_W C17_W     40         40           2   30   30    80000
# C18_K C18_K     38         38           2   30   30    76000
# C18_W C18_W     37         34           2   21   21    68000
# C19_K C19_K    299        299           3   66  135   150000
# C19_W C19_W    334        334           4   61  116   150000
# C20_K C20_K    341        341           3   87  130   150000
# C20_W C20_W    367        367           4   38  158   150000
# C21_K C21_K    418        418           3  102  203   150000
# C21_W C21_W    780        780           4  107  272   150000
# C22_K C22_K     97         97           3   31   35   150000
# C22_W C22_W    212        212           4   40   65   150000
# C23_K C23_K    663        663           3  108  314   150000
# C23_W C23_W    385        377           3   40  214   150000
# C24_K C24_K    222        222           3   51  104   150000
# C24_W C24_W    195        193           3   39  106   150000
# C25_K C25_K     75         75           1   75   75   150000
# C25_W C25_W    458        458           4   81  175   150000
# C26_K C26_K    230        230           3   36  148   150000
# C26_W C26_W    249        248           3   56  113   150000
# C27_K C27_K    322        322           3   63  150   150000
# C27_W C27_W    361        353           3   89  145   150000

# 2021-03-21 16:19:38 : Batching Peak Calls!, 0.008 mins elapsed.
# 2021-03-21 16:19:38 : Batch Execution w/ safelapply!, 0 mins elapsed.
# 2021-03-21 17:57:54 : Identifying Reproducible Peaks!, 98.266 mins elapsed.
# 2021-03-21 18:00:43 : Creating Union Peak Set!, 101.093 mins elapsed.
# Converged after 14 iterations!
#   Plotting Ggplot!
#   2021-03-21 18:01:26 : Finished Creating Union Peak Set (481439)!, 101.802 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(Tbx1_E825)
merged_peak_setGR <- getPeakSet(Tbx1_E825)

# add peakmatrix to object
Tbx1_E825 <- addPeakMatrix(Tbx1_E825, force = T, binarize = F, ceiling = 4)
getAvailableMatrices(Tbx1_E825)




#
markerTest_C1 <- getMarkerFeatures(
  ArchRProj = Tbx1_E825, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
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
  groupBy = "Clusters4",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C19_W",
  bgdGroups = "C19_K"
)
pma_C19 <- plotMarkers(seMarker = markerTest_C19, name = "C19_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

