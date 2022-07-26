###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
#### Chapters 7 & 8, scRNA-seq integration and cluster identification

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 0.9.5


######## 
# Re-load ArchR project!
Tbx1_E825 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_E825")
Tbx1_E825
######## 

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
  embedding = "UMAP",
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
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = Tbx1_E825, 
        addDOC = FALSE, width = 5, height = 5)

######## 7.5 Marker Genes Imputation with MAGIC (only for GAS not GIS)
Tbx1_E825 <- addImputeWeights(Tbx1_E825, reducedDims = "IterativeLSI")

p13 <- plotEmbedding(
  ArchRProj = Tbx1_E825, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
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
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
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
  embedding = "UMAP",
  pal = pal)
p15
plotPDF(p15, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825", load = TRUE)
Tbx1_E825

########################
########################
########################

### go to temp intermediate script for manual spot-checking and reclustering, otherwise continue to step 4 after saving

########################
########################
########################


