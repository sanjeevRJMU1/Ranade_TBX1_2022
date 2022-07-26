###### Tbx1 WT v KO mesoderm subset1
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
#### Chapters 7 & 8, scRNA-seq integration and cluster identification

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset1_210517")
mesoderm_subset
######## 


#### 7.3 Identifying Marker Genes
markersGS <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  binarize = F,
  maxCells = 1000
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1.25")

#### Visualize as heatmap
markerGenes <- c("Tnnt2","Dlx5","Wt1","Isl1","Foxf1","Fgf8","Tbx1","Mab21l2","Meox1","Osr1","Fgf10","Wnt5a","Lix1")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
library(ComplexHeatmap)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = mesoderm_subset, addDOC = FALSE)

#### Visualize as featureplots
p11 <- plotEmbedding(
  ArchRProj = mesoderm_subset, 
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
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = mesoderm_subset, 
        addDOC = FALSE, width = 5, height = 5)

######## 7.5 Marker Genes Imputation with MAGIC (only for GAS not GIS)
mesoderm_subset <- addImputeWeights(mesoderm_subset, reducedDims = "Harmony")

p13 <- plotEmbedding(
  ArchRProj = mesoderm_subset, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(mesoderm_subset)
)

#Rearrange for grid plotting
p14 <- lapply(p11, function(x){
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
        ArchRProj = mesoderm_subset, 
        addDOC = FALSE, width = 5, height = 5)


###### Browser Tracks
p_browser <- plotBrowserTrack(
  ArchRProj = mesoderm_subset, 
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
        ArchRProj = mesoderm_subset, 
        addDOC = FALSE, width = 5, height = 5)

## Launch browser
#ArchRBrowser(mesoderm_subset)

##8.3 GAS spot sheck and manual assignment
markerGenes_spotmarkerGenes_spot <-c("Smoc1","Smoc2")
p_spotcheck <- plotEmbedding(
  ArchRProj = mesoderm_subset,
  colorBy = "GeneScoreMatrix",
  name = markerGenes_spotmarkerGenes_spot,
  continuousSet = "horizonExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(mesoderm_subset)
)
p_spotcheck$Smoc1

###################################################################################################
## at this point go to 3b and re-clustering BEFORE integrating with scRNAseq
## then, when happy-ish with clustering results, run lines 21-119 again then save:
saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/mesoderm_subset", load = TRUE)
###################################################################################################




# ################### Chapter 8 Defining Cluster Identity with scRNA-seq
# ## This now uses the scRNAseq from the same tissue used for scATACseq 
# library(Seurat)
# scRNA <- readRDS("/Users/sranade/scRNA-seq/2019_scRNA-seq/WT_B6/E775_E105/WT_atlas_11idents_201028.RDS")
# scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
# scRNA$clusters <- scRNA@active.ident
# table(scRNA@active.ident)
# # CNC      SHF_Progenitors             Endoderm           Epicardium        Cardiomyocyte     ParaxialMesoderm 
# # 2557                 7067                 9407                 2099                 8194                 1816 
# # Endocardium             Ectoderm LateralPlateMesoderm                EndMT                Blood 
# # 2618                 5101                 1378                  560                  827
# 
# 
# seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
# colnames(colData(seRNA))
# table(colData(seRNA)$clusters)
# # spot check that its the same as line 93
# 
# #### 8.1.1 Unconstrained Integration, no imputation. 
# mesoderm_subset <- addGeneIntegrationMatrix(
#   ArchRProj = mesoderm_subset, 
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = seRNA,
#   addToArrow = TRUE,
#   groupRNA = "clusters",
#   nameCell = "predictedCell",
#   nameGroup = "predictedGroup",
#   nameScore = "predictedScore",
#   useImputation = F,
#   force = T
# )
# 
# pal <- paletteDiscrete(values = seRNA$clusters)
# p15 <- plotEmbedding(
#   mesoderm_subset, 
#   colorBy = "cellColData",
#   name = "predictedGroup",
#   embedding = "UMAP",
#   pal = pal)
# p15
# plotPDF(p15, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = mesoderm_subset, addDOC = FALSE, width = 5, height = 5)
# 
# saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/mesoderm_subset", load = TRUE)
# mesoderm_subset

########################
########################
########################

### go to temp intermediate script for manual spot-checking and reclustering, otherwise continue to step 4 after saving

########################
########################
########################


