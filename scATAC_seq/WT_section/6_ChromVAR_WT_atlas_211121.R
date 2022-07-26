#### Chapter 13: ChromVAR --> the featureplot of motifs is a super nice visual
#### 

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

# load object
wt_atlas <- loadArchRProject(path = "/Users/sranade/scATAC-seq/new_wt_section_211111")
wt_atlas


#### Chapter 13 ChromVAR Deviatons Enrichment with ArchR --> tried multiple motifsets. stick with homer.
wt_atlas <- addMotifAnnotations(ArchRProj = wt_atlas, motifSet = "homer", name = "Motif", force = T)

wt_atlas <- addBgdPeaks(wt_atlas, force = T)

wt_atlas <- addDeviationsMatrix(
  ArchRProj = wt_atlas, 
  peakAnnotation = "Motif",
  force = TRUE
)
#2022-07-01 21:22:43 : Completed Computing Deviations!, 372.057 mins elapsed.

plotVarDev <- getVarDeviations(wt_atlas, name = "MotifMatrix", plot = TRUE)
plotVarDev

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = wt_atlas, addDOC = FALSE)

# extract subset of motifs for downstream analysis
markerMotifs <- c("z:AP.2alpha.AP2_3",
                  "z:Isl1.Homeobox_139",
                  "z:RXR.NR..DR1_251")

########## 
## Feature plots of motif deviations. VERY USEFUL!!
########## 
p <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(wt_atlas)
)

p2 <- lapply(p, function(x){
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
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(p, name = "Plot-Motif-Featureplots-w-Imputation", width = 5, height = 5, ArchRProj = wt_atlas, addDOC = FALSE)


##### Draw Ridge plots
p <- plotGroups(ArchRProj = wt_atlas, 
                groupBy = "Clusters2", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(wt_atlas)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = wt_atlas, addDOC = FALSE)

## To see how these TF deviation z-scores compare to the inferred gene expression via gene scores of the 
## corresponding TF genes, we can overlay the gene scores for each of these TFs on the UMAPHarmony embedding.
markerRNA <- getFeatures(wt_atlas, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")

p <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "GeneScoreMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(wt_atlas)
)
p2 <- lapply(p, function(x){
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
do.call(cowplot::plot_grid, c(list(ncol = 2),p2))

markerRNA <- getFeatures(wt_atlas, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
p <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "GeneIntegrationMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAPHarmony",
  continuousSet = "blueYellow",
  imputeWeights = NULL
)
p2 <- lapply(p, function(x){
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
do.call(cowplot::plot_grid, c(list(ncol = 2),p2))

saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)
# numberOfCells(1): 64956
# medianTSS(1): 12.193
# medianFrags(1): 48811.5





