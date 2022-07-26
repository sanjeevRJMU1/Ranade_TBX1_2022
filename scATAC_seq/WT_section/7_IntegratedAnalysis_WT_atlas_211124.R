###### All WT atlas time points: 2 x E775, 3 x E825, 2 x E925
#### Chapter 15: Integrative Analysis with ArchR: Cis-co-accessibility and peak-to-gene linkages

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

# load object
wt_atlas <- loadArchRProject(path = "/Users/sranade/scATAC-seq/new_wt_section_211111")
wt_atlas


## Cis-coaccessibility analysis
wt_atlas <- addCoAccessibility(
  ArchRProj = wt_atlas,
  dimsToUse = 2:25,
  reducedDims = "Harmony"
)

cA <- getCoAccessibility(
  ArchRProj = wt_atlas,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = T
)
cA[[1]]
#GRanges object with 354,981 ranges and 9 metadata columns:


## playing with the resolution will decrease the number of loops (presumably higher filter? or just illustrative?)
# cA <- getCoAccessibility(
#   ArchRProj = wt_atlas,
#   corCutOff = 0.5,
#   resolution = 10000,
#   returnLoops = TRUE
# )
# cA[[1]]

## 15.2.1 Plotting browser tracks of Co-accessibility
markerGenes  <- c("Sh3bgr","S1pr1","Foxd2")

p <- plotBrowserTrack(
  ArchRProj = wt_atlas, 
  groupBy = "Clusters2", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(wt_atlas, corCutOff = 0.8)
)

dev.off()
grid::grid.newpage()
grid::grid.draw(p$Tnnt2)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
        ArchRProj = wt_atlas, 
        addDOC = FALSE, width = 5, height = 5)

##### Peak to gene linkages --> prediction score cutoff = 0.4, consistent with ML pred score cutoff
wt_atlas <- addPeak2GeneLinks(
  ArchRProj = wt_atlas,
  dimsToUse = 2:25,
  reducedDims = "IterativeLSI",
  predictionCutoff = 0.4,
  addEmpiricalPval = T
)

p2g <- getPeak2GeneLinks(
  ArchRProj = wt_atlas,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = T
)

#metadata(p2g[[2]])

## 15.3.1 Plotting browser tracks with peak-to-gene links --> play with this...
markerGenes  <- c("Srrm1","Sh3bgr","Gipc2")

p <- plotBrowserTrack(
  ArchRProj = wt_atlas, 
  groupBy = "Clusters2", 
  geneSymbol = markerGenes, 
  upstream = 5000,
  downstream = 20000,normMethod = "none"
)

grid::grid.newpage()
grid::grid.draw(p$Gipc2)


saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/new_wt_section_211111 
# samples(18): o_E105_3 p_E105_4 ... a_E775_1 c_E825_1
# sampleColData names(1): ArrowFiles
# cellColData names(22): Sample TSSEnrichment ... predictedGroup predictedScore
# numberOfCells(1): 64956
# medianTSS(1): 12.193
# medianFrags(1): 48811.5


#### browser tracks of Gipc2, Sh3bgr and Srrm1 
p <- plotBrowserTrack(
  ArchRProj = wt_atlas, 
  groupBy = "Clusters2", 
  geneSymbol = markerGenes, 
  upstream = 5000,
  downstream = 20000,normMethod = "none")

  grid::grid.newpage()
  grid::grid.draw(p$Sh3bgr)
  

##
  p <- plotBrowserTrack(
    ArchRProj = wt_atlas, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 1000,
    loops = p2g,normMethod = "none"
  )
  
  grid::grid.newpage()
  grid::grid.draw(p$Srrm1)
  
##
markerGenes  <- c("Srrm1","Sh3bgr","Gipc2")
  
p <- plotBrowserTrack(
    ArchRProj = wt_atlas, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 40000,
    downstream = 5000,normMethod = "none"
  )
  
grid::grid.newpage()
grid::grid.draw(p$Gipc2)
  
########
# neural crest genes
markerGenes  <- c("Prr30")

p <- plotBrowserTrack(
  ArchRProj = wt_atlas, 
  groupBy = "Clusters2", 
  geneSymbol = markerGenes, 
  upstream = 500000,
  downstream = 5000,normMethod = "ReadsInTSS"
)

grid::grid.newpage()
grid::grid.draw(p$Prr30)

markerGenes  <- c("Epha4")

p <- plotBrowserTrack(
  ArchRProj = wt_atlas, 
  groupBy = "Clusters2", 
  geneSymbol = markerGenes, 
  upstream = 200000,
  downstream = 5000,normMethod = "ReadsInTSS"
)

grid::grid.newpage()
grid::grid.draw(p$Epha4)


markerGenes  <- c("Fam43a")

p <- plotBrowserTrack(
  ArchRProj = wt_atlas, 
  groupBy = "Clusters2", 
  geneSymbol = markerGenes, 
  upstream = 2000,
  downstream = 35000,normMethod = "ReadsInTSS"
)

grid::grid.newpage()
grid::grid.draw(p$Fam43a)




## 15.3.2 Plotting a heatmap of peak-to-gene links
# p <- plotPeak2GeneHeatmap(ArchRProj = wt_atlas, groupBy = "Clusters2",k=10)
# p
# plotPDF(p, 
#         name = "plotPeak2GeneHeatmap.pdf", 
#         ArchRProj = wt_atlas, 
#         addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)









