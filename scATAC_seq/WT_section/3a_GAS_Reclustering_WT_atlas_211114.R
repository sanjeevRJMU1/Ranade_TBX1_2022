#### Chapters 7 & 8, scRNA-seq integration and cluster identification

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

######## 
# Re-load ArchR project!
wt_atlas <- loadArchRProject(path = "/Users/sranade/scATAC-seq/new_wt_section_211111")
wt_atlas
######## 


# #### 7.3 Identifying Marker Genes
markersGS <- getMarkerFeatures(ArchRProj = wt_atlas,
                               useMatrix = "GeneScoreMatrix",
                               groupBy = "Clusters",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon",
                               binarize = F,
                               maxCells = 1000)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")

markerGenes <- c("Tnnt2","Dlx5","Dlx2","Epcam","Sox2","Sox10","Pecam1","Hba-x","Phox2a","Phox2b",
                 "Lix1","Alx1","Msx1","Msx2","Gsc","Shh","Fgf8","Fgf10","Sema3c","Wnt2","Wnt4",
                 "Tbx2","Tbx3","Bmp2","Rspo3","Bmp4","Foxf1","Osr1","Hes5","Crabp1","Pax3","Pax9","Pax5",
                 "Meox1","Ebf1","Otx2","Gbx1","Foxd1","Foxc1","Foxc2","Foxa1","Smoc1","Hand1","Hand2",
                 "Wt1","Tbx18","Mab21l2","Six1","Six2","Isl1","Tbx1","Tagln","Twist1","Snai1","Tdgf1","Meox2","Mitf","Myf5","Sox8")

######## 7.5 Marker Genes Imputation with MAGIC (only for GAS not GIS)
wt_atlas<-addImputeWeights(wt_atlas,dimsToUse = 2:25)

p11 <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(wt_atlas)
)

plotPDF(plotList = p11, 
        name = "Plot-UMAPHarmony-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = wt_atlas, 
        addDOC = FALSE, width = 5, height = 5)

marker_genes_2 <- c("Tnnt2","Isl1","Osr1","Fgf8","Tbx1","Wt1","Mab21l2","Meox1","Hand2",
                    "Tdgf1","Rgs5","Fgf10","Irx4","Vsnl1","Cited1","Nr2f1","Foxd1","Foxc2",
                    "Alx1","Tlx1","Aplnr","Nrg1","Hand1","Tbx18","Foxf1","Bmp2","Rspo3",
                    "Dlx1","Dlx2","Dlx3","Dlx4","Dlx5","Dlx6",
                    "Tdgf1","Irx3","Irx1","Tbx2","Tbx3","Nkx2-5","Gata4",
                    "Lix1","Tcf21","Lhx1","Lhx2","Gbx1","Otx2","Sox10","Sox2","Myf5","Pax3","Pitx2",
                    "Twist1","Prrx1","Prrx2","Col1a1","Pdgfra","Ebf1","Ebf2","Lum","Sox9","Msx1","Msx2",
                    "Pax1","Pax9","Pax5","Wnt4","Hoxa4","Hoxb4","Hoxc4","Hoxd4","Bdnf","Asb4","Zic1","Foxg1","Foxp2",
                    "Tagln","Acta2","Des","Vcl")

p11_2 <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "GeneScoreMatrix", 
  name = marker_genes_2, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(wt_atlas)
)

plotPDF(plotList = p11_2, 
        name = "Plot-UMAPHarmony-Marker-Genes-W-Imputation_2.pdf", 
        ArchRProj = wt_atlas, 
        addDOC = FALSE, width = 5, height = 5)


# ###### Browser Tracks
p_browser <- plotBrowserTrack(ArchRProj = wt_atlas,
                              groupBy = "Clusters",
                              geneSymbol = markerGenes,
                              upstream = 10000,
                              downstream = 10000)
dev.off()
#grid::grid.newpage()
#grid::grid.draw(p_browser$Tbx20)
plotPDF(plotList = p_browser,
        name = "Plot-Tracks-Marker-Genes.pdf",
        ArchRProj = wt_atlas,addDOC = FALSE, width = 5, height = 5)

## Launch browser
#ArchRBrowser(wt_atlas)

###################################################################################################
## Consider re-clustering?
### 4. Iterative LSI (TF-IDF and SVD)
wt_atlas <- addIterativeLSI(
  ArchRProj = wt_atlas,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 5, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1, 0.2, 0.3,0.4), 
    n.start = 10, maxClusters = NULL), 
  varFeatures = 25000, 
  dimsToUse = 2:30,
  LSIMethod = 2,
  seed = 1,
  force = T
)

# # Batch correction with Harmony --> from now on, use "Harmony" for reducedDims and "UMAP" for embedding
wt_atlas <- addHarmony(
  ArchRProj = wt_atlas,
  reducedDims = "IterativeLSI",
  dimsToUse = 2:30,
  name = "Harmony",
  groupBy = "Sample", force = T
)

#### 5. Clustering
wt_atlas <- addClusters(
  input = wt_atlas,
  reducedDims = "Harmony",
  dimsToUse = 2:30,
  method = "Seurat",
  name = "Clusters",
  resolution = 1,maxClusters = NULL, force = T
)
wt_atlas <- addUMAP(
  ArchRProj = wt_atlas, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony",
  metric = "cosine",
  force = T
)

p9 <- plotEmbedding(ArchRProj = wt_atlas, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p10 <- plotEmbedding(ArchRProj = wt_atlas, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
plotPDF(p9,p10, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = wt_atlas, addDOC = FALSE, width = 5, height = 5)

table(wt_atlas$Clusters,wt_atlas$Sample)



p11 <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(wt_atlas)
)

plotPDF(plotList = p11, 
        name = "Plot-UMAPHarmony-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = wt_atlas, 
        addDOC = FALSE, width = 5, height = 5)

p11_2 <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "GeneScoreMatrix", 
  name = marker_genes_2, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(wt_atlas)
)

plotPDF(plotList = p11_2, 
        name = "Plot-UMAPHarmony-Marker-Genes-W-Imputation_2.pdf", 
        ArchRProj = wt_atlas, 
        addDOC = FALSE, width = 5, height = 5)


marker_genes_3 <- c("Tnnt2","Tbx1","Wt1","Mab21l2","Meox1","Pecam1", "Sox10","Sox2","Epcam",
                    "Alx1","Dlx1","Dlx2","Dlx3","Dlx4","Dlx5","Hoxa4","Hoxb4","Hoxc4","Hoxd4")

p11_3 <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "GeneScoreMatrix", 
  name = marker_genes_3, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(wt_atlas)
)

plotPDF(plotList = p11_3, 
        name = "Plot-UMAPHarmony-Marker-Genes-W-Imputation_3.pdf", 
        ArchRProj = wt_atlas, 
        addDOC = FALSE, width = 5, height = 5)



####
saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)



