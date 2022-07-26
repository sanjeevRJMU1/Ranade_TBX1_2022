#### Chapters 11 and 12: Marker Peaks and Motif Enrichment
#### 

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

# load object
wt_atlas <- loadArchRProject(path = "/Users/sranade/scATAC-seq/new_wt_section_211111")
wt_atlas


## 11. Marker Peaks (can use either PeakMatrix or GeneScoreMatrix!)
markersPeaks <- getMarkerFeatures(
  ArchRProj = wt_atlas, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
# class: SummarizedExperiment 
# dim: 513172 10 
# metadata(2): MatchInfo Params
# assays(7): Log2FC Mean ... AUC MeanBGD
# rownames(513172): 1 2 ... 513171 513172
# rowData names(4): seqnames idx start end
# colnames(10): Blood Cardiomyocyte ... MesodermProgenitor NeuralCrest
# colData names(0):


markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1.5", returnGR = TRUE)
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.5",
  transpose = F,
  nLabel = 1,
  nPrint = 1)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = wt_atlas, addDOC = FALSE)

#### export the markerlist genomic ranges object results for a specific cluster to a data frame
ahf_df <- as.data.frame(markerList@listData[["EndMT"]])
write.csv(ahf_df, file = "/Users/sranade/scATAC-seq/wt_atlas/jeeves_exports/ahf.csv")

## export all peaks from peak heatmap
peak_matrix_heatmap <- heatmapPeaks@ht_list[["Row Z-Scores
113643 features
PeakMatrix"]]@matrix
peak_matrix_heatmap<-as.data.frame(peak_matrix_heatmap)
write.csv(peak_matrix_heatmap, file = "/Users/sranade/scATAC-seq/new_wt_section_211111/ExportsforFigures/peakmatrix_heatmap.csv",row.names = T)


##### 11.2.3 Marker Peaks in Browser Tracks
p <- plotBrowserTrack(
  ArchRProj = wt_atlas,
  groupBy = "Clusters2",
  geneSymbol = c("Wnt2"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 2", returnGR = TRUE),
  upstream = 100000,
  downstream = 100000
)
grid::grid.newpage()
grid::grid.draw(p$Wnt2)

# visual for only plotting specific clusters
p <- plotBrowserTrack(
  ArchRProj = wt_atlas,
  groupBy = "Clusters2",
  geneSymbol = c("Tshz2"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[c("Ventricle_Myocardium","Endocardium")],
  upstream = 10000,
  downstream = 10000
)
grid::grid.newpage()
grid::grid.draw(p$Tshz2)

#### 11.3 Pairwise Testing Between Groups
# See script 5b


## 12: Motif Enrichment in Differential Peaks 
# actually I think 12.2 is more intuitive first but 12.1 is more what people want
## 12.1 Enriched motifs in differentially enriched peaks in cluster a vs b
# wt_atlas <- addMotifAnnotations(ArchRProj = wt_atlas, motifSet = "homer", name = "Motif", force = T)
# 
# motifsUp <- peakAnnoEnrichment(
#   seMarker = markerTest,
#   ArchRProj = wt_atlas,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
# motifsUp
# 
# # plot
# df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
# df <- df[order(df$mlog10Padj, decreasing = TRUE),]
# df$rank <- seq_len(nrow(df))
# head(df)
# 
# 
# ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
#   geom_point(size = 1) +
#   ggrepel::geom_label_repel(
#     data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
#     size = 1.5,
#     nudge_x = 2,
#     color = "black"
#   ) + theme_ArchR() + 
#   ylab("-log10(P-adj) Motif Enrichment") + 
#   xlab("Rank Sorted TFs Enriched") +
#   scale_color_gradientn(colors = paletteContinuous(set = "comet"))
# 
# ggUp
# 
# motifsDo <- peakAnnoEnrichment(
#   seMarker = markerTest,
#   ArchRProj = wt_atlas,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
# )
# motifsDo
# df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
# df <- df[order(df$mlog10Padj, decreasing = TRUE),]
# df$rank <- seq_len(nrow(df))
# head(df)
# 
# ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
#   geom_point(size = 1) +
#   ggrepel::geom_label_repel(
#     data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
#     size = 1.5,
#     nudge_x = 2,
#     color = "black"
#   ) + theme_ArchR() + 
#   ylab("-log10(FDR) Motif Enrichment") +
#   xlab("Rank Sorted TFs Enriched") +
#   scale_color_gradientn(colors = paletteContinuous(set = "comet"))
# 
# ggDo
# 
# plotPDF(ggUp, ggDo, name = "Ventricle-vs-OFT-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = wt_atlas, addDOC = FALSE)

## 12.2 what motifs are enriched per cluster?
wt_atlas <- addMotifAnnotations(ArchRProj = wt_atlas, motifSet = "homer", name = "Motif", force = T)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = wt_atlas,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = wt_atlas, addDOC = FALSE)

saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)











