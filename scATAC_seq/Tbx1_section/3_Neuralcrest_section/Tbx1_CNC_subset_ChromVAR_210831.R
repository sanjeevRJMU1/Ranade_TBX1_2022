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

## 11. Marker Peaks (can use either PeakMatrix or GeneScoreMatrix!)
markersPeaks <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters6",
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


markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.5",
  transpose = F,
  nLabel = 1,
  nPrint = 1)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

CNC_subset <- addMotifAnnotations(ArchRProj = CNC_subset, motifSet = "cisbp", name = "Motif", force = T)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CNC_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = CNC_subset, addDOC = FALSE)
motif_df <- as.data.frame(enrichMotifs@assays$data$mlog10Padj)
write.table(motif_df, file="/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/cluster_motifs/enrichmotif_logP.txt", quote=F, sep="\t", row.names=T, col.names=T)


##
### TODO: add support for partial rownames
object2subset <- enrichMotifs
## create vector with FULL name matches to rows of interest
rownames2match <- c("GATA.SCL.Zf.bHLH_111",
                    "EKLF.Zf_65",
                    "Hnf1.Homeobox_125",
                    "Bach1.bZIP_16",
                    "Mef2d.MADS_164",
                    "Tbx20.T.box_284",
                    "Gata4.Zf_108",
                    "Nkx2.1.Homeobox_182",
                    "Sox2.HMG_261",
                    "OCT4.SOX2.TCF.NANOG.POU.Homeobox.HMG_198",
                    "Brn1.POU.Homeobox_27",
                    "HOXA2.Homeobox_128",
                    "c.Jun.CRE.bZIP_142",
                    "Fosl2.bZIP_86",
                    "Atf2.bZIP_11",
                    "Etv2.ETS_82",
                    "GRHL2.CP2_119",
                    "AP.2gamma.AP2_2",
                    "FOXA1.Forkhead_87",
                    "Six2.Homeobox_255",
                    "ERG.ETS_72",
                    "EWS.FLI1.fusion.ETS_84",
                    "ETV1.ETS_81",
                    "Fli1.ETS_85",
                    "Gata2.Zf_102",
                    "TEAD.TEA_294",
                    "NF1.CTF_176",
                    "Tlx..NR_297",
                    "TEAD2.TEA_292",
                    "Gata6.Zf_109",
                    "Tcf3.HMG_288",
                    "Foxo3.Forkhead_96",
                    "Pitx1.Ebox.Homeobox.bHLH_224",
                    "ZBTB18.Zf_307",
                    "NeuroD1.bHLH_173",
                    "Olig2.bHLH_202",
                    "AP.2alpha.AP2_3",
                    "Sox9.HMG_265",
                    "Lhx2.Homeobox_151",
                    "RXR.NR..DR1_251"
)
## pull indices where rownames match - full string matches only
indices4subset <- which(rownames(object2subset) %in% rownames2match, arr.ind=TRUE)
## subset SummarizedExperiment by row index
SE_subset <- object2subset[indices4subset, 1:ncol(object2subset)]

### sanity checks
object2subset
SE_subset

heatmapEM <- plotEnrichHeatmap(SE_subset, transpose = F)
order_rows <- c("Cardiomyocyte","SHF_Progenitorsenitor","LPM","Epicardium","EndMT","Endothelium","NeuralCrest","Endoderm","Ectoderm","Blood")
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = wt_atlas, addDOC = FALSE)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot", row_order= order_rows)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = wt_atlas, addDOC = FALSE)

























############ ChromVAR




CNC_subset <- addBgdPeaks(CNC_subset, force = T)

CNC_subset <- addDeviationsMatrix(
  ArchRProj = CNC_subset, 
  peakAnnotation = "Motif",
  force = TRUE
)
## 2020-09-09 08:38:39 : Completed Computing Deviations!, 19.982 mins elapsed.

plotVarDev <- getVarDeviations(CNC_subset, name = "MotifMatrix", plot = TRUE)
plotVarDev

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = CNC_subset, addDOC = FALSE)

# extract subset of motifs for downstream analysis
motifs <- c("Gata4","Tbx20","Pitx1","Foxa1","ERG","Hoxa2")
markerMotifs <- getFeatures(CNC_subset, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

##### Draw Ridge plots
p <- plotGroups(ArchRProj = CNC_subset, 
                groupBy = "Clusters2", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(CNC_subset)
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

plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = CNC_subset, addDOC = FALSE)

########## 
## Feature plots of motif deviations. VERY USEFUL!!
########## 
p <- plotEmbedding(
  ArchRProj = CNC_subset, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(CNC_subset)
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

## To see how these TF deviation z-scores compare to the inferred gene expression via gene scores of the 
## corresponding TF genes, we can overlay the gene scores for each of these TFs on the UMAP embedding.
markerRNA <- getFeatures(CNC_subset, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")

p <- plotEmbedding(
  ArchRProj = CNC_subset, 
  colorBy = "GeneScoreMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(CNC_subset)
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

markerRNA <- getFeatures(CNC_subset, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
p <- plotEmbedding(
  ArchRProj = CNC_subset, 
  colorBy = "GeneIntegrationMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP",
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







