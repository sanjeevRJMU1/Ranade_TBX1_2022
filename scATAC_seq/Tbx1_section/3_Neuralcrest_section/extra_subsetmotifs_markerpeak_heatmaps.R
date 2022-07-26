## Marker Peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = CNC_subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters6",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks


heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.5",
  transpose = F,invert = F
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")





# motif enrichment
CNC_subset <- addMotifAnnotations(ArchRProj = CNC_subset, motifSet = "homer", name = "Motif", force = T)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CNC_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 2"
)
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")




plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = CNC_subset, addDOC = FALSE)

motif_df <- as.data.frame(enrichMotifs@assays$data$mlog10Padj)
write.table(motif_df, file="/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/cluster_motifs/enrichmotif_logP_homer_new.txt", quote=F, sep="\t", row.names=T, col.names=T)

## subset motifs for custom heatmap
object2subset <- enrichMotifs
## create vector with FULL name matches to rows of interest
rownames2match <- c("FoxL2.Forkhead_93",
                    "Gata6.Zf_109",
                    "Tgif2.Homeobox_296",
                    "Smad2.MAD_256",
                    "Pitx1.Ebox.Homeobox.bHLH_224",
                    "ZBTB18.Zf_307",
                    "Olig2.bHLH_202",
                    "Brn1.POU.Homeobox_27",
                    "AP.2alpha.AP2_3",
                    "Sox10.HMG_259",
                    "RXR.NR..DR1_251",
                    "HOXA2.Homeobox_128",
                    "Pdx1.Homeobox_219")
## pull indices where rownames match - full string matches only
indices4subset <- which(rownames(object2subset) %in% rownames2match, arr.ind=TRUE)
## subset SummarizedExperiment by row index
SE_subset <- object2subset[indices4subset, 1:ncol(object2subset)]

### sanity checks
object2subset
SE_subset

heatmapEM <- plotEnrichHeatmap(SE_subset, transpose = F)
order_rows <- c("Cardiac_NeuralCrest","Craniofacial_NeuralCrest","Migratory_NeuralCrest","PA3_Cardiac_NeuralCrest")
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap_subset", width = 8, height = 6, ArchRProj = CNC_subset, addDOC = FALSE)

