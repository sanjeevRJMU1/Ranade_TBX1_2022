library(viridis)
library(ComplexHeatmap)
library(circlize)
# #### 7.3 Identifying Marker Genes
markersGS <- getMarkerFeatures(ArchRProj = mesoderm_subset,
                               useMatrix = "GeneScoreMatrix",
                               groupBy = "Clusters6",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon",
                               binarize = F,
                               maxCells = 1000)
#markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")


markerGenes <- c(
  "Foxf1","Osr1",
  "Foxd1","Meox2","Zic1",
  "Tbx1","Foxc1","Foxc2",
  "Lix1","Msx2","Alx1",
  "Dlx2","Dlx5",
  "Bdnf","Ptx3",
  "Rgs5","Isl1","Hand2",
  "Tbx18","Wt1",
  "Myog","Myf5","Eya4","Pax3",
  "Actc1","Tnnt2"
)
  
length(markerGenes)


##Making a heatmap:
mark_stats <- as.matrix(assay(markersGS))
row.names(mark_stats) <-markersGS@elementMetadata@listData[["name"]]
mark_stats.plot <- mark_stats[markerGenes,]
## Cap to 2 /-1
mark_stats.plot[mark_stats.plot > 2] <- 2
mark_stats.plot[mark_stats.plot < (-2)] <- -2
marker_labels = rowAnnotation(foo = anno_mark(
  at = grep(paste(markerGenes, collapse="|"),row.names(mark_stats.plot)),
  labels = row.names(mark_stats.plot)[grepl(paste(markerGenes, collapse="|"), row.names(mark_stats.plot))]),
  annotation_legend_param = list(direction = "vertical",
                                 nrow = 1, labels_gp = gpar(fontsize = 0.1),
                                 padding = unit(10, "mm")))
my_order <-c(
  "Posterior_Second_Heart_Field",
  "Paraxial_Mesoderm",
  "Cardiopharyngeal_Mesoderm",
  "Cardiopharyngeal_Mesenchyme",
  "Neural_Crest",
  "Cranial_Mesenchyme",
  "Smooth_Muscle",
  "Epicardium",
  "Anterior_Second_Heart_Field",
  "Cardiomyocyte")
mark_stats.plot<-mark_stats.plot[,my_order]

draw(Heatmap(matrix = mark_stats.plot,
             cluster_rows = F,
             cluster_columns = F,
             clustering_method_columns = "ward.D2",
             clustering_method_rows = "ward.D2",
             show_row_names = F,
             show_column_names = T,
             col=colorRamp2(c(-2, 0, 1), c("grey", "white", "red")),
             right_annotation = marker_labels,
             border = TRUE,
             heatmap_legend_param = list(direction = "vertical", title = "Z-score GS")), heatmap_legend_side = "right")

