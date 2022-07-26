### 
# export GAS heatmap



library(viridis)
library(circlize)
# #### 7.3 Identifying Marker Genes
markersGS <- getMarkerFeatures(ArchRProj = CNC_subset,
                               useMatrix = "GeneScoreMatrix",
                               groupBy = "Clusters6",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon",
                               binarize = F,
                               maxCells = 1000)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")


markerGenes <- c("Erbb3","Sox10","Sox2","Pitx1","Dlx2","Dlx3","Dlx5","Pou3f3","Smoc1","Smoc2","Hand2","Dlk1","Isl1","Foxf1","Rgs5","Hoxa3","Hoxb3","Hoxd4","Hoxa4","Hoxb4","Hoxc4","Six2")
length(markerGenes)


##Making a heatmap:
mark_stats <- as.matrix(assay(markersGS))
row.names(mark_stats) <-markersGS@elementMetadata@listData[["name"]]
mark_stats.plot <- mark_stats[markerGenes,]
## Cap to 2 /-1
mark_stats.plot[mark_stats.plot > 2] <- 2
mark_stats.plot[mark_stats.plot < (-1)] <- -1
marker_labels = rowAnnotation(foo = anno_mark(
  at = grep(paste(markerGenes, collapse="|"),row.names(mark_stats.plot)),
  labels = row.names(mark_stats.plot)[grepl(paste(markerGenes, collapse="|"), row.names(mark_stats.plot))]),
  annotation_legend_param = list(direction = "vertical",
                                 nrow = 1, labels_gp = gpar(fontsize = 0.1),
                                 padding = unit(10, "mm")))
my_order <-c("Migratory_NeuralCrest","Craniofacial_NeuralCrest","Cardiac_NeuralCrest","PA3_Cardiac_NeuralCrest")
mark_stats.plot<-mark_stats.plot[,my_order]

draw(Heatmap(matrix = mark_stats.plot,
             cluster_rows = F,
             cluster_columns = F,
             clustering_method_columns = "ward.D2",
             clustering_method_rows = "ward.D2",
             show_row_names = F,
             show_column_names = T,
             col=colorRamp2(c(-1, 0, 1), c("grey", "white", "red")),
             right_annotation = marker_labels,
             border = TRUE,
             heatmap_legend_param = list(direction = "vertical", title = "Z-score GS")), heatmap_legend_side = "right")

