###### pma 
CNC_subset <- addMotifAnnotations(ArchRProj = CNC_subset, motifSet = "homer", name = "Motif", force = T)

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest_PA3_Cardiac_NeuralCrest_E115,
  ArchRProj = CNC_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 1"
)
motifsUp
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

#plotPDF(ggUp, name = "E115_pa3_up", width = 4, height = 7, ArchRProj = CNC_subset, addDOC = FALSE)


motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest_PA3_Cardiac_NeuralCrest_E115,
  ArchRProj = CNC_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -1"
)
motifsDo
df_gained <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df_gained <- df_gained[order(df_gained$mlog10Padj, decreasing = TRUE),]
df_gained$rank <- seq_len(nrow(df_gained))
head(df_gained)
ggDo <- ggplot(df_gained, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df_gained[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

#plotPDF(ggDo, name = "E115_pma_down", width = 4, height = 7, ArchRProj = CNC_subset, addDOC = FALSE)
