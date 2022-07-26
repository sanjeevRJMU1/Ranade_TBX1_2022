#### Export files for presentation
plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Clusters6", embedding = "UMAPHarmony")
p1<-plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p2<-plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = CNC_subset, addDOC = FALSE, width = 5, height = 5)


# table
table(CNC_subset$Clusters6)
table(CNC_subset$Clusters,CNC_subset$Sample)
# GAS
CNC_subset <- addImputeWeights(CNC_subset, reducedDims = "Harmony",dimsToUse = 2:20)
markerGenes<-c("Emx2","Smoc1","Ebf1","Dlx1","Dlx2","Dlx3","Dlx5","Dlx6",
               "Hand2","Hand1","Dlk1","Msx2","Foxf1","Gata3","Rgs5","Isl1",
               "Pou3f3","Barx1","Sox2","Sox10","Cdh19",
               "Hoxb3","Hoxd4","Six2","Hoxa3","Tbx2","Meox1","Hoxa4","Mef2c","Hoxb4")

markerGenes_2<-c("Lrp1","Nkx2-5","Erbb3","Cdh19")
p22 <- plotEmbedding(
  ArchRProj = CNC_subset, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes_2, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(CNC_subset)
)
#p22$Lrp1
plotPDF(plotList = p22, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation_p22.pdf", 
        ArchRProj = CNC_subset, 
        addDOC = FALSE, width = 5, height = 5)





####
# get the groupBW files
getGroupBW(ArchRProj = CNC_subset, groupBy = "Clusters7", tileSize = 100, ceiling = 4, normMethod = "ReadsInTSS")















