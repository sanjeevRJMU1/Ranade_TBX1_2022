################### Chapter 8 Defining Cluster Identity with scRNA-seq
## This now uses the scRNAseq from the same tissue used for scATACseq
library(Seurat)
scRNA <- readRDS("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_CNC_annotated.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)
# CNC      SHF_Progenitors             Endoderm           Epicardium        Cardiomyocyte     ParaxialMesoderm
# 2557                 7067                 9407                 2099                 8194                 1816
# Endocardium             Ectoderm LateralPlateMesoderm                EndMT                Blood
# 2618                 5101                 1378                  560                  827


seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)
# spot check that its the same as line 93

#### 8.1.1 Unconstrained Integration, no imputation.
CNC_subset <- addGeneIntegrationMatrix(
  ArchRProj = CNC_subset,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seRNA,
  addToArrow = F,
  groupRNA = "clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  useImputation = F,
  force = T
)

pal <- paletteDiscrete(values = seRNA$clusters)
p15 <- plotEmbedding(
  CNC_subset,
  colorBy = "cellColData",
  name = "predictedGroup",
  embedding = "UMAPHarmony",
  pal = pal)
p15
plotPDF(p15, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = CNC_subset, addDOC = FALSE, width = 5, height = 5)

################
### jaccard index
library(tidyr)
library(RColorBrewer)
#Now let's look at the overlap between the annotations.

#We will calculate a Jaccard index between all possible overlaps. For plotting, we only keep sets showing at least JI >= 0.15 with another set.

#### Jeeves edit: I renamed columuns and changed color of heatmap (& i wanted to show numbers on heatmap)
# requires you to have already run RNA integration on your atac-seq object


atac_CNC_subset_ji <- Reduce(bind_rows,lapply(unique(CNC_subset$Clusters6), function(SR_label) {
  df <- Reduce(bind_rows,lapply(unique(CNC_subset$predictedGroup), function(CCA_label) {
    ji <- sum(CNC_subset$Clusters6 == SR_label & CNC_subset$predictedGroup == CCA_label)/sum(CNC_subset$Clusters6 == SR_label | CNC_subset$predictedGroup ==CCA_label)
    df <- data.frame(CCA_cell_type=CCA_label, JI=ji, stringsAsFactors = F)
    return(df)
  }))
  df$SR_label_cell_type <- SR_label
  return(df)
}))

atac_CNC_subset_ji_mat <- spread(atac_CNC_subset_ji, key = SR_label_cell_type, value = JI)
row.names(atac_CNC_subset_ji_mat) <- atac_CNC_subset_ji_mat$CCA_cell_type
atac_CNC_subset_ji_mat <- atac_CNC_subset_ji_mat[,2:ncol(atac_CNC_subset_ji_mat)]
summary(apply(atac_CNC_subset_ji_mat, 1, max))
summary(apply(atac_CNC_subset_ji_mat, 2, max))
pheatmap::pheatmap(atac_CNC_subset_ji_mat, clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Blues")))(100))
i <- apply(atac_CNC_subset_ji_mat, 1, max) >= 0.15
j <- apply(atac_CNC_subset_ji_mat, 2, max) >= 0.15
pheatmap::pheatmap(atac_CNC_subset_ji_mat[i,j], clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Reds")))(100), display_numbers = T, number_color = "White")


#Overall, the two annotations agree well. 
##### Peak to gene linkages
CNC_subset <- addPeak2GeneLinks(
  ArchRProj = CNC_subset,
  maxDist = 250000,
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = CNC_subset,
  corCutOff = 0.4,
  resolution = 1,
  returnLoops = TRUE
)

p2g[[1]]

## representative marker genes and p2g links --> save all pdf as 5 in v 4 in
# migratory/progenitor NC 
markerGenes  <- c("Sox10")

p <- plotBrowserTrack(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters6", 
  geneSymbol = markerGenes, 
  upstream = 15000,
  downstream = 65000,
  loops = getPeak2GeneLinks(CNC_subset, corCutOff = 0.7)
)

grid::grid.newpage()
grid::grid.draw(p$Sox10)

# craniofacial
markerGenes  <- c("Smoc1")

p <- plotBrowserTrack(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters6", 
  geneSymbol = markerGenes, 
  upstream = 35000,
  downstream = 35000,
  loops = getPeak2GeneLinks(CNC_subset, corCutOff = 0.7)
)

grid::grid.newpage()
grid::grid.draw(p$Smoc1)

# cardiac
markerGenes  <- c("Foxp1")

p <- plotBrowserTrack(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters6", 
  geneSymbol = markerGenes, 
  upstream = 200000,
  downstream = 2000,
  loops = getPeak2GeneLinks(CNC_subset, corCutOff = 0.3)
)

grid::grid.newpage()
grid::grid.draw(p$Foxp1)


# PA3
markerGenes  <- c("Rxra")

p <- plotBrowserTrack(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters6", 
  geneSymbol = markerGenes, 
  upstream = 70000,
  downstream = 10000,
  loops = getPeak2GeneLinks(CNC_subset, corCutOff = 0.2)
)

grid::grid.newpage()
grid::grid.draw(p$Rxra)


##################
# PA3
markerGenes  <- c("Klf12")

p <- plotBrowserTrack(
  ArchRProj = CNC_subset, 
  groupBy = "Clusters6", 
  geneSymbol = markerGenes, 
  upstream = 500000,
  downstream = 1100000,
  normMethod = "none"
)

grid::grid.newpage()
grid::grid.draw(p$Klf12)


ArchRBrowser(CNC_subset)
