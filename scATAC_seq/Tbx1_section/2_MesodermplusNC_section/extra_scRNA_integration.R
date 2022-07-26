################### Chapter 8 Defining Cluster Identity with scRNA-seq
## This now uses the scRNAseq from the same tissue used for scATACseq
library(Seurat)
scRNA <- readRDS("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_mesoderm_cnc_annotated_210830.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)

seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)
# spot check that its the same as line 93

#### 8.1.1 Unconstrained Integration, no imputation.
mesoderm_subset <- addGeneIntegrationMatrix(
  ArchRProj = mesoderm_subset,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  useImputation = F,
  force = T
)

pal <- paletteDiscrete(values = seRNA$clusters)
p15 <- plotEmbedding(
  mesoderm_subset,
  colorBy = "cellColData",
  name = "predictedGroup",
  embedding = "UMAPHarmony",
  pal = pal)
p15
plotPDF(p15, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = mesoderm_subset, addDOC = FALSE, width = 5, height = 5)



### jaccard index
library(tidyr)
library(RColorBrewer)
#Now let's look at the overlap between the annotations.

#We will calculate a Jaccard index between all possible overlaps. For plotting, we only keep sets showing at least JI >= 0.15 with another set.

#### Jeeves edit: I renamed columuns and changed color of heatmap (& i wanted to show numbers on heatmap)
# requires you to have already run RNA integration on your atac-seq object


atac_mesoderm_subset_ji <- Reduce(bind_rows,lapply(unique(mesoderm_subset$Clusters6), function(SR_label) {
  df <- Reduce(bind_rows,lapply(unique(mesoderm_subset$predictedGroup), function(CCA_label) {
    ji <- sum(mesoderm_subset$Clusters6 == SR_label & mesoderm_subset$predictedGroup == CCA_label)/sum(mesoderm_subset$Clusters6 == SR_label | mesoderm_subset$predictedGroup ==CCA_label)
    df <- data.frame(CCA_cell_type=CCA_label, JI=ji, stringsAsFactors = F)
    return(df)
  }))
  df$SR_label_cell_type <- SR_label
  return(df)
}))

atac_mesoderm_subset_ji_mat <- spread(atac_mesoderm_subset_ji, key = SR_label_cell_type, value = JI)
row.names(atac_mesoderm_subset_ji_mat) <- atac_mesoderm_subset_ji_mat$CCA_cell_type
atac_mesoderm_subset_ji_mat <- atac_mesoderm_subset_ji_mat[,2:ncol(atac_mesoderm_subset_ji_mat)]
summary(apply(atac_mesoderm_subset_ji_mat, 1, max))
summary(apply(atac_mesoderm_subset_ji_mat, 2, max))
pheatmap::pheatmap(atac_mesoderm_subset_ji_mat, clustering_method = NULL, color = colorRampPalette((brewer.pal(n = 9, name ="Blues")))(100))
i <- apply(atac_mesoderm_subset_ji_mat, 1, max) >= 0.15
j <- apply(atac_mesoderm_subset_ji_mat, 2, max) >= 0.15
pheatmap::pheatmap(atac_mesoderm_subset_ji_mat[i,j], cluster_rows = F, cluster_cols = F, color = colorRampPalette((brewer.pal(n = 9, name ="Reds")))(100), display_numbers = T, number_color = "White")
