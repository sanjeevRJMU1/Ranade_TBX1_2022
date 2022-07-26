###### Tbx1 E9.25 scATAC-seq. 2 WT (Embryos 1 and 4) and 2 KO (Embryos 3 and 5)
#### Chapter 15: Integrative Analysis with ArchR: Cis-co-accessibility, peak-to-gene linkages and TF regulators

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

######## 
# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517")
mesoderm_subset
######## 
# re-run RNA integration
library(Seurat)
scRNA <- readRDS("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_mesoderm_cnc_annotated_210830.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)

seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)

#
mesoderm_subset <- addGeneIntegrationMatrix(
  ArchRProj = mesoderm_subset,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seRNA,
  addToArrow = T,
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
#plotPDF(p15, name = "Plot-UMAP-RNA-Integration_220124.pdf", ArchRProj = mesoderm_subset, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517", load = TRUE)





##### Peak to gene linkages
mesoderm_subset <- addPeak2GeneLinks(
  ArchRProj = mesoderm_subset,
  dimsToUse = 2:20,
  addEmpiricalPval = T,
  maxDist = 250000,
  reducedDims = "Harmony"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = mesoderm_subset,
  corCutOff = 0.4,
  resolution = 1,
  returnLoops = TRUE
)

p2g[[1]]

## 15.3.1 Plotting browser tracks with peak-to-gene links --> play with this...
markerGenes  <- c("Aplnr")

p <- plotBrowserTrack(
  ArchRProj = mesoderm_subset, 
  groupBy = "Clusters3", 
  geneSymbol = markerGenes, 
  upstream = 200000,
  downstream = 200000,
  loops = getPeak2GeneLinks(mesoderm_subset)
)

grid::grid.newpage()
grid::grid.draw(p$Aplnr)

saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517", load = TRUE)

######################################################################################## 
#### export peaks
# Filtered Low Prediction Score Cells (421 of 64956, 0.006), 0.024 mins elapsed.
# set peak to gene corcutoff at 0.5

p2g <- getPeak2GeneLinks(
  ArchRProj = mesoderm_subset,
  corCutOff = 0.7,
  FDRCutOff = 1e-04,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1,
  returnLoops = F
)

p2g
## when resolution set to 1, number of links are 67606 and changing resolution has little effect. It's the curCutOff that makes a big deal
## corcutoff 0.4, links = 113395 (wt_atlas), 81522 (mesoderm)
## corcutoff 0.5, links = 67606, 41404 (mesoderm)
## corcutoff 0.6, links = 35299, 20129 (mesoderm)
## corcutoff 0.7, links = 17103, 8558 (mesoderm)

# set to 0.4 since mesoderm subset is not as predictive as whole atlas


### exporting p2g dataframe
# if setting returnloops = T, you get a Granges object with 67k ranges of different lengths, which = the peak to TSS of linked gene
p2g <- getPeak2GeneLinks(
  ArchRProj = mesoderm_subset,
  corCutOff = 0.4,
  FDRCutOff = 1e-04,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1,
  returnLoops = T
)

p2g[[1]]
#GRanges object with 81484 ranges and 2 metadata columns:

# if setting returnloops = F, you get a data frame with 67k peaks and an idxATAC value for each
p2g_F <- getPeak2GeneLinks(
  ArchRProj = mesoderm_subset,
  corCutOff = 0.4,
  FDRCutOff = 1e-04,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1,
  returnLoops = F
)

p2g_F
#DataFrame with 81522 rows and 8 columns

metadata(p2g_F)[[1]]
# GRanges object with 475450 ranges and 0 metadata columns: (this is all the peaks in the peak matrix)
## the key is to match the idxATAC in p2g_F with the metadata object of p2g_F

p2g_F_df <- metadata(p2g_F)[[1]]

p2g_F_peaks_export <- data.frame(seqnames=seqnames(p2g_F_df),
                                 starts=start(p2g_F_df)-1,
                                 ends=end(p2g_F_df),
                                 names=c(rep(".", length(p2g_F_df))),
                                 strands=p2g_F_df@strand@values)

p2g_F_list <-p2g_F@listData[["idxATAC"]]
number_uni <-length(unique(p2g_F$idxATAC))
# 53,448
### this shows that the same peak can sometimes be linked to a different gene, explaining the lower number of unique values
subset_df <- p2g_F_peaks_export[rownames(p2g_F_peaks_export) %in% p2g_F$idxATAC,]
# this is a subsetted dataframe of unique p2g links (>0.4 cor) (but does not include the correlation value as a column!)


# format in same manner as you would use for homer 
subset_df <- subset_df %>% mutate(id = row_number())
subset_df$blankVar <- NA
subset_df_homerformat <- subset_df[,c(1,2,3,6,7,5)]
write.table(subset_df_homerformat, file="/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/p2g_exports_220203/p2g_linkages.bed", quote=F, sep="\t", row.names=F, col.names=F)
















