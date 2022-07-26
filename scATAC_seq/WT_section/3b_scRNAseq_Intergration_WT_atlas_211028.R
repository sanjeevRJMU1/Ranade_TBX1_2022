###################################################################################################

################### Chapter 8 Defining Cluster Identity with scRNA-seq
library(Seurat)
scRNA <- readRDS("/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119/wt_atlas_rna_annotated.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)

# Endoderm      Neural_Crest     Cardiomyocyte    SHF_Progenitor Paraxial_Mesoderm          Ectoderm        Epicardium       Endothelium 
# 10968              8060              7860             14795              1982              7707              2695              3499 
# LPM             EndMT             Blood 
# 1894               699              1371 


seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)
# spot check that its the same as line 93

#### 8.1.1 Unconstrained Integration, no imputation. 
wt_atlas <- addGeneIntegrationMatrix(
  ArchRProj = wt_atlas, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupRNA = "clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  dimsToUse = 2:25,
  nameScore = "predictedScore",
  useImputation = F,
  force = T
)

pal <- paletteDiscrete(values = seRNA$clusters)
p13 <- plotEmbedding(
  wt_atlas, 
  colorBy = "cellColData",
  name = "predictedGroup",
  embedding = "UMAPHarmony",
  pal = pal)
p13
plotPDF(p13, name = "Plot-UMAPHarmony-RNA-Integration.pdf", ArchRProj = wt_atlas, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)
wt_atlas

########################
########################
########################

### go to temp intermediate script for manual spot-checking and reclustering, otherwise continue to step 4 after saving

########################
########################
########################