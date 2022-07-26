#### This is a temp script that fits between #2 and #3 ONLY for re-clustering
## run scripts 2 and 3 until the spot checking part, then run this temp script and then continue with #4
## bc its a temp script, I do not include the standard loading and closing lines. This is not meant to be the part
## where you start/re-start. Also, all plot saving lines are moved to the end so that you only save it once.

### 4. Iterative LSI (TF-IDF and SVD)
Tbx1_E825 <- addIterativeLSI(
  ArchRProj = Tbx1_E825,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 3, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2,0.5), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:12,
  LSIMethod = 2,
  seed = 1,
  force = T
)

#### 5. Clustering
Tbx1_E825 <- addClusters(
  input = Tbx1_E825,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.9, force = T
)


#### 6. Single-cell Embeddings
### UMAP with IterativeLSI
Tbx1_E825 <- addUMAP(
  ArchRProj = Tbx1_E825, 
  reducedDims = "IterativeLSI", 
  name = "UMAP",
  metric = "cosine",
  force = T
)
plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")




p16 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p17 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p16, p17, type = "h")


table(Tbx1_E825$Clusters,Tbx1_E825$Sample)

########################
########################
########################


##8.3 GAS spot sheck and manual assignment
markerGenes_spotmarkerGenes_spot <-c("Tbx1","Wt1","Meox1","Pmp22","Fgf8","Foxf1","Tbx18","Tnnt2","Irx4","Vsnl1","Tdgf1","Six2","Ebf1","Foxd1","Tlx1","Pecam1","Sox2","Epcam","Dlx2")
p_spotcheck <- plotEmbedding(
  ArchRProj = Tbx1_E825,
  colorBy = "GeneScoreMatrix",
  name = markerGenes_spotmarkerGenes_spot,
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(Tbx1_E825)
)
p_spotcheck$Pecam1


# clusters
p10

## My manual labeling --> make sure this is a separate metadata column from existing clusters
library(plyr)
temp_cluster_names <- as.character(revalue(Tbx1_E825$Clusters, c(
  "C1"="Endoderm",
  "C2"="Endoderm",
  "C3"="Endoderm",
  "C4"="Endoderm",
  "C5"="Endoderm",
  "C6"="Ectoderm",
  "C7"="CNC",
"C8"="CNC",
"C9"="Ectoderm",
"C10"="Ectoderm",
"C11"="Ectoderm",
"C12"="Ectoderm",
"C13"="Endoderm",
"C14"="Endoderm",
"C15"="LPM",
"C16"="CM_Ventricle",
"C17"="CM_OFT",
"C18"="CM_Atrium",
  "C19"="Epicardium",
  "C20"="Epicardium",
  "C21"="AHF",
  "C22"="pSHF",
  "C23"="Endocardium",
  "C24"="PM",
  "C25"="PharyngealMesoderm")))
length(unique(temp_cluster_names))
Tbx1_E825$Clusters2 <- temp_cluster_names
Tbx1_E825@cellColData@listData[["Clusters2"]] <- temp_cluster_names
table(Tbx1_E825$Clusters2,Tbx1_E825$Sample)

#                   a_KO_1 b_WT_1 c_KO_2 d_KO_3 e_WT_2 f_WT_3 g_WT_4
# AHF                    39    143     94     49     21     68     65
# CM_Atrium              90    179    156    122     63     71     89
# CM_OFT                 34     96     58     35     82     92     85
# CM_Ventricle          106    267    174    117    182    107    202
# CNC                    22     23    198     39      0      1     16
# Ectoderm              445    764   1200    644     97    448     57
# Endocardium            68    119    135     97     62     79     79
# Endoderm              293    556    817    383    210    399    320
# Epicardium            117    252    231    166     27    150    183
# LPM                    36     76    138     42      0     52    109
# PharyngealMesoderm     98    113    301     97     31     88     43
# PM                    128    322    296    206     90    182      1
# pSHF                  114    277    294    234     47    151     29

p18 <- plotEmbedding(
  Tbx1_E825,
  colorBy = "cellColData",
  name = "Clusters2",
  embedding = "UMAP",
  pal = pal)

p18

### save all intermediate plots
p16 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p17 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p16, p17, type = "h")
plotPDF(p16,p17, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)



plotPDF(p18,
        name = "Plot-UMAP-Manual-Annotations.pdf",
        ArchRProj = Tbx1_E825,
        addDOC = FALSE, width = , height = 5)

########################
## In this workflow, separate WT and KO per cluster and add a new metadata column called Clusters3

table(Tbx1_E825$Clusters2,Tbx1_E825$Sample)

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters3
condition.vec.01 <- substr(Tbx1_E825$Sample,3,3)
unique(condition.vec.01)
condition.vec.02 <- as.vector(Tbx1_E825$Clusters2)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
# condition.vec
# AHF_K                AHF_W          CM_Atrium_K          CM_Atrium_W             CM_OFT_K             CM_OFT_W 
# 182                  297                  368                  402                  127                  355 
# CM_Ventricle_K       CM_Ventricle_W                CNC_K                CNC_W           Ectoderm_K           Ectoderm_W 
# 397                  758                  259                   40                 2289                 1366 
# Endocardium_K        Endocardium_W           Endoderm_K           Endoderm_W         Epicardium_K         Epicardium_W 
# 300                  339                 1493                 1485                  514                  612 
# LPM_K                LPM_W PharyngealMesoderm_K PharyngealMesoderm_W                 PM_K                 PM_W 
# 216                  237                  496                  275                  630                  595 
# pSHF_K               pSHF_W 
# 642                  504 

Tbx1_E825$Clusters3 <- as.character(condition.vec)
table(Tbx1_E825$Clusters3)


saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/SR39_Tbx1_E825", load = TRUE)
Tbx1_E825


########################
########################
# Proceed to step 4


