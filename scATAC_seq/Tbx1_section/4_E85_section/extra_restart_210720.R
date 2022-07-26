###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
##### Re-start project on 05/24

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


######## 
# Re-load ArchR project!
Tbx1_E825 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_E825_only/Tbx1_E825")
Tbx1_E825
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825 
# samples(7): c_KO_2 b_WT_1 ... f_WT_3 e_WT_2
# sampleColData names(1): ArrowFiles
# cellColData names(23): Sample TSSEnrichment ... ReadsInPeaks FRIP
# numberOfCells(1): 15178
# medianTSS(1): 11.38
# medianFrags(1): 51252
######## 

## Table of Clusters
# Clusters = simple clustering all cells
table(Tbx1_E825$Clusters)
# C1  C10  C11  C12  C13  C14  C15  C16  C17  C18  C19   C2  C20  C21  C22  C23  C24  C25   C3   C4   C5   C6   C7   C8   C9 
# 625 1333  131  296  103  295  234 1309  444  314  150  289 1038  660 1220  153  699  840 1162  714  205 1103  735  559  567

# Clusters2 = Clusters1 separated by WT and KO (this got over-written at some point, it used to be just the 
## manually annotated clusters)
table(Tbx1_E825$Clusters2)
# C1_K  C1_W C10_K C10_W C11_K C11_W C12_K C12_W C13_K C13_W C14_K C14_W C15_K C15_W C16_K C16_W C17_K C17_W C18_K C18_W C19_K 
# 292   333   471   862    56    75   157   139    47    56   253    42   111   123   654   655   216   228   202   112   112 
# C19_W  C2_K  C2_W C20_K C20_W C21_K C21_W C22_K C22_W C23_K C23_W C24_K C24_W C25_K C25_W  C3_K  C3_W  C4_K  C4_W  C5_K  C5_W 
# 38    96   193   611   427   421   239   823   397    50   103   274   425   393   447   581   581   453   261   116    89 
# C6_K  C6_W  C7_K  C7_W  C8_K  C8_W  C9_K  C9_W 
# 678   425   323   412   240   319   283   284 

# Clusters3 = deprecated, do NOT use this info --> i ended up over-writing this with manual cluster annotation (07/21)
table(Tbx1_E825$Clusters3)
### don't use

# Clusters3 = deprecated, do NOT use this info
table(Tbx1_E825$Clusters4)
### don't use

########################################################################################################################
## re-re-re-vist on 210720 for exporting figures for supp data
########################################################################################################################



# 1. re-visit cluster annotation.
p1 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p1, p2, type = "h")
# 25 clusters, appropriately resolved
#################
#### Visualize GAS
Tbx1_E825 <- addImputeWeights(Tbx1_E825, reducedDims = "Harmony", dimsToUse = 2:20)

markerGenes <- c("Tnnt2","Epcam","Sox2","Isl1","Osr1","Fgf8","Tbx1","Pecam1","Wt1","Dlx2","Hba-x","Mab21l2","Meox1",
                 "Tdgf1","Rgs5","Fgf10","Irx4","Vsnl1","Cited1","Nr2f1","Foxd1","Foxc2","Alx1","Tlx1","Aplnr","Nrg1","Hand1","Tbx18","Foxf1","Bmp2","Rspo3")

p13 <- plotEmbedding(
  ArchRProj = Tbx1_E825, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(Tbx1_E825)
)
plotPDF(plotList = p13, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = Tbx1_E825, 
        addDOC = FALSE, width = 5, height = 5)

#################
# library(Seurat)
# scRNA <- readRDS("/Users/sranade/scRNA-seq/2019_scRNA-seq/WT_B6/E85_aggr_201022/E85_aggr_mnn_annotated_210524.RDS")
# scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 3000)
# scRNA$clusters <- scRNA@active.ident
# table(scRNA@active.ident)
# # Ectoderm        VentricleCM   ParaxialMesoderm           AtrialCM           Endoderm PharyngealMesoderm 
# # 1792                700               1041                513               1403                435 
# # OFT_CM                JCF        Endocardium               pSHF                AHF                CNC 
# # 422                385                379                343                258                188 
# # AVC_CM 
# # 160  
# 
# seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
# colnames(colData(seRNA))
# table(colData(seRNA)$clusters)
# # spot check that its the same as line 93
# 
# #### 8.1.1 Unconstrained Integration, no imputation. 
# Tbx1_E825 <- addGeneIntegrationMatrix(
#   ArchRProj = Tbx1_E825, 
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = seRNA,
#   addToArrow = TRUE,
#   groupRNA = "clusters",
#   nameCell = "predictedCell",
#   nameGroup = "predictedGroup",
#   nameScore = "predictedScore",
#   nGenes = 3000,
#   useImputation = F,
#   force = T
# )
# 
# pal <- paletteDiscrete(values = seRNA$clusters)
# p15 <- plotEmbedding(
#   Tbx1_E825, 
#   colorBy = "cellColData",
#   name = "predictedGroup",
#   embedding = "UMAPHarmony",
#   pal = pal)
# p15
# plotPDF(p15, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)

##########################################
# call peaks on clusters1, wt and ko --> deleted this whole section, but kept it in previous iteration of 
# script. refer back to that. not needed for exporting figures so i cut it.

##########################################
# 
# table(Tbx1_E825$Clusters,Tbx1_E825$Sample)
# 
# # Break up each cluster into WT and KO, then add this as a new cellcoldata column called "Clusters3
# condition.vec.01 <- substr(Tbx1_E825$Sample,3,3)
# unique(condition.vec.01)
# condition.vec.02 <- as.vector(Tbx1_E825$Clusters)
# condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
# table(condition.vec)
# # condition.vec
# # C1_K  C1_W C10_K C10_W C11_K C11_W C12_K C12_W C13_K C13_W C14_K C14_W C15_K C15_W C16_K C16_W C17_K C17_W C18_K C18_W C19_K 
# # 292   333   471   862    56    75   157   139    47    56   253    42   111   123   654   655   216   228   202   112   112 
# # C19_W  C2_K  C2_W C20_K C20_W C21_K C21_W C22_K C22_W C23_K C23_W C24_K C24_W C25_K C25_W  C3_K  C3_W  C4_K  C4_W  C5_K  C5_W 
# # 38    96   193   611   427   421   239   823   397    50   103   274   425   393   447   581   581   453   261   116    89 
# # C6_K  C6_W  C7_K  C7_W  C8_K  C8_W  C9_K  C9_W 
# # 678   425   323   412   240   319   283   284  
# 
# Tbx1_E825$Clusters2 <- as.character(condition.vec)
# table(Tbx1_E825$Clusters2)
# 

############################
############################
# subset scATACseq object to include only mesoderm population

# Subset clusters of interest in the hackiest fucking way 
# clone object
mesoderm_subset <- Tbx1_E825

# assign Y/N values to only the clusters you want
library(plyr)
temp_cluster_names <- as.character(revalue(mesoderm_subset$Clusters, c(
  "C1"="False",
  "C2"="False",
  "C3"="True",
  "C4"="True",
  "C5"="False",
  "C6"="False",
  "C7"="False",
  "C8"="False",
  "C9"="True",
  "C10"="True",
  "C11"="True",
  "C12"="True",
  "C13"="False",
  "C14"="False",
  "C15"="True",
  "C16"="True",
  "C17"="True",
  "C18"="False",
  "C19"="False",
  "C20"="False",
  "C21"="False",
  "C22"="False",
  "C23"="False",
  "C24"="True",
  "C25"="True"
)))
length(unique(temp_cluster_names))
# Clusters3 = t/f. Clusters2 = WT/KO
mesoderm_subset$Clusters3 <- temp_cluster_names
mesoderm_subset@cellColData@listData[["Clusters3"]] <- temp_cluster_names
table(mesoderm_subset$Clusters3,mesoderm_subset$Sample)
# a_KO_1 b_WT_1 c_KO_2 d_KO_3 e_WT_2 f_WT_3 g_WT_4
# False    819   1428   2301   1144    365    920    472
# True     771   1759   1791   1087    547    968    806

## subset
idxSample <- BiocGenerics::which(mesoderm_subset$Clusters3 %in% "True")
cellsSample <- mesoderm_subset$cellNames[idxSample]
subsetArchRProject(
  ArchRProj = mesoderm_subset,
  cells = cellsSample,
  outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210524", force = T
)
## Dropping ImputeWeights Since You Are Subsetting Cells! ImputeWeights is a cell-x-cell Matrix!
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210516 
# samples(7): c_KO_2 b_WT_1 ... a_KO_1 e_WT_2
# sampleColData names(1): ArrowFiles
# cellColData names(25): Sample TSSEnrichment ... Clusters4 Clusters5
# numberOfCells(1): 7417
# medianTSS(1): 11.359
# medianFrags(1): 54455


####################################################################################
# coarsely annotate clusters
####################################################################################
library(plyr)
temp_cluster_names <- as.character(revalue(Tbx1_E825$Clusters, c(
  "C1"="Endothelium",
  "C2"="Endoderm",
  "C3"="Mes_Prog",
  "C4"="Mes_Prog",
  "C5"="Endoderm",
  "C6"="Endoderm",
  "C7"="Endoderm",
  "C8"="Endoderm",
  "C9"="Cardiomyocyte",
  "C10"="Cardiomyocyte",
  "C11"="Cardiomyocyte",
  "C12"="Mes_Prog",
  "C13"="Endoderm",
  "C14"="Ectoderm",
  "C15"="Mes_Prog",
  "C16"="Mes_Prog",
  "C17"="Mes_Prog",
  "C18"="Ectoderm",
  "C19"="Ectoderm",
  "C20"="Ectoderm",
  "C21"="Ectoderm",
  "C22"="Ectoderm",
  "C23"="Ectoderm",
  "C24"="Mes_Prog",
  "C25"="Mes_Prog"
)))
length(unique(temp_cluster_names))
# [1] 5
Tbx1_E825$Clusters3 <- temp_cluster_names
Tbx1_E825@cellColData@listData[["Clusters3"]] <- temp_cluster_names
table(Tbx1_E825$Clusters3,Tbx1_E825$Sample)


p4<-plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters3", embedding = "UMAPHarmony")
plotPDF(p4, name = "Plot-UMAP-Manual_Integration.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E825_only/Tbx1_E825", load = TRUE)






##############
# 4. Make new project, new scripts and go back essentially to step 2
##############






