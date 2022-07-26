###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
##### Re-start project on 05/24

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


######## 
# Re-load ArchR project!
Tbx1_E825 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_E825")
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
# C1  C10  C11  C12  C13  C14  C15  C16  C17  C18  C19   C2   C3   C4   C5   C6   C7   C8   C9 
# 290  176  112 1238  309 1441 1384  458 1969 1365  665  635  206 1036 1382  324   84 1157  947
# Clusters2 = manually annotated clusters
table(Tbx1_E825$Clusters2)
# AHF          CM_Atrium             CM_OFT       CM_Ventricle                CNC           Ectoderm 
# 479                770                482               1155                299               3655 
# Endocardium           Endoderm         Epicardium                LPM PharyngealMesoderm                 PM 
# 639               2978               1126                453                771               1225 
# pSHF 
# 1146 

# Clusters3 = manually annotated clusters, divided by WT v KO
table(Tbx1_E825$Clusters3)
# AHF_K                AHF_W          CM_Atrium_K          CM_Atrium_W             CM_OFT_K 
# 182                  297                  368                  402                  127 
# CM_OFT_W       CM_Ventricle_K       CM_Ventricle_W                CNC_K                CNC_W 
# 355                  397                  758                  259                   40 
# Ectoderm_K           Ectoderm_W        Endocardium_K        Endocardium_W           Endoderm_K 
# 2289                 1366                  300                  339                 1493 
# Endoderm_W         Epicardium_K         Epicardium_W                LPM_K                LPM_W 
# 1485                  514                  612                  216                  237 
# PharyngealMesoderm_K PharyngealMesoderm_W                 PM_K                 PM_W               pSHF_K 
# 496                  275                  630                  595                  642 
# pSHF_W 
# 504


########################################################################################################################
########################################################################################################################
# What I did was re-run the entire workflow from script #2 up to #5 (DAR). I tried different versions of 
# clustering and liked the harmony + dims 2:30. The peak calls right now are for harmony etc for Clusters 4
# overall, the conclusions are that the mesoderm progenitors are fucked, the myocytes are spared and one 
# population of endoderm cells (that have GAS for Tbx1) do show slight defects...
# also, i am really liking the whole aggr object + subset'ing. In the end, that may be the winner here!
# save project and move to hub on 03/22.
########################################################################################################################
########################################################################################################################
# re-start this on 05/16 and re-visit conclusions from E8.25 stage. Key question, include in paper? We have no 
# scRNA-seq for this time point. but maybe its ok. maybe the key conclusion is that mesoderm subsets are fucked
# but that does not mean any given cell type does not form...(unlike hand2 e.g.)

### as mentioned before, key may be to zoom in on mesoderm subset and be solid there. Then bring in neural crest at the very end
########################################################################################################################
########################################################################################################################
# 1. re-visit cluster annotation.
p9 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p10 <- plotEmbedding(ArchRProj = Tbx1_E825, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony", rastr = F)
ggAlignPlots(p9, p10, type = "h")
# 19 clusters, looks nicely resolved, maybe a touch under-clustered but otherwise no issues.
#################
# 2. subset E8.25 RNAseq to only include mesoderm
library(Seurat)
scRNA <- readRDS("/Users/sranade/scRNA-seq/2019_scRNA-seq/WT_B6/E85_aggr_201022/E85_aggr_mnn_annotated_201022.RDS")
DimPlot(scRNA, reduction = "umap", label = T)
scRNA <- RunUMAP(scRNA, dims = 1:15, reduction = "mnn", verbose = TRUE)
scRNA <- FindNeighbors(scRNA,reduction = "mnn", dims = 1:15, verbose = TRUE, force.recalc = T)
scRNA <- FindClusters(scRNA, verbose = TRUE, resolution = 1.2)
DimPlot(scRNA, reduction = "umap", label = T)

scRNA <- SubsetData(scRNA, ident.use = c(8,1,3,10,11,16,5,0,12,13))

scRNA <- RunUMAP(scRNA, dims = 1:12, reduction = "mnn", verbose = TRUE)
scRNA <- FindNeighbors(scRNA,reduction = "mnn", dims = 1:12, verbose = TRUE, force.recalc = T)
scRNA <- FindClusters(scRNA, verbose = TRUE, resolution = 1)
DimPlot(scRNA, reduction = "umap", label = T)
VlnPlot(scRNA, c("nFeature_RNA"))
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.4)

cluster01.markers <- FindMarkers(scRNA, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
cluster912.markers <- FindMarkers(scRNA, ident.1 = 9, ident.2 = 12, min.pct = 0.25)
cluster610.markers <- FindMarkers(scRNA, ident.1 = 6, ident.2 = 10, min.pct = 0.25)
cluster78.markers <- FindMarkers(scRNA, ident.1 = 7, ident.2 = 8, min.pct = 0.25)

new.cluster.ids <- c("ParaxialMesoderm",
                     "ParaxialMesoderm",
                     "JCF",
                     "pSHF",
                     "AVC_CM",
                     "Atrial_CM",
                     "Ventricle_CM",
                     "AHF",
                     "AHF",
                     "OFT_CM",
                     "Ventricle_CM",
                     "SHFProg",
                     "OFT_CM",
                     "PharyngealMesoderm")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
DimPlot(scRNA, reduction = "umap", label = T) + NoLegend()
saveRDS(scRNA, file = "/Users/sranade/scATAC-seq/Tbx1_E825/scRNA_new_annotation_210516.RDS")
write.csv(scRNA.markers, file = "/Users/sranade/scATAC-seq/Tbx1_E825/scRNA_markers.csv")

##############
# 3. subset scATACseq object to include only mesoderm population

# Subset clusters of interest in the hackiest fucking way 
# clone object
mesoderm_subset <- Tbx1_E825

# assign Y/N values to only the clusters you want
library(plyr)
temp_cluster_names <- as.character(revalue(mesoderm_subset$Clusters, c(
  "C1"="False",
  "C2"="False",
  "C3"="False",
  "C4"="False",
  "C5"="False",
  "C6"="False",
  "C7"="False",
  "C8"="True",
  "C9"="True",
  "C10"="False",
  "C11"="False",
  "C12"="False",
  "C13"="False",
  "C14"="True",
  "C15"="True",
  "C16"="True",
  "C17"="False",
  "C18"="True",
  "C19"="True"
)))
length(unique(temp_cluster_names))
# put some random clusters number for now. clusters5
mesoderm_subset$Clusters5 <- temp_cluster_names
mesoderm_subset@cellColData@listData[["Clusters5"]] <- temp_cluster_names
table(mesoderm_subset$Clusters5,mesoderm_subset$Sample)
# a_KO_1 b_WT_1 c_KO_2 d_KO_3 e_WT_2 f_WT_3 g_WT_4
# False    847   1506   2396   1186    381    955    490
# True     743   1681   1696   1045    531    933    788

## subset
idxSample <- BiocGenerics::which(mesoderm_subset$Clusters5 %in% "True")
cellsSample <- mesoderm_subset$cellNames[idxSample]
subsetArchRProject(
  ArchRProj = mesoderm_subset,
  cells = cellsSample,
  outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210516", force = T
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


##############
# 4. Make new project, new scripts and go back essentially to step 2
##############






