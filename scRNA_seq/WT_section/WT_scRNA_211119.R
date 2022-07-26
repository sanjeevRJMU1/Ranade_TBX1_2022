# E7.75 - E11.5

setwd("/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119")

library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)

# Load the wt_atlas_rna dataset
wt_atlas_rna.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")


# Initialize the Seurat object with the raw (non-normalized data)
wt_atlas_rna <- CreateSeuratObject(counts = wt_atlas_rna.data, project = "wt_atlas_rna", min.cells = 3, min.features = 200)
wt_atlas_rna 
# An object of class Seurat 
# 27173 features across 99434 samples within 1 assay 
# Active assay: RNA (27173 features, 0 variable features)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
wt_atlas_rna[["percent.mt"]] <- PercentageFeatureSet(wt_atlas_rna, pattern = "^mt-")

# Add sample IDs to metadata. Gem group numbers correspond to order of samples in aggregation.csv
classification.vec <- as.numeric(gsub(".*-","", (colnames(x = wt_atlas_rna))))
names(classification.vec) <- colnames(x = wt_atlas_rna)
classification.vec[classification.vec=="1"] <- "E775_1"
classification.vec[classification.vec=="2"] <- "E775_2"
classification.vec[classification.vec=="3"] <- "E775_3"
classification.vec[classification.vec=="4"] <- "E775_4"
classification.vec[classification.vec=="5"] <- "E775_5"
classification.vec[classification.vec=="6"] <- "E85_1"
classification.vec[classification.vec=="7"] <- "E85_2"
classification.vec[classification.vec=="8"] <- "E925_1"
classification.vec[classification.vec=="9"] <- "E925_2"
classification.vec[classification.vec=="10"] <- "E925_3"
classification.vec[classification.vec=="11"] <- "E925_4"
classification.vec[classification.vec=="12"] <- "E105_1"
classification.vec[classification.vec=="13"] <- "E105_2"
classification.vec[classification.vec=="14"] <- "E105_3"
classification.vec[classification.vec=="15"] <- "E105_4"
classification.vec[classification.vec=="16"] <- "E115_1"
classification.vec[classification.vec=="17"] <- "E115_2"

wt_atlas_rna$"gem.group" <- classification.vec
head(wt_atlas_rna@meta.data) # gem groups now are listed as numbers, based on order used for aggregation.csv

# Visualize QC metrics as a violin plot
VlnPlot(wt_atlas_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize QC metrics as a 2d scatter plot
plot1 <- FeatureScatter(wt_atlas_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wt_atlas_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# In selecting cells, strict criteria for percent mito, # features and # UMI are used. Selection 
# is done iteratively, re-running the featurescatter plot for sanity check after each subset.


# 1. Filter cells by percent.mito 
wt_atlas_rna <- subset(wt_atlas_rna, subset = percent.mt < 15)
plot1 <- FeatureScatter(wt_atlas_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wt_atlas_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# 2. Filter cells by nCount_RNA (UMI counts)
wt_atlas_rna <- subset(wt_atlas_rna, subset = nCount_RNA < 80000)
plot1 <- FeatureScatter(wt_atlas_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wt_atlas_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


# 3. Filter cells by nFeature_RNA (measured genes)
wt_atlas_rna <- subset(wt_atlas_rna, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500)
plot1 <- FeatureScatter(wt_atlas_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wt_atlas_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# See how many cells are left
wt_atlas_rna
# An object of class Seurat 
# 27173 features across 69215 samples within 1 assay 
# Active assay: RNA (27173 features, 0 variable features)

# cell cycle regression 
cc.genes <- readLines(con = "/Users/sranade/scRNA-seq/projects/cell_cycle_vignette_files/cell_cycle_ssr.txt")
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:45]
g2m.genes <- cc.genes[46:100]

# assign cell cycle scores
wt_atlas_rna <- CellCycleScoring(wt_atlas_rna, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = wt_atlas_rna@meta.data)

##
saveRDS(wt_atlas_rna, file = "WT_pre-SCT.RDS")


# SCTranform replaces Normalize Data, Find variable genes and ScaleData. Regression can be done here for things
# like percent.mito, gene of interest, etc. I usually only do cell cycle regression, but again, feel free to play around with this.
wt_atlas_rna <- SCTransform(wt_atlas_rna, vars.to.regress = c("S.Score", "G2M.Score"), verbose = TRUE)

# Run PCA and batch correction (using MNN)
wt_atlas_rna <- RunPCA(wt_atlas_rna, verbose = TRUE)
wt_atlas_rna_mnn <- RunFastMNN(object.list = SplitObject(wt_atlas_rna, split.by = "gem.group"))

# Cluster
wt_atlas_rna_mnn <- RunUMAP(wt_atlas_rna_mnn, dims = 1:20, reduction = "mnn", verbose = TRUE)
wt_atlas_rna_mnn <- FindNeighbors(wt_atlas_rna_mnn,reduction = "mnn", dims = 1:20, verbose = TRUE)
wt_atlas_rna_mnn <- FindClusters(wt_atlas_rna_mnn, verbose = TRUE, resolution = 0.8)
p5 <- DimPlot(wt_atlas_rna_mnn, reduction = "umap", label = T)
p6 <- DimPlot(wt_atlas_rna_mnn, reduction = "umap", group.by = "gem.group")
CombinePlots(plots = list(p5, p6))

saveRDS(wt_atlas_rna_mnn, file = "/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119/wt_atlas_rna_unannotated.RDS")
wt_atlas_rna_mnn.markers <- FindAllMarkers(wt_atlas_rna_mnn, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.6)
write.csv(wt_atlas_rna_mnn.markers, file = "/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119/wt_atlas_rna_unannotated_all-markers.csv")



# full interation of qc
VlnPlot(wt_atlas_rna_mnn, c("nFeature_RNA"), pt.size = 0.1)
table(wt_atlas_rna_mnn@active.ident,wt_atlas_rna_mnn$gem.group)
wt_atlas_rna_mnn.markers <- FindAllMarkers(wt_atlas_rna_mnn, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

wt_atlas_rna_mnn <- SubsetData(wt_atlas_rna_mnn, ident.remove = c(30,27,25,13,26))
wt_atlas_rna_mnn_trim <- subset(wt_atlas_rna_mnn,idents = c(30,27,25,13,26), invert=TRUE)
table(wt_atlas_rna_mnn_trim@active.ident,wt_atlas_rna_mnn_trim$gem.group)
wt_atlas_rna_mnn_trim@active.ident <- droplevels(wt_atlas_rna_mnn_trim@active.ident)

wt_atlas_rna_mnn_trim <- RunUMAP(wt_atlas_rna_mnn_trim, dims = 1:50, reduction = "mnn", verbose = TRUE)
wt_atlas_rna_mnn_trim <- FindNeighbors(wt_atlas_rna_mnn_trim,reduction = "mnn", dims = 1:50, verbose = TRUE)
wt_atlas_rna_mnn_trim <- FindClusters(wt_atlas_rna_mnn_trim, verbose = TRUE, resolution = 2)
DimPlot(wt_atlas_rna_mnn_trim, reduction = "umap", label = T) + NoLegend()

VlnPlot(wt_atlas_rna_mnn_trim, c("nFeature_RNA"), pt.size = 0.1)
## some of the clusters may be doublets...remove
wt_atlas_rna_mnn_trim <- subset(wt_atlas_rna_mnn_trim, idents = c(50,51,52),invert=TRUE)
wt_atlas_rna_mnn_trim@active.ident<-droplevels(wt_atlas_rna_mnn_trim@active.ident)


wt_atlas_rna_mnn_trim <- RunUMAP(wt_atlas_rna_mnn_trim, dims = 1:40, reduction = "mnn", verbose = TRUE)
wt_atlas_rna_mnn_trim <- FindNeighbors(wt_atlas_rna_mnn_trim,reduction = "mnn", dims = 1:40, verbose = TRUE)
wt_atlas_rna_mnn_trim <- FindClusters(wt_atlas_rna_mnn_trim, verbose = TRUE, resolution = 1.5)
DimPlot(wt_atlas_rna_mnn_trim, reduction = "umap", label = T) + NoLegend()


wt_atlas_rna_mnn_trim <- FindClusters(wt_atlas_rna_mnn_trim, verbose = TRUE, resolution = 1.8)
DimPlot(wt_atlas_rna_mnn_trim, reduction = "umap", label = T) + NoLegend()

# poor qc clusters are an issue, raising the lower limit of #genes/cell
wt_atlas_rna_mnn_trim <- subset(wt_atlas_rna_mnn_trim, subset = nFeature_RNA > 3000)


wt_atlas_rna_mnn_trim <- RunUMAP(wt_atlas_rna_mnn_trim, dims = 1:50, reduction = "mnn", verbose = TRUE)
wt_atlas_rna_mnn_trim <- FindNeighbors(wt_atlas_rna_mnn_trim,reduction = "mnn", dims = 1:50, verbose = TRUE)
wt_atlas_rna_mnn_trim <- FindClusters(wt_atlas_rna_mnn_trim, verbose = TRUE, resolution = 1.6)
DimPlot(wt_atlas_rna_mnn_trim, reduction = "umap", label = T) + NoLegend()


saveRDS(wt_atlas_rna_mnn_trim, file = "/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119/wt_atlas_rna_unannotated.RDS")

# Based on marker gene classification, rename idents, rerun graph
new.cluster.ids <- c("Endoderm",
                     "Endoderm",
                     "Neural_Crest",
                     "Cardiomyocyte",
                     "Neural_Crest",
                     "SHF_Progenitor",
                     "SHF_Progenitor",
                     "SHF_Progenitor",
                     "Paraxial_Mesoderm",
                     "Ectoderm",
                     "Epicardium",
                     "Neural_Crest",
                     "Cardiomyocyte",
                     "Endothelium",
                     "Endothelium",
                     "Cardiomyocyte",
                     "SHF_Progenitor",
                     "SHF_Progenitor",
                     "Ectoderm",
                     "Endoderm",
                     "Ectoderm",
                     "SHF_Progenitor",
                     "LPM",
                     "SHF_Progenitor",
                     "SHF_Progenitor",
                     "Neural_Crest",
                     "Ectoderm",
                     "SHF_Progenitor",
                     "Cardiomyocyte",
                     "Ectoderm",
                     "Epicardium",
                     "Endoderm",
                     "Endoderm",
                     "EndMT",
                     "Endoderm",
                     "Endoderm",
                     "Blood",
                     "Endoderm",
                     "Cardiomyocyte",
                     "Blood",
                     "Ectoderm",
                     "LPM",
                     "Endoderm",
                     "Cardiomyocyte",
                     "Endoderm",
                     "Blood",
                     "Neural_Crest")
names(new.cluster.ids) <- levels(wt_atlas_rna_mnn_trim)
unique(new.cluster.ids)
wt_atlas_rna_mnn_trim <- RenameIdents(wt_atlas_rna_mnn_trim, new.cluster.ids)
p7 <- DimPlot(wt_atlas_rna_mnn_trim, reduction = "umap", label = T)


saveRDS(wt_atlas_rna_mnn_trim, file = "/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119/wt_atlas_rna_annotated.RDS")
wt_atlas_rna_mnn_trim.markers <- FindAllMarkers(wt_atlas_rna_mnn_trim, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(wt_atlas_rna_mnn_trim.markers, file = "/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119/wt_atlas_rna_annotated_all-markers.csv")


# upon restart
setwd("/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119")
library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)

wt_atlas_rna_mnn <- readRDS(file = "wt_atlas_rna_annotated.RDS")
DimPlot(wt_atlas_rna_mnn, reduction = "umap", label = T) + NoLegend()
DimPlot(wt_atlas_rna_mnn, reduction = "umap", label = F) + NoLegend()
DimPlot(wt_atlas_rna_mnn, reduction = "umap", label = F) 

table(wt_atlas_rna_mnn$gem.group)
condition.vec.01 <- substr(wt_atlas_rna_mnn$gem.group,1,3)
table(condition.vec.01)
# condition.vec.01
# E10   E11   E77   E85   E92 
# 14228 10021 10483  8537 18261 

wt_atlas_rna_mnn$TimePoint <- as.character(condition.vec.01)
table(wt_atlas_rna_mnn$TimePoint)
# E10   E11   E77   E85   E92 
# 14228 10021 10483  8537 18261 

DimPlot(wt_atlas_rna_mnn, reduction = "umap",group.by = "TimePoint", label = F) + NoLegend()
DimPlot(wt_atlas_rna_mnn, reduction = "umap",group.by = "TimePoint", label = F)


#### plot gem groups by themselves
DimPlot(scRNA, reduction = "umap",group.by = "TimePoint", label = F) + NoLegend()


#######


#gene_list <- grep("^Grh",rownames(scRNA@assays$SCT@counts),value = TRUE)
gene_list2 <- c("Pmp22","Hand1","Mab21l2",
                "Tbx18","Wt1",
                "Actc1","Tnnt2",
                "Hba-x","Klf1",
                "Dlx5","Dlx2","Sox10",
                "Meox1",
                "Osr1","Isl1","Tbx1","Foxd1","Foxc2",
                "Epcam","Cldn6","Spint2",
                "Sox2","Pou3f1","Pou5f1",
                "Pecam1","Cdh5","Plvap",
                "Twist1","Snai1","Postn")
gene_list_supfig1 <- c("Dlx5","Dlx2","Epcam","Cldn6","Sox2","Pou3f1","Hba-x","Hbb-y","Hand1","Mab21l2","Tbx18","Wt1","Actc1","Tnnt2","Pecam1","Cdh5","Twist1","Tbx20")
pdf(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/1.pdf", width = 6, height = 8)
DotPlot(scRNA, features = gene_list2, col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 6) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend()
dev.off()


DoHeatmap(subset(scRNA, downsample=100), features = gene_list2,slot = "scale.data",assay = "SCT", label = F, raster = F, disp.max = 1.5) + scale_fill_gradientn(colors = c("cornflowerblue", "white", "red"))


###################################  
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

###################
gene_list <- grep("^Rfx",rownames(scRNA@assays$SCT@counts),value = TRUE)
DotPlot(scRNA, features = gene_list, col.max = 2.5, col.min = -2.5, cols = c("grey","red"))+coord_flip() + RotatedAxis() + xlab('') + ylab('')


#################
table(scRNA@active.ident)
# Endoderm      Neural_Crest     Cardiomyocyte    SHF_Progenitor Paraxial_Mesoderm          Ectoderm        Epicardium       Endothelium 
# 10968              8060              7860             14795              1982              7707              2695              3499 
# LPM             EndMT             Blood 
# 1894               699              1371

my_order <-c(
  "Ectoderm",
  "SHF_Progenitor",
  "Blood",
  "EndMT",
  "Endothelium",
  "Paraxial_Mesoderm",
  "Endoderm",
  "Neural_Crest",
  "Cardiomyocyte",
  "Epicardium",
  "LPM"
)
scRNA@active.ident <- factor(x = scRNA@active.ident, levels = my_order)


gene_list_fig1 <- c(
  "Tfap2a","Tfap4",
  "Foxa1","Foxa2","Gata4",
  "Tbx20","Mef2c",
  "Erg","Etv2",
  "Atf3",
  "Klf1","Gata1",
  "Hoxa2","Foxf1",
  "Pou5f1","Sox2")  
pdf(file = "/Users/sranade/Dropbox (Gladstone)/Temp_mail_attach_dump/scriptsyouneed/scRNA/pdf/TF_dotplot_new.pdf", width = 4, height = 6)
DotPlot(scRNA, features = gene_list_fig1, col.max = 1.5, col.min = -2.5,dot.scale = 6, cols = c("grey","red"))+coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend()
dev.off()


### new order 220616

my_order <-c(
  "Blood",
  "Cardiomyocyte",
  "Ectoderm",
  "EndMT",
  "Endoderm",
  "Endothelium",
  "Epicardium",
  "LPM",
  "Neural_Crest",
  "Paraxial_Mesoderm",
  "SHF_Progenitor"
)
scRNA@active.ident <- factor(x = scRNA@active.ident, levels = my_order)



# pdf(file = "/Users/sranade/Dropbox (Gladstone)/NewManuscriptFigures_211128/Figure_2/Gipc2.pdf", width = 6, height = 6)
# FeaturePlot(scRNA, c("Gipc2"))
# dev.off()
# 
# pdf(file = "/Users/sranade/Dropbox (Gladstone)/NewManuscriptFigures_211128/Figure_2/Srrm1.pdf", width = 6, height = 6)
# FeaturePlot(scRNA, c("Srrm1"))
# dev.off()

## heart model
VlnPlot(scRNA, c("Gipc2"), pt.size = 0.1) + NoLegend()

VlnPlot(scRNA, c("Srrm1"), pt.size = 0.1) + NoLegend()

VlnPlot(scRNA, c("Sh3bgr"), pt.size = 0.1) + NoLegend()

VlnPlot(scRNA, c("Fam43a"), pt.size = 0.1) + NoLegend()

VlnPlot(scRNA, c("Epha4"), pt.size = 0.1) + NoLegend()

VlnPlot(scRNA, c("Prr30"), pt.size = 0.1) + NoLegend()





table(scRNA@active.ident, scRNA$gem.group)

#####################
# featureplots for specific TFs
FeaturePlot(wt_atlas_rna_mnn, c("Dlx5"))
FeaturePlot(wt_atlas_rna_mnn, c("Sox10"))
FeaturePlot(wt_atlas_rna_mnn, c("Twist1"))
FeaturePlot(wt_atlas_rna_mnn, c("Rbpj"))

#######
# upon restart
library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)

wt_atlas_rna_mnn <- readRDS(file = "wt_atlas_rna_annotated.RDS")
DimPlot(wt_atlas_rna_mnn, reduction = "umap", label = T) + NoLegend()
DimPlot(wt_atlas_rna_mnn, reduction = "umap", label = F) + NoLegend()
DimPlot(wt_atlas_rna_mnn, reduction = "umap", label = F) 

table(wt_atlas_rna_mnn$gem.group)
condition.vec.01 <- substr(wt_atlas_rna_mnn$gem.group,1,3)
table(condition.vec.01)
# condition.vec.01
# E10   E11   E77   E85   E92 
# 14228 10021 10483  8537 18261 

wt_atlas_rna_mnn$TimePoint <- as.character(condition.vec.01)
table(wt_atlas_rna_mnn$TimePoint)
# E10   E11   E77   E85   E92 
# 14228 10021 10483  8537 18261 

DimPlot(wt_atlas_rna_mnn, reduction = "umap",group.by = "TimePoint", label = F) + NoLegend()
DimPlot(wt_atlas_rna_mnn, reduction = "umap",group.by = "TimePoint", label = F)


#### plot gem groups by themselves
DimPlot(scRNA, reduction = "umap",group.by = "TimePoint", label = F) + NoLegend()



