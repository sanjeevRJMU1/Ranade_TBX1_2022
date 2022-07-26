######## Dot plot of AHF DGE for manuscript
#####################
# re-starting on 220103

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)


Tbx1_mesoderm_cnc <- readRDS(file = "Tbx1_mesoderm_cnc_annotated_210830.RDS")
table(Tbx1_mesoderm_cnc@active.ident,Tbx1_mesoderm_cnc$gem.group)

clust_AHF <- SubsetData(Tbx1_mesoderm_cnc, ident.use = "Anterior_Second_Heart_Field")
clust_AHF_WT <- OldWhichCells(clust_AHF, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_AHF_KO <- OldWhichCells(clust_AHF, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_AHF, cells= clust_AHF_WT) <- "WT_E925"
Idents(object = clust_AHF, cells= clust_AHF_KO) <- "KO_E925"
clust_AHF_findmark <- FindMarkers(clust_AHF, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_AHF_findmark,file="/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE_mesoderm_CNC/annotated/E925/clust_AHF_findmark_220127.csv")

#### re-run DAR/DGE
#### after running homer annotate on motif1 (Tbox motif)
AHF_motif1 <- read.csv(file="/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_220110/AHF/homer_annotate/AHF_E95_down_motif1_outputfile.txt", sep = '\t')
# filter rows without a motif
AHF_motif1_filtered <- AHF_motif1[AHF_motif1$X1.GTGHTGACAGAT.BestGuess.Tbx20.T.box..Heart.Tbx20.ChIP.Seq.GSE29636..Homer.0.874..Distance.From.Peak.sequence.strand.conservation.!="",]
write.csv(AHF_motif1_filtered, file = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_220110/AHF/DAR_DGE/AHF_motif1_filtered.csv", row.names = F)

####################################################
### Intersecting DAR with DGE
library(dplyr)
library(data.table)
library(purrr)


# Read in your list of genes that are differentially expressed per cluster
AHF_DGE_list <- read.csv("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE_mesoderm_CNC/annotated/E925/DAR_DGE_AHF_E95/AHF_E95_DGE_220127.txt", header = F)
AHF_DGE_list<-dplyr::pull(AHF_DGE_list)


# subset where "Genes" is the column in the df that contains the gene names

AHF_Tbx1_DAR_DGE_match <- subset(AHF_motif1_filtered, AHF_motif1_filtered$Gene.Name %in% AHF_DGE_list)
AHF_Tbx1_DAR_DGE_match <- AHF_Tbx1_DAR_DGE_match[,c(2,3,4,6,5,9,10,16)]

# sort by genes with multiple peaks
AHF_Tbx1_DAR_DGE_match_sorted <- as.data.frame(sort(table(AHF_Tbx1_DAR_DGE_match$Gene.Name)))
AHF_Tbx1_DAR_DGE_match_sorted<-AHF_Tbx1_DAR_DGE_match_sorted[!(AHF_Tbx1_DAR_DGE_match_sorted$Freq==0),]
AHF_Tbx1_DAR_DGE_match_sorted <- map_dfr(AHF_Tbx1_DAR_DGE_match_sorted, rev)


## write
write.csv(AHF_Tbx1_DAR_DGE_match,file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE_mesoderm_CNC/annotated/E925/DAR_DGE_AHF_E95/AHF_DARDGE_220127.csv")



# sort by genes with multiple peaks
AHF_Tbx1_DAR_DGE_match_sorted <- as.data.frame(sort(table(AHF_Tbx1_DAR_DGE_match$Gene.Name)))
AHF_Tbx1_DAR_DGE_match_sorted<-AHF_Tbx1_DAR_DGE_match_sorted[!(AHF_Tbx1_DAR_DGE_match_sorted$Freq==0),]
AHF_Tbx1_DAR_DGE_match_sorted <- map_dfr(AHF_Tbx1_DAR_DGE_match_sorted, rev)



AHF_Tbx1_DAR_DGE_match_sorted_trim<-AHF_Tbx1_DAR_DGE_match_sorted[,1]
write.csv(AHF_Tbx1_DAR_DGE_match_sorted_trim, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE_mesoderm_CNC/annotated/E925/DAR_DGE_AHF_E95/AHF_DARDGE_220127.csv", quote = F, row.names = F, sep = '\t')


##################

## generate venn diagrams of overlaps
library(VennDiagram)
grid.newpage()
venn.plot <- draw.pairwise.venn(area1           = 4362,
                                area2           = 563,
                                cross.area      = 257,
                                category        = c("DAR", "TboxMotifs"),
                                fill            = c("cornflowerblue", "yellow"),
                                lty             = "blank",
                                cex             = 2,
                                cat.cex         = 2,
                                cat.pos         = c(285, 105),
                                cat.dist        = 0.09,
                                cat.just        = list(c(-1, -1), c(1, 1)),
                                ext.pos         = 30,
                                ext.dist        = -0.05,
                                ext.length      = 0.85,
                                ext.line.lwd    = 2,
                                ext.line.lty    = "dashed"
)

##################
AHF_subset <- SubsetData(Tbx1_mesoderm_cnc, ident.use = c("Anterior_Second_Heart_Field"))
AHF_subset_WT <- OldWhichCells(AHF_subset, subset.name = "gem.group", accept.value = c("WT_E925"))
AHF_subset_KO <- OldWhichCells(AHF_subset, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = AHF_subset, cells= AHF_subset_WT) <- "WT_E925"
Idents(object = AHF_subset, cells= AHF_subset_KO) <- "KO_E925"

AHF_subset <- SubsetData(AHF_subset, ident.use = c("KO_E925","WT_E925"))
AHF_subset@active.ident<-droplevels(AHF_subset@active.ident)
table(AHF_subset@active.ident)


my_levels<-c("WT_E925","KO_E925")
AHF_subset@active.ident <- factor(x = AHF_subset@active.ident, levels = my_levels)





c7_dot <- c("Sema3c","Sema6d","Sema6a","Sema3e","Fgf8","Fgf10","Fgf1","Fgf13","Wnt5a","Tgfb1","Aldh1a2","Crabp1","Efnb2","Plxna4","Plxna2","Slit1","Lef1","Hoxb1","Nkx2-6","Isl1","Gata6","Ebf1","Ebf2","Pitx2","Nr2f2","Six1")
DotPlot(AHF_subset,cols = c("grey","red") ,features = c7_dot) + coord_flip()

no_tf_dot <- c("Sema3c","Sema6d","Sema6a","Sema3e","Plxna4","Plxna2","Efnb2","Slit1","Hey1","Nrg1","Pdgfc","Fgf1","Fgf8","Fgf10","Wnt5a","Tgfb1","Aldh1a2","Crabp1")
DotPlot(AHF_subset,cols = c("grey","red") ,features = no_tf_dot) + coord_flip()

up_dot <- c("Pitx2","Nr2f1","Nr2f2","Osr1","Ebf2","Six1")
DotPlot(AHF_subset,cols = c("grey","red") ,features = up_dot) + coord_flip()






new_down_genes<- c("Fgf10",
                   "Fgf1",
                   "Fgf8",
                   "Sema3c",
                   "Sema6d",
                   "Sema6a",
                   "Sema3e",
                   "Plxna2",
                   "Plxna4",
                   "Efnb2",
                   "Nrg1",
                   "Gfra1",
                   "Cdh13",
                   "Halr1",
                   "Wnt5a",
                   "Lef1",
                   "Nrxn1",
                   "Igfbpl1",
                   "Chrdl1",
                   "Srgap1",
                   "Crabp1",
                   "Crabp2",
                   "Aldh1a2",
                   "Pdgfc",
                   "Adamts9",
                   "Edn1",
                   "Bmp7",
                   "Slit1",
                   "Hey1",
                   "Tgfb1")

my_levels<-c("KO_E925","WT_E925")
AHF_subset@active.ident <- factor(x = AHF_subset@active.ident, levels = my_levels)

DotPlot(AHF_subset,cols = c("grey","red") ,features = new_down_genes) + RotatedAxis()

###
# patterning genes
E95_ahf_patterning <- c("Pitx2","Nr2f2","Nr2f1","Cyp26a1","Cyp26c1","Fgf8","Isl1","Mef2c")
DotPlot(AHF_subset,cols = c("grey","red") ,features = E95_ahf_patterning) + RotatedAxis()

my_levels<-c("WT_E925","KO_E925")
AHF_subset@active.ident <- factor(x = AHF_subset@active.ident, levels = my_levels)

VlnPlot(AHF_subset, c("Pitx2","Nr2f2","Osr1","Nr2f1","Meis1","Fgf8","Isl1","Mef2c"), ncol = 4, pt.size = 0.1)

VlnPlot(AHF_subset, c("Pitx2","Nr2f1","Osr1","Nr2f2","Ebf2","Six1"), ncol = 3, pt.size = 0.1)
VlnPlot(AHF_subset, c("Pitx2","Nr2f1","Osr1","Nr2f2","Hoxa1","Hoxb1"), ncol = 3, pt.size = 0.1)


##################
# Tbx1 expression in all cell types over time
FeaturePlot(Tbx1_mesoderm_cnc, c("Tbx1"))
### subset just the WT at 3 time points
wt_subset_e95 <- SubsetData(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E925"))
e95_dot <- DotPlot(wt_subset_e95, features = c("Tbx1"))

wt_subset_e105 <- SubsetData(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E105"))
e105_dot <- DotPlot(wt_subset_e105, features = c("Tbx1"))

wt_subset_e115 <- SubsetData(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E115"))
e115_dot <- DotPlot(wt_subset_e115, features = c("Tbx1"))

# plots
e95_dot
e105_dot
e115_dot

#### 
wt_cells <- OldWhichCells(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E925","WT_E105","WT_E115"))
wt_subset_tbx1_expr <- SubsetData(Tbx1_mesoderm_cnc, cells =wt_cells)
table(wt_subset_tbx1_expr@active.ident)
table(wt_subset_tbx1_expr$gem.group)
e115_dot <- DotPlot(wt_subset_tbx1_expr, features = c("Tbx1"), split.by = "active.ident")
e115_dot



#########
Tbx1_mesoderm_cnc_E115 <- Tbx1_mesoderm_cnc
clust_Cardiomyocyte <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Cardiomyocyte")
table(clust_Cardiomyocyte@active.ident)


clust_Cardiomyocyte_WT <- OldWhichCells(clust_Cardiomyocyte, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Cardiomyocyte_KO <- OldWhichCells(clust_Cardiomyocyte, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Cardiomyocyte, cells= clust_Cardiomyocyte_WT) <- "WT_E115"
Idents(object = clust_Cardiomyocyte, cells= clust_Cardiomyocyte_KO) <- "KO_E115"

clust_Cardiomyocyte <- SubsetData(clust_Cardiomyocyte, ident.use = c("WT_E115","KO_E115"))
clust_Cardiomyocyte@active.ident<-droplevels(clust_Cardiomyocyte@active.ident)
table(clust_Cardiomyocyte@active.ident)
clust_Cardiomyocyte_findmark <- FindMarkers(clust_Cardiomyocyte, ident.1 = "WT_E115", ident.2 = "KO_E115")

cm_e115_dge<- c("Nppa",
                "Vsnl1",
                "Stard10",
                "Nr2f2",
                "Nr2f1",
                "Wnt2",
                "Tbx5",
                "Mybpc1",
                "Tnc",
                "Sema3c",
                "Rspo3",
                "Bmp4",
                "Isl1",
                "Barx1")

DotPlot(clust_Cardiomyocyte,cols = c("grey","red") ,features = cm_e115_dge) + RotatedAxis()
my_order <- c("WT_E115","KO_E115")
clust_Cardiomyocyte@active.ident <- factor(x = clust_Cardiomyocyte@active.ident, levels = my_order)
VlnPlot(clust_Cardiomyocyte, c("Rspo3","Bmp4","Sema3c","Isl1","Tbx5","Nr2f1","Nr2f2","Wnt2"), ncol = 4, pt.size = 0.1)






