## CNC subset re-start


setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)


##
###
Tbx1_CNC <-readRDS(file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_CNC_annotated.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()
table(Tbx1_CNC@active.ident)

####
# Dotplots of DGE per cluster
Tbx1_CNC_E115 <- Tbx1_CNC
Craniofacial_NeuralCrest_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = "PA3_Cardiac_NeuralCrest")
Craniofacial_NeuralCrest_E115_WT <- OldWhichCells(Craniofacial_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
Craniofacial_NeuralCrest_E115_KO <- OldWhichCells(Craniofacial_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = Craniofacial_NeuralCrest_E115, cells= Craniofacial_NeuralCrest_E115_WT) <- "WT_E115"
Idents(object = Craniofacial_NeuralCrest_E115, cells= Craniofacial_NeuralCrest_E115_KO) <- "KO_E115"

Craniofacial_NeuralCrest_E115 <- SubsetData(Craniofacial_NeuralCrest_E115, ident.use = c("KO_E115","WT_E115"))
table(Craniofacial_NeuralCrest_E115@active.ident)
Craniofacial_NeuralCrest_E115@active.ident <- droplevels(Craniofacial_NeuralCrest_E115@active.ident)

my_levels<-c("WT_E115","KO_E115")
Craniofacial_NeuralCrest_E115@active.ident <- factor(x = Craniofacial_NeuralCrest_E115@active.ident, levels = my_levels)



c7_dot <- c("3110099E03Rik","Emx2","Smoc1","Dlx6","Dlx5","Dlx2","Dlx3","Kctd12","Zfp536","Foxd1","Sox9","Pcdh19","Efnb2","Cntfr",
            "Tbx15","Cxcl14","Foxc2","Hes1")
DotPlot(Craniofacial_NeuralCrest_E115, features = c7_dot) + coord_flip()


####
# Dotplots of DGE per cluster
Tbx1_CNC_E115 <- Tbx1_CNC
PA3_Cardiac_NeuralCrest_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = "PA3_Cardiac_NeuralCrest")
PA3_Cardiac_NeuralCrest_E115_WT <- OldWhichCells(PA3_Cardiac_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
PA3_Cardiac_NeuralCrest_E115_KO <- OldWhichCells(PA3_Cardiac_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = PA3_Cardiac_NeuralCrest_E115, cells= PA3_Cardiac_NeuralCrest_E115_WT) <- "WT_E115"
Idents(object = PA3_Cardiac_NeuralCrest_E115, cells= PA3_Cardiac_NeuralCrest_E115_KO) <- "KO_E115"

PA3_Cardiac_NeuralCrest_E115 <- SubsetData(PA3_Cardiac_NeuralCrest_E115, ident.use = c("KO_E115","WT_E115"))
table(PA3_Cardiac_NeuralCrest_E115@active.ident)
PA3_Cardiac_NeuralCrest_E115@active.ident <- droplevels(PA3_Cardiac_NeuralCrest_E115@active.ident)

my_levels<-c("WT_E115","KO_E115")
PA3_Cardiac_NeuralCrest_E115@active.ident <- factor(x = PA3_Cardiac_NeuralCrest_E115@active.ident, levels = my_levels)

#
c7_dot <- c("Bdnf","Irx3","Bnc2","Six1","Six2","Ebf2","Ebf1","Prrx1","Twist1","Barx1","Dlx1","Dlx5","Dlx6","Foxc1","Foxc2","Foxd2","Foxf1","Foxd1","Shox2",
            "Meox1","Etv1","Etv4","Etv5","Hoxb5","Hes1","Hey1","Tfap2b","Hand2","Arid3a")
#DotPlot(PA3_Cardiac_NeuralCrest_E115,cols = c("grey", "red") ,features = c7_dot) + coord_flip()+NoLegend()
DotPlot(PA3_Cardiac_NeuralCrest_E115,cols = c("grey", "red") ,features = c7_dot) + RotatedAxis() + coord_flip()

### only down
c7_dot <- c("Barx1","Dlx5","Dlx6","Foxc1","Foxc2","Foxd2","Foxf1","Foxd1","Shox2",
            "Meox1","Etv1","Etv4","Etv5","Hes1","Hey1","Tfap2b","Hand2","Arid3a")
#DotPlot(PA3_Cardiac_NeuralCrest_E115,cols = c("grey", "red") ,features = c7_dot) + coord_flip()+NoLegend()
DotPlot(PA3_Cardiac_NeuralCrest_E115,cols = c("grey", "red") ,features = c7_dot) + RotatedAxis() + coord_flip()

up_dge <- c("Bdnf","Irx3","Bnc2","Six1","Six2","Ebf2","Ebf1","Prrx1","Twist1")
DotPlot(PA3_Cardiac_NeuralCrest_E115,cols = c("grey", "red") ,features = up_dge) + RotatedAxis() + coord_flip()

VlnPlot(PA3_Cardiac_NeuralCrest_E115, c("Bdnf","Irx3","Six1","Six2","Ebf2","Ebf1","Prrx1","Twist1"), ncol = 4, pt.size = 0.1)



c7_dot_2 <- c("Etv1","Etv4","Etv5","Tfap2b","Hand2","Foxc1","Foxc2","Foxd2","Foxf1","Foxd1","Hoxb5","Hoxa5","Hoxa4","Hoxb4","Hoxc4","Hoxd4","Dlx6","Dlx5","Dlx1")
DotPlot(PA3_Cardiac_NeuralCrest_E115,cols = c("grey", "blue") ,features = c7_dot_2) + coord_flip()

VlnPlot(PA3_Cardiac_NeuralCrest_E115, c("Nrp1", "Rarb","Satb2","Dkk1"), ncol = 2, pt.size = 0.1)


###############
DotPlot(PA3_Cardiac_NeuralCrest_E115, features=c("Hoxa4","Hoxb4","Hoxc4","Hoxd4","Hoxa5","Hoxb5"), cols = c("grey","red"))+coord_flip()
VlnPlot(PA3_Cardiac_NeuralCrest_E115, c("Hoxa4","Hoxb4","Hoxc4","Hoxd4"), ncol = 2, pt.size = 0.1)


