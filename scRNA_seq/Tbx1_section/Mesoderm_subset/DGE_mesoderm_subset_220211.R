####### re-run DGE on mesoderm subset using latest scRNA-seq object

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

## re-start 08/29
Tbx1_mesoderm_cnc <- readRDS(file = "Tbx1_mesoderm_cnc_annotated_210830.RDS")

##############
Tbx1_mesoderm_cnc_E925 <- Tbx1_mesoderm_cnc
table(Tbx1_mesoderm_cnc_E925@active.ident,Tbx1_mesoderm_cnc_E925$gem.group)


setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE_mesoderm_CNC/annotated/re-run_220211/E925")


#0
clust_Posterior_Second_Heart_Field <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Posterior_Second_Heart_Field")
clust_Posterior_Second_Heart_Field_WT <- OldWhichCells(clust_Posterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Posterior_Second_Heart_Field_KO <- OldWhichCells(clust_Posterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Posterior_Second_Heart_Field, cells= clust_Posterior_Second_Heart_Field_WT) <- "WT_E925"
Idents(object = clust_Posterior_Second_Heart_Field, cells= clust_Posterior_Second_Heart_Field_KO) <- "KO_E925"
clust_Posterior_Second_Heart_Field_findmark <- FindMarkers(clust_Posterior_Second_Heart_Field, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Posterior_Second_Heart_Field_findmark, file="clust_Posterior_Second_Heart_Field_findmark.csv")

#1
clust_Paraxial_Mesoderm <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Paraxial_Mesoderm")
clust_Paraxial_Mesoderm_WT <- OldWhichCells(clust_Paraxial_Mesoderm, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Paraxial_Mesoderm_KO <- OldWhichCells(clust_Paraxial_Mesoderm, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Paraxial_Mesoderm, cells= clust_Paraxial_Mesoderm_WT) <- "WT_E925"
Idents(object = clust_Paraxial_Mesoderm, cells= clust_Paraxial_Mesoderm_KO) <- "KO_E925"
clust_Paraxial_Mesoderm_findmark <- FindMarkers(clust_Paraxial_Mesoderm, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Paraxial_Mesoderm_findmark,file="clust_Paraxial_Mesoderm_findmark.csv")

#2
clust_Cardiopharyngeal_Mesoderm  <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Cardiopharyngeal_Mesoderm")
clust_Cardiopharyngeal_Mesoderm_WT <- OldWhichCells(clust_Cardiopharyngeal_Mesoderm , subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Cardiopharyngeal_Mesoderm_KO <- OldWhichCells(clust_Cardiopharyngeal_Mesoderm , subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Cardiopharyngeal_Mesoderm , cells= clust_Cardiopharyngeal_Mesoderm_WT) <- "WT_E925"
Idents(object = clust_Cardiopharyngeal_Mesoderm , cells= clust_Cardiopharyngeal_Mesoderm_KO) <- "KO_E925"
clust_Cardiopharyngeal_Mesoderm_findmark <- FindMarkers(clust_Cardiopharyngeal_Mesoderm , ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Cardiopharyngeal_Mesoderm_findmark,file="clust_Cardiopharyngeal_Mesoderm _findmark.csv")

#3
clust_Cardiopharyngeal_Mesenchyme <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Cardiopharyngeal_Mesenchyme")
clust_Cardiopharyngeal_Mesenchyme_WT <- OldWhichCells(clust_Cardiopharyngeal_Mesenchyme, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Cardiopharyngeal_Mesenchyme_KO <- OldWhichCells(clust_Cardiopharyngeal_Mesenchyme, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Cardiopharyngeal_Mesenchyme, cells= clust_Cardiopharyngeal_Mesenchyme_WT) <- "WT_E925"
Idents(object = clust_Cardiopharyngeal_Mesenchyme, cells= clust_Cardiopharyngeal_Mesenchyme_KO) <- "KO_E925"
clust_Cardiopharyngeal_Mesenchyme_findmark <- FindMarkers(clust_Cardiopharyngeal_Mesenchyme, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Cardiopharyngeal_Mesenchyme_findmark,file="clust_Cardiopharyngeal_Mesenchyme_findmark.csv")

#4
clust_Cranial_Mesenchyme <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Cranial_Mesenchyme")
clust_Cranial_Mesenchyme_WT <- OldWhichCells(clust_Cranial_Mesenchyme, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Cranial_Mesenchyme_KO <- OldWhichCells(clust_Cranial_Mesenchyme, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Cranial_Mesenchyme, cells= clust_Cranial_Mesenchyme_WT) <- "WT_E925"
Idents(object = clust_Cranial_Mesenchyme, cells= clust_Cranial_Mesenchyme_KO) <- "KO_E925"
clust_Cranial_Mesenchyme_findmark <- FindMarkers(clust_Cranial_Mesenchyme, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Cranial_Mesenchyme_findmark,file="clust_Cranial_Mesenchyme_findmark.csv")

#5
clust_Anterior_Second_Heart_Field <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Anterior_Second_Heart_Field")
clust_Anterior_Second_Heart_Field_WT <- OldWhichCells(clust_Anterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Anterior_Second_Heart_Field_KO <- OldWhichCells(clust_Anterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Anterior_Second_Heart_Field, cells= clust_Anterior_Second_Heart_Field_WT) <- "WT_E925"
Idents(object = clust_Anterior_Second_Heart_Field, cells= clust_Anterior_Second_Heart_Field_KO) <- "KO_E925"
clust_Anterior_Second_Heart_Field_findmark <- FindMarkers(clust_Anterior_Second_Heart_Field, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Anterior_Second_Heart_Field_findmark,file="clust_Anterior_Second_Heart_Field_findmark.csv")


#6
clust_Epicardium <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Epicardium")
clust_Epicardium_WT <- OldWhichCells(clust_Epicardium, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Epicardium_KO <- OldWhichCells(clust_Epicardium, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Epicardium, cells= clust_Epicardium_WT) <- "WT_E925"
Idents(object = clust_Epicardium, cells= clust_Epicardium_KO) <- "KO_E925"
clust_Epicardium_findmark <- FindMarkers(clust_Epicardium, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Epicardium_findmark,file="clust_Epicardium_findmark.csv")

#7
clust_Cardiomyocyte <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Cardiomyocyte")
clust_Cardiomyocyte_WT <- OldWhichCells(clust_Cardiomyocyte, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Cardiomyocyte_KO <- OldWhichCells(clust_Cardiomyocyte, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Cardiomyocyte, cells= clust_Cardiomyocyte_WT) <- "WT_E925"
Idents(object = clust_Cardiomyocyte, cells= clust_Cardiomyocyte_KO) <- "KO_E925"
clust_Cardiomyocyte_findmark <- FindMarkers(clust_Cardiomyocyte, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Cardiomyocyte_findmark,file="clust_Cardiomyocyte_findmark.csv")

#8
clust_Smooth_Muscle <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "Smooth_Muscle")
clust_Smooth_Muscle_WT <- OldWhichCells(clust_Smooth_Muscle, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_Smooth_Muscle_KO <- OldWhichCells(clust_Smooth_Muscle, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_Smooth_Muscle, cells= clust_Smooth_Muscle_WT) <- "WT_E925"
Idents(object = clust_Smooth_Muscle, cells= clust_Smooth_Muscle_KO) <- "KO_E925"
clust_Smooth_Muscle_findmark <- FindMarkers(clust_Smooth_Muscle, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_Smooth_Muscle_findmark,file="clust_Smooth_Muscle_findmark.csv")


#############################################
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE_mesoderm_CNC/annotated/re-run_220211/E105")

Tbx1_mesoderm_cnc_E105 <- Tbx1_mesoderm_cnc
table(Tbx1_mesoderm_cnc_E105@active.ident,Tbx1_mesoderm_cnc_E105$gem.group)



#0
clust_Posterior_Second_Heart_Field <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Posterior_Second_Heart_Field")
clust_Posterior_Second_Heart_Field_WT <- OldWhichCells(clust_Posterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Posterior_Second_Heart_Field_KO <- OldWhichCells(clust_Posterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Posterior_Second_Heart_Field, cells= clust_Posterior_Second_Heart_Field_WT) <- "WT_E105"
Idents(object = clust_Posterior_Second_Heart_Field, cells= clust_Posterior_Second_Heart_Field_KO) <- "KO_E105"
clust_Posterior_Second_Heart_Field_findmark <- FindMarkers(clust_Posterior_Second_Heart_Field, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Posterior_Second_Heart_Field_findmark, file="clust_Posterior_Second_Heart_Field_findmark.csv")

#1
clust_Paraxial_Mesoderm <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Paraxial_Mesoderm")
clust_Paraxial_Mesoderm_WT <- OldWhichCells(clust_Paraxial_Mesoderm, subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Paraxial_Mesoderm_KO <- OldWhichCells(clust_Paraxial_Mesoderm, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Paraxial_Mesoderm, cells= clust_Paraxial_Mesoderm_WT) <- "WT_E105"
Idents(object = clust_Paraxial_Mesoderm, cells= clust_Paraxial_Mesoderm_KO) <- "KO_E105"
clust_Paraxial_Mesoderm_findmark <- FindMarkers(clust_Paraxial_Mesoderm, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Paraxial_Mesoderm_findmark,file="clust_Paraxial_Mesoderm_findmark.csv")

#2
clust_Cardiopharyngeal_Mesoderm  <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Cardiopharyngeal_Mesoderm")
clust_Cardiopharyngeal_Mesoderm_WT <- OldWhichCells(clust_Cardiopharyngeal_Mesoderm , subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Cardiopharyngeal_Mesoderm_KO <- OldWhichCells(clust_Cardiopharyngeal_Mesoderm , subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Cardiopharyngeal_Mesoderm , cells= clust_Cardiopharyngeal_Mesoderm_WT) <- "WT_E105"
Idents(object = clust_Cardiopharyngeal_Mesoderm , cells= clust_Cardiopharyngeal_Mesoderm_KO) <- "KO_E105"
clust_Cardiopharyngeal_Mesoderm_findmark <- FindMarkers(clust_Cardiopharyngeal_Mesoderm , ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Cardiopharyngeal_Mesoderm_findmark,file="clust_Cardiopharyngeal_Mesoderm _findmark.csv")

#3
clust_Cardiopharyngeal_Mesenchyme <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Cardiopharyngeal_Mesenchyme")
clust_Cardiopharyngeal_Mesenchyme_WT <- OldWhichCells(clust_Cardiopharyngeal_Mesenchyme, subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Cardiopharyngeal_Mesenchyme_KO <- OldWhichCells(clust_Cardiopharyngeal_Mesenchyme, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Cardiopharyngeal_Mesenchyme, cells= clust_Cardiopharyngeal_Mesenchyme_WT) <- "WT_E105"
Idents(object = clust_Cardiopharyngeal_Mesenchyme, cells= clust_Cardiopharyngeal_Mesenchyme_KO) <- "KO_E105"
clust_Cardiopharyngeal_Mesenchyme_findmark <- FindMarkers(clust_Cardiopharyngeal_Mesenchyme, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Cardiopharyngeal_Mesenchyme_findmark,file="clust_Cardiopharyngeal_Mesenchyme_findmark.csv")

#4
clust_Cranial_Mesenchyme <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Cranial_Mesenchyme")
clust_Cranial_Mesenchyme_WT <- OldWhichCells(clust_Cranial_Mesenchyme, subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Cranial_Mesenchyme_KO <- OldWhichCells(clust_Cranial_Mesenchyme, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Cranial_Mesenchyme, cells= clust_Cranial_Mesenchyme_WT) <- "WT_E105"
Idents(object = clust_Cranial_Mesenchyme, cells= clust_Cranial_Mesenchyme_KO) <- "KO_E105"
clust_Cranial_Mesenchyme_findmark <- FindMarkers(clust_Cranial_Mesenchyme, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Cranial_Mesenchyme_findmark,file="clust_Cranial_Mesenchyme_findmark.csv")

#5
clust_Anterior_Second_Heart_Field <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Anterior_Second_Heart_Field")
clust_Anterior_Second_Heart_Field_WT <- OldWhichCells(clust_Anterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Anterior_Second_Heart_Field_KO <- OldWhichCells(clust_Anterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Anterior_Second_Heart_Field, cells= clust_Anterior_Second_Heart_Field_WT) <- "WT_E105"
Idents(object = clust_Anterior_Second_Heart_Field, cells= clust_Anterior_Second_Heart_Field_KO) <- "KO_E105"
clust_Anterior_Second_Heart_Field_findmark <- FindMarkers(clust_Anterior_Second_Heart_Field, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Anterior_Second_Heart_Field_findmark,file="clust_Anterior_Second_Heart_Field_findmark.csv")


#6
clust_Epicardium <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Epicardium")
clust_Epicardium_WT <- OldWhichCells(clust_Epicardium, subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Epicardium_KO <- OldWhichCells(clust_Epicardium, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Epicardium, cells= clust_Epicardium_WT) <- "WT_E105"
Idents(object = clust_Epicardium, cells= clust_Epicardium_KO) <- "KO_E105"
clust_Epicardium_findmark <- FindMarkers(clust_Epicardium, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Epicardium_findmark,file="clust_Epicardium_findmark.csv")

#7
clust_Cardiomyocyte <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Cardiomyocyte")
clust_Cardiomyocyte_WT <- OldWhichCells(clust_Cardiomyocyte, subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Cardiomyocyte_KO <- OldWhichCells(clust_Cardiomyocyte, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Cardiomyocyte, cells= clust_Cardiomyocyte_WT) <- "WT_E105"
Idents(object = clust_Cardiomyocyte, cells= clust_Cardiomyocyte_KO) <- "KO_E105"
clust_Cardiomyocyte_findmark <- FindMarkers(clust_Cardiomyocyte, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Cardiomyocyte_findmark,file="clust_Cardiomyocyte_findmark.csv")

#8
clust_Smooth_Muscle <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = "Smooth_Muscle")
clust_Smooth_Muscle_WT <- OldWhichCells(clust_Smooth_Muscle, subset.name = "gem.group", accept.value = c("WT_E105"))
clust_Smooth_Muscle_KO <- OldWhichCells(clust_Smooth_Muscle, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust_Smooth_Muscle, cells= clust_Smooth_Muscle_WT) <- "WT_E105"
Idents(object = clust_Smooth_Muscle, cells= clust_Smooth_Muscle_KO) <- "KO_E105"
clust_Smooth_Muscle_findmark <- FindMarkers(clust_Smooth_Muscle, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust_Smooth_Muscle_findmark,file="clust_Smooth_Muscle_findmark.csv")


###############################3
## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE_mesoderm_CNC/annotated/re-run_220211/E115")

Tbx1_mesoderm_cnc_E115 <- Tbx1_mesoderm_cnc
table(Tbx1_mesoderm_cnc_E115@active.ident,Tbx1_mesoderm_cnc_E115$gem.group)



#0
clust_Posterior_Second_Heart_Field <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Posterior_Second_Heart_Field")
clust_Posterior_Second_Heart_Field_WT <- OldWhichCells(clust_Posterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Posterior_Second_Heart_Field_KO <- OldWhichCells(clust_Posterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Posterior_Second_Heart_Field, cells= clust_Posterior_Second_Heart_Field_WT) <- "WT_E115"
Idents(object = clust_Posterior_Second_Heart_Field, cells= clust_Posterior_Second_Heart_Field_KO) <- "KO_E115"
clust_Posterior_Second_Heart_Field_findmark <- FindMarkers(clust_Posterior_Second_Heart_Field, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Posterior_Second_Heart_Field_findmark, file="clust_Posterior_Second_Heart_Field_findmark.csv")

#1
clust_Paraxial_Mesoderm <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Paraxial_Mesoderm")
clust_Paraxial_Mesoderm_WT <- OldWhichCells(clust_Paraxial_Mesoderm, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Paraxial_Mesoderm_KO <- OldWhichCells(clust_Paraxial_Mesoderm, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Paraxial_Mesoderm, cells= clust_Paraxial_Mesoderm_WT) <- "WT_E115"
Idents(object = clust_Paraxial_Mesoderm, cells= clust_Paraxial_Mesoderm_KO) <- "KO_E115"
clust_Paraxial_Mesoderm_findmark <- FindMarkers(clust_Paraxial_Mesoderm, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Paraxial_Mesoderm_findmark,file="clust_Paraxial_Mesoderm_findmark.csv")

#2
clust_Cardiopharyngeal_Mesoderm  <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Cardiopharyngeal_Mesoderm")
clust_Cardiopharyngeal_Mesoderm_WT <- OldWhichCells(clust_Cardiopharyngeal_Mesoderm , subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Cardiopharyngeal_Mesoderm_KO <- OldWhichCells(clust_Cardiopharyngeal_Mesoderm , subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Cardiopharyngeal_Mesoderm , cells= clust_Cardiopharyngeal_Mesoderm_WT) <- "WT_E115"
Idents(object = clust_Cardiopharyngeal_Mesoderm , cells= clust_Cardiopharyngeal_Mesoderm_KO) <- "KO_E115"
clust_Cardiopharyngeal_Mesoderm_findmark <- FindMarkers(clust_Cardiopharyngeal_Mesoderm , ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Cardiopharyngeal_Mesoderm_findmark,file="clust_Cardiopharyngeal_Mesoderm _findmark.csv")

#3
clust_Cardiopharyngeal_Mesenchyme <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Cardiopharyngeal_Mesenchyme")
clust_Cardiopharyngeal_Mesenchyme_WT <- OldWhichCells(clust_Cardiopharyngeal_Mesenchyme, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Cardiopharyngeal_Mesenchyme_KO <- OldWhichCells(clust_Cardiopharyngeal_Mesenchyme, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Cardiopharyngeal_Mesenchyme, cells= clust_Cardiopharyngeal_Mesenchyme_WT) <- "WT_E115"
Idents(object = clust_Cardiopharyngeal_Mesenchyme, cells= clust_Cardiopharyngeal_Mesenchyme_KO) <- "KO_E115"
clust_Cardiopharyngeal_Mesenchyme_findmark <- FindMarkers(clust_Cardiopharyngeal_Mesenchyme, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Cardiopharyngeal_Mesenchyme_findmark,file="clust_Cardiopharyngeal_Mesenchyme_findmark.csv")

#4
clust_Cranial_Mesenchyme <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Cranial_Mesenchyme")
clust_Cranial_Mesenchyme_WT <- OldWhichCells(clust_Cranial_Mesenchyme, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Cranial_Mesenchyme_KO <- OldWhichCells(clust_Cranial_Mesenchyme, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Cranial_Mesenchyme, cells= clust_Cranial_Mesenchyme_WT) <- "WT_E115"
Idents(object = clust_Cranial_Mesenchyme, cells= clust_Cranial_Mesenchyme_KO) <- "KO_E115"
clust_Cranial_Mesenchyme_findmark <- FindMarkers(clust_Cranial_Mesenchyme, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Cranial_Mesenchyme_findmark,file="clust_Cranial_Mesenchyme_findmark.csv")

#5
clust_Anterior_Second_Heart_Field <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Anterior_Second_Heart_Field")
clust_Anterior_Second_Heart_Field_WT <- OldWhichCells(clust_Anterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Anterior_Second_Heart_Field_KO <- OldWhichCells(clust_Anterior_Second_Heart_Field, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Anterior_Second_Heart_Field, cells= clust_Anterior_Second_Heart_Field_WT) <- "WT_E115"
Idents(object = clust_Anterior_Second_Heart_Field, cells= clust_Anterior_Second_Heart_Field_KO) <- "KO_E115"
clust_Anterior_Second_Heart_Field_findmark <- FindMarkers(clust_Anterior_Second_Heart_Field, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Anterior_Second_Heart_Field_findmark,file="clust_Anterior_Second_Heart_Field_findmark.csv")


#6
clust_Epicardium <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Epicardium")
clust_Epicardium_WT <- OldWhichCells(clust_Epicardium, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Epicardium_KO <- OldWhichCells(clust_Epicardium, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Epicardium, cells= clust_Epicardium_WT) <- "WT_E115"
Idents(object = clust_Epicardium, cells= clust_Epicardium_KO) <- "KO_E115"
clust_Epicardium_findmark <- FindMarkers(clust_Epicardium, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Epicardium_findmark,file="clust_Epicardium_findmark.csv")

#7
clust_Cardiomyocyte <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Cardiomyocyte")
clust_Cardiomyocyte_WT <- OldWhichCells(clust_Cardiomyocyte, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Cardiomyocyte_KO <- OldWhichCells(clust_Cardiomyocyte, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Cardiomyocyte, cells= clust_Cardiomyocyte_WT) <- "WT_E115"
Idents(object = clust_Cardiomyocyte, cells= clust_Cardiomyocyte_KO) <- "KO_E115"
clust_Cardiomyocyte_findmark <- FindMarkers(clust_Cardiomyocyte, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Cardiomyocyte_findmark,file="clust_Cardiomyocyte_findmark.csv")

#8
clust_Smooth_Muscle <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = "Smooth_Muscle")
clust_Smooth_Muscle_WT <- OldWhichCells(clust_Smooth_Muscle, subset.name = "gem.group", accept.value = c("WT_E115"))
clust_Smooth_Muscle_KO <- OldWhichCells(clust_Smooth_Muscle, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust_Smooth_Muscle, cells= clust_Smooth_Muscle_WT) <- "WT_E115"
Idents(object = clust_Smooth_Muscle, cells= clust_Smooth_Muscle_KO) <- "KO_E115"
clust_Smooth_Muscle_findmark <- FindMarkers(clust_Smooth_Muscle, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust_Smooth_Muscle_findmark,file="clust_Smooth_Muscle_findmark.csv")



