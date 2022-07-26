### prop test
# scProportionTest - https://github.com/rpolicastro/scProportionTest


library("scProportionTest")

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

Craniofacial_NeuralCrest_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = "Craniofacial_NeuralCrest")
Craniofacial_NeuralCrest_E115_WT <- OldWhichCells(Craniofacial_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
Craniofacial_NeuralCrest_E115_KO <- OldWhichCells(Craniofacial_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = Craniofacial_NeuralCrest_E115, cells= Craniofacial_NeuralCrest_E115_WT) <- "WT_E115"
Idents(object = Craniofacial_NeuralCrest_E115, cells= Craniofacial_NeuralCrest_E115_KO) <- "KO_E115"


table(Tbx1_CNC_E115@active.ident,Tbx1_CNC_E115$gem.group)

#                           KO_E105 KO_E115 KO_E925 WT_E105 WT_E115 WT_E925
# Craniofacial_NeuralCrest     845     172     279     997    1802     131
# Cardiac_NeuralCrest          935     527       1     780    1162       0
# Migratory_NeuralCrest        304     406     226     197     178     149
# PA3_Cardiac_NeuralCrest      124     313       1     135     561       6


Tbx1_CNC@meta.data$annotated_idents <- Idents(Tbx1_CNC)
#
prop_test <- sc_utils(Tbx1_CNC)

prop_test <- permutation_test(
  prop_test, cluster_identity = "annotated_idents",
  sample_1 = "WT_E115", sample_2 = "KO_E115",
  sample_identity = "gem.group"
)
permutation_plot(prop_test)


### run on E11.5 subset only
CNC_E115
CNC_E115 <- SubsetData(CNC_E115_E115, subset.name = "gem.group", accept.value =c("WT_E115","KO_E115"))
table(CNC_E115@active.ident,CNC_E115$gem.group)

CNC_E115@meta.data$annotated_idents <- Idents(CNC_E115)
#
prop_test <- sc_utils(CNC_E115)

prop_test <- permutation_test(
  prop_test, cluster_identity = "annotated_idents",
  sample_1 = "WT_E115", sample_2 = "KO_E115",
  sample_identity = "gem.group"
)
permutation_plot(prop_test)

# same result
