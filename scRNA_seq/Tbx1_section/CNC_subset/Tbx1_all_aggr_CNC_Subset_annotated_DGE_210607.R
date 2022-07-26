## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/annotated/E925")
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

Tbx1_CNC_E925 <- Tbx1_CNC

#
Craniofacial_NeuralCrest_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = "Craniofacial_NeuralCrest")
Craniofacial_NeuralCrest_E925_WT <- OldWhichCells(Craniofacial_NeuralCrest_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
Craniofacial_NeuralCrest_E925_KO <- OldWhichCells(Craniofacial_NeuralCrest_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = Craniofacial_NeuralCrest_E925, cells= Craniofacial_NeuralCrest_E925_WT) <- "WT_E925"
Idents(object = Craniofacial_NeuralCrest_E925, cells= Craniofacial_NeuralCrest_E925_KO) <- "KO_E925"
Craniofacial_NeuralCrest_E925_findmark <- FindMarkers(Craniofacial_NeuralCrest_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(Craniofacial_NeuralCrest_E925_findmark, file="Craniofacial_NeuralCrest_E925_findmark.csv")

#
Cardiac_NeuralCrest_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = "Cardiac_NeuralCrest")
Cardiac_NeuralCrest_E925_WT <- OldWhichCells(Cardiac_NeuralCrest_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
Cardiac_NeuralCrest_E925_KO <- OldWhichCells(Cardiac_NeuralCrest_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = Cardiac_NeuralCrest_E925, cells= Cardiac_NeuralCrest_E925_WT) <- "WT_E925"
Idents(object = Cardiac_NeuralCrest_E925, cells= Cardiac_NeuralCrest_E925_KO) <- "KO_E925"
Cardiac_NeuralCrest_E925_findmark <- FindMarkers(Cardiac_NeuralCrest_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(Cardiac_NeuralCrest_E925_findmark, file="Cardiac_NeuralCrest_E925_findmark.csv")

#
Migratory_NeuralCrest_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = "Migratory_NeuralCrest")
Migratory_NeuralCrest_E925_WT <- OldWhichCells(Migratory_NeuralCrest_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
Migratory_NeuralCrest_E925_KO <- OldWhichCells(Migratory_NeuralCrest_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = Migratory_NeuralCrest_E925, cells= Migratory_NeuralCrest_E925_WT) <- "WT_E925"
Idents(object = Migratory_NeuralCrest_E925, cells= Migratory_NeuralCrest_E925_KO) <- "KO_E925"
Migratory_NeuralCrest_E925_findmark <- FindMarkers(Migratory_NeuralCrest_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(Migratory_NeuralCrest_E925_findmark, file="Migratory_NeuralCrest_E925_findmark.csv")

#
PA3_Cardiac_NeuralCrest_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = "PA3_Cardiac_NeuralCrest")
PA3_Cardiac_NeuralCrest_E925_WT <- OldWhichCells(PA3_Cardiac_NeuralCrest_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
PA3_Cardiac_NeuralCrest_E925_KO <- OldWhichCells(PA3_Cardiac_NeuralCrest_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = PA3_Cardiac_NeuralCrest_E925, cells= PA3_Cardiac_NeuralCrest_E925_WT) <- "WT_E925"
Idents(object = PA3_Cardiac_NeuralCrest_E925, cells= PA3_Cardiac_NeuralCrest_E925_KO) <- "KO_E925"
PA3_Cardiac_NeuralCrest_E925_findmark <- FindMarkers(PA3_Cardiac_NeuralCrest_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(PA3_Cardiac_NeuralCrest_E925_findmark, file="PA3_Cardiac_NeuralCrest_E925_findmark.csv")

###########
## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/annotated/E105")
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

Tbx1_CNC_E105 <- Tbx1_CNC

#
Craniofacial_NeuralCrest_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = "Craniofacial_NeuralCrest")
Craniofacial_NeuralCrest_E105_WT <- OldWhichCells(Craniofacial_NeuralCrest_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
Craniofacial_NeuralCrest_E105_KO <- OldWhichCells(Craniofacial_NeuralCrest_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = Craniofacial_NeuralCrest_E105, cells= Craniofacial_NeuralCrest_E105_WT) <- "WT_E105"
Idents(object = Craniofacial_NeuralCrest_E105, cells= Craniofacial_NeuralCrest_E105_KO) <- "KO_E105"
Craniofacial_NeuralCrest_E105_findmark <- FindMarkers(Craniofacial_NeuralCrest_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(Craniofacial_NeuralCrest_E105_findmark, file="Craniofacial_NeuralCrest_E105_findmark.csv")

#
Cardiac_NeuralCrest_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = "Cardiac_NeuralCrest")
Cardiac_NeuralCrest_E105_WT <- OldWhichCells(Cardiac_NeuralCrest_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
Cardiac_NeuralCrest_E105_KO <- OldWhichCells(Cardiac_NeuralCrest_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = Cardiac_NeuralCrest_E105, cells= Cardiac_NeuralCrest_E105_WT) <- "WT_E105"
Idents(object = Cardiac_NeuralCrest_E105, cells= Cardiac_NeuralCrest_E105_KO) <- "KO_E105"
Cardiac_NeuralCrest_E105_findmark <- FindMarkers(Cardiac_NeuralCrest_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(Cardiac_NeuralCrest_E105_findmark, file="Cardiac_NeuralCrest_E105_findmark.csv")

#
Migratory_NeuralCrest_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = "Migratory_NeuralCrest")
Migratory_NeuralCrest_E105_WT <- OldWhichCells(Migratory_NeuralCrest_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
Migratory_NeuralCrest_E105_KO <- OldWhichCells(Migratory_NeuralCrest_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = Migratory_NeuralCrest_E105, cells= Migratory_NeuralCrest_E105_WT) <- "WT_E105"
Idents(object = Migratory_NeuralCrest_E105, cells= Migratory_NeuralCrest_E105_KO) <- "KO_E105"
Migratory_NeuralCrest_E105_findmark <- FindMarkers(Migratory_NeuralCrest_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(Migratory_NeuralCrest_E105_findmark, file="Migratory_NeuralCrest_E105_findmark.csv")

#
PA3_Cardiac_NeuralCrest_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = "PA3_Cardiac_NeuralCrest")
PA3_Cardiac_NeuralCrest_E105_WT <- OldWhichCells(PA3_Cardiac_NeuralCrest_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
PA3_Cardiac_NeuralCrest_E105_KO <- OldWhichCells(PA3_Cardiac_NeuralCrest_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = PA3_Cardiac_NeuralCrest_E105, cells= PA3_Cardiac_NeuralCrest_E105_WT) <- "WT_E105"
Idents(object = PA3_Cardiac_NeuralCrest_E105, cells= PA3_Cardiac_NeuralCrest_E105_KO) <- "KO_E105"
PA3_Cardiac_NeuralCrest_E105_findmark <- FindMarkers(PA3_Cardiac_NeuralCrest_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(PA3_Cardiac_NeuralCrest_E105_findmark, file="PA3_Cardiac_NeuralCrest_E105_findmark.csv")

###############
## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/annotated/E115")
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

Tbx1_CNC_E115 <- Tbx1_CNC

#
Craniofacial_NeuralCrest_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = "Craniofacial_NeuralCrest")
Craniofacial_NeuralCrest_E115_WT <- OldWhichCells(Craniofacial_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
Craniofacial_NeuralCrest_E115_KO <- OldWhichCells(Craniofacial_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = Craniofacial_NeuralCrest_E115, cells= Craniofacial_NeuralCrest_E115_WT) <- "WT_E115"
Idents(object = Craniofacial_NeuralCrest_E115, cells= Craniofacial_NeuralCrest_E115_KO) <- "KO_E115"
Craniofacial_NeuralCrest_E115_findmark <- FindMarkers(Craniofacial_NeuralCrest_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(Craniofacial_NeuralCrest_E115_findmark, file="Craniofacial_NeuralCrest_E115_findmark.csv")

#
Cardiac_NeuralCrest_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = "Cardiac_NeuralCrest")
Cardiac_NeuralCrest_E115_WT <- OldWhichCells(Cardiac_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
Cardiac_NeuralCrest_E115_KO <- OldWhichCells(Cardiac_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = Cardiac_NeuralCrest_E115, cells= Cardiac_NeuralCrest_E115_WT) <- "WT_E115"
Idents(object = Cardiac_NeuralCrest_E115, cells= Cardiac_NeuralCrest_E115_KO) <- "KO_E115"
Cardiac_NeuralCrest_E115_findmark <- FindMarkers(Cardiac_NeuralCrest_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(Cardiac_NeuralCrest_E115_findmark, file="Cardiac_NeuralCrest_E115_findmark.csv")

#
Migratory_NeuralCrest_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = "Migratory_NeuralCrest")
Migratory_NeuralCrest_E115_WT <- OldWhichCells(Migratory_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
Migratory_NeuralCrest_E115_KO <- OldWhichCells(Migratory_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = Migratory_NeuralCrest_E115, cells= Migratory_NeuralCrest_E115_WT) <- "WT_E115"
Idents(object = Migratory_NeuralCrest_E115, cells= Migratory_NeuralCrest_E115_KO) <- "KO_E115"
Migratory_NeuralCrest_E115_findmark <- FindMarkers(Migratory_NeuralCrest_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(Migratory_NeuralCrest_E115_findmark, file="Migratory_NeuralCrest_E115_findmark.csv")

#
PA3_Cardiac_NeuralCrest_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = "PA3_Cardiac_NeuralCrest")
PA3_Cardiac_NeuralCrest_E115_WT <- OldWhichCells(PA3_Cardiac_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
PA3_Cardiac_NeuralCrest_E115_KO <- OldWhichCells(PA3_Cardiac_NeuralCrest_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = PA3_Cardiac_NeuralCrest_E115, cells= PA3_Cardiac_NeuralCrest_E115_WT) <- "WT_E115"
Idents(object = PA3_Cardiac_NeuralCrest_E115, cells= PA3_Cardiac_NeuralCrest_E115_KO) <- "KO_E115"
PA3_Cardiac_NeuralCrest_E115_findmark <- FindMarkers(PA3_Cardiac_NeuralCrest_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(PA3_Cardiac_NeuralCrest_E115_findmark, file="PA3_Cardiac_NeuralCrest_E115_findmark.csv")




