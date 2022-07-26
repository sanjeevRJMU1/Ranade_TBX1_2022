## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE/annotated/E925")

Tbx1_mesoderm_cnc_E925 <- Tbx1_mesoderm_cnc
table(Tbx1_mesoderm_cnc_E925@active.ident,Tbx1_mesoderm_cnc_E925$gem.group)



#0
clust_pSHF <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "pSHF")
clust_pSHF_WT <- OldWhichCells(clust_pSHF, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_pSHF_KO <- OldWhichCells(clust_pSHF, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_pSHF, cells= clust_pSHF_WT) <- "WT_E925"
Idents(object = clust_pSHF, cells= clust_pSHF_KO) <- "KO_E925"
clust_pSHF_findmark <- FindMarkers(clust_pSHF, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_pSHF_findmark, file="clust_pSHF_findmark.csv")

#1
clust_ParaxialMesoderm <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "ParaxialMesoderm")
clust_ParaxialMesoderm_WT <- OldWhichCells(clust_ParaxialMesoderm, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_ParaxialMesoderm_KO <- OldWhichCells(clust_ParaxialMesoderm, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_ParaxialMesoderm, cells= clust_ParaxialMesoderm_WT) <- "WT_E925"
Idents(object = clust_ParaxialMesoderm, cells= clust_ParaxialMesoderm_KO) <- "KO_E925"
clust_ParaxialMesoderm_findmark <- FindMarkers(clust_ParaxialMesoderm, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_ParaxialMesoderm_findmark,file="clust_ParaxialMesoderm_findmark.csv")

#2
clust_PharyngealMesoderm  <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "PharyngealMesoderm")
clust_PharyngealMesoderm_WT <- OldWhichCells(clust_PharyngealMesoderm , subset.name = "gem.group", accept.value = c("WT_E925"))
clust_PharyngealMesoderm_KO <- OldWhichCells(clust_PharyngealMesoderm , subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_PharyngealMesoderm , cells= clust_PharyngealMesoderm_WT) <- "WT_E925"
Idents(object = clust_PharyngealMesoderm , cells= clust_PharyngealMesoderm_KO) <- "KO_E925"
clust_PharyngealMesoderm_findmark <- FindMarkers(clust_PharyngealMesoderm , ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_PharyngealMesoderm_findmark,file="clust_PharyngealMesoderm _findmark.csv")

#3
clust_SHF_Mesenchyme <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "SHF_Mesenchyme")
clust_SHF_Mesenchyme_WT <- OldWhichCells(clust_SHF_Mesenchyme, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_SHF_Mesenchyme_KO <- OldWhichCells(clust_SHF_Mesenchyme, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_SHF_Mesenchyme, cells= clust_SHF_Mesenchyme_WT) <- "WT_E925"
Idents(object = clust_SHF_Mesenchyme, cells= clust_SHF_Mesenchyme_KO) <- "KO_E925"
clust_SHF_Mesenchyme_findmark <- FindMarkers(clust_SHF_Mesenchyme, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_SHF_Mesenchyme_findmark,file="clust_SHF_Mesenchyme_findmark.csv")

#4
clust_NeuralCrest <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "NeuralCrest")
clust_NeuralCrest_WT <- OldWhichCells(clust_NeuralCrest, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_NeuralCrest_KO <- OldWhichCells(clust_NeuralCrest, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_NeuralCrest, cells= clust_NeuralCrest_WT) <- "WT_E925"
Idents(object = clust_NeuralCrest, cells= clust_NeuralCrest_KO) <- "KO_E925"
clust_NeuralCrest_findmark <- FindMarkers(clust_NeuralCrest, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_NeuralCrest_findmark,file="clust_NeuralCrest_findmark.csv")

#5
clust_AHF <- SubsetData(Tbx1_mesoderm_cnc_E925, ident.use = "AHF")
clust_AHF_WT <- OldWhichCells(clust_AHF, subset.name = "gem.group", accept.value = c("WT_E925"))
clust_AHF_KO <- OldWhichCells(clust_AHF, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust_AHF, cells= clust_AHF_WT) <- "WT_E925"
Idents(object = clust_AHF, cells= clust_AHF_KO) <- "KO_E925"
clust_AHF_findmark <- FindMarkers(clust_AHF, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust_AHF_findmark,file="clust_AHF_findmark.csv")


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

