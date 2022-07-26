## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/E925")
table(Tbx1_CNC@active.ident,Tbx1_CNC$gem.group)

Tbx1_CNC_E925 <- Tbx1_CNC

#0
clust0_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 0)
clust0_E925_WT <- OldWhichCells(clust0_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
clust0_E925_KO <- OldWhichCells(clust0_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust0_E925, cells= clust0_E925_WT) <- "WT_E925"
Idents(object = clust0_E925, cells= clust0_E925_KO) <- "KO_E925"
clust0_E925_findmark <- FindMarkers(clust0_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust0_E925_findmark, file="clust0_E925_findmark.csv")

# #1
# clust1_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 1)
# clust1_E925_WT <- OldWhichCells(clust1_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
# clust1_E925_KO <- OldWhichCells(clust1_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
# Idents(object = clust1_E925, cells= clust1_E925_WT) <- "WT_E925"
# Idents(object = clust1_E925, cells= clust1_E925_KO) <- "KO_E925"
# clust1_E925_findmark <- FindMarkers(clust1_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
# write.csv(clust1_E925_findmark,file="clust1_E925_findmark.csv")

#2
clust2_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 2)
clust2_E925_WT <- OldWhichCells(clust2_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
clust2_E925_KO <- OldWhichCells(clust2_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust2_E925, cells= clust2_E925_WT) <- "WT_E925"
Idents(object = clust2_E925, cells= clust2_E925_KO) <- "KO_E925"
clust2_E925_findmark <- FindMarkers(clust2_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust2_E925_findmark,file="clust2_E925_findmark.csv")

#3
clust3_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 3)
clust3_E925_WT <- OldWhichCells(clust3_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
clust3_E925_KO <- OldWhichCells(clust3_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = clust3_E925, cells= clust3_E925_WT) <- "WT_E925"
Idents(object = clust3_E925, cells= clust3_E925_KO) <- "KO_E925"
clust3_E925_findmark <- FindMarkers(clust3_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
write.csv(clust3_E925_findmark,file="clust3_E925_findmark.csv")

# #4
# clust4_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 4)
# clust4_E925_WT <- OldWhichCells(clust4_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
# clust4_E925_KO <- OldWhichCells(clust4_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
# Idents(object = clust4_E925, cells= clust4_E925_WT) <- "WT_E925"
# Idents(object = clust4_E925, cells= clust4_E925_KO) <- "KO_E925"
# clust4_E925_findmark <- FindMarkers(clust4_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
# write.csv(clust4_E925_findmark,file="clust4_E925_findmark.csv")
# 
# #5
# clust5_E925 <- SubsetData(Tbx1_CNC_E925, ident.use = 5)
# clust5_E925_WT <- OldWhichCells(clust5_E925, subset.name = "gem.group", accept.value = c("WT_E925"))
# clust5_E925_KO <- OldWhichCells(clust5_E925, subset.name = "gem.group", accept.value = c("KO_E925"))
# Idents(object = clust5_E925, cells= clust5_E925_WT) <- "WT_E925"
# Idents(object = clust5_E925, cells= clust5_E925_KO) <- "KO_E925"
# clust5_E925_findmark <- FindMarkers(clust5_E925, ident.1 = "WT_E925", ident.2 = "KO_E925")
# write.csv(clust5_E925_findmark,file="clust5_E925_findmark.csv")



## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/E105")
Tbx1_CNC_E105 <- Tbx1_CNC

#0
clust0_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 0)
clust0_E105_WT <- OldWhichCells(clust0_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust0_E105_KO <- OldWhichCells(clust0_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust0_E105, cells= clust0_E105_WT) <- "WT_E105"
Idents(object = clust0_E105, cells= clust0_E105_KO) <- "KO_E105"
clust0_E105_findmark <- FindMarkers(clust0_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust0_E105_findmark, file="clust0_E105_findmark.csv")

#1
clust1_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 1)
clust1_E105_WT <- OldWhichCells(clust1_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust1_E105_KO <- OldWhichCells(clust1_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust1_E105, cells= clust1_E105_WT) <- "WT_E105"
Idents(object = clust1_E105, cells= clust1_E105_KO) <- "KO_E105"
clust1_E105_findmark <- FindMarkers(clust1_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust1_E105_findmark,file="clust1_E105_findmark.csv")

#2
clust2_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 2)
clust2_E105_WT <- OldWhichCells(clust2_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust2_E105_KO <- OldWhichCells(clust2_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust2_E105, cells= clust2_E105_WT) <- "WT_E105"
Idents(object = clust2_E105, cells= clust2_E105_KO) <- "KO_E105"
clust2_E105_findmark <- FindMarkers(clust2_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust2_E105_findmark,file="clust2_E105_findmark.csv")

#3
clust3_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 3)
clust3_E105_WT <- OldWhichCells(clust3_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust3_E105_KO <- OldWhichCells(clust3_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust3_E105, cells= clust3_E105_WT) <- "WT_E105"
Idents(object = clust3_E105, cells= clust3_E105_KO) <- "KO_E105"
clust3_E105_findmark <- FindMarkers(clust3_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust3_E105_findmark,file="clust3_E105_findmark.csv")

#4
clust4_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 4)
clust4_E105_WT <- OldWhichCells(clust4_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust4_E105_KO <- OldWhichCells(clust4_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust4_E105, cells= clust4_E105_WT) <- "WT_E105"
Idents(object = clust4_E105, cells= clust4_E105_KO) <- "KO_E105"
clust4_E105_findmark <- FindMarkers(clust4_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust4_E105_findmark,file="clust4_E105_findmark.csv")

#5
clust5_E105 <- SubsetData(Tbx1_CNC_E105, ident.use = 5)
clust5_E105_WT <- OldWhichCells(clust5_E105, subset.name = "gem.group", accept.value = c("WT_E105"))
clust5_E105_KO <- OldWhichCells(clust5_E105, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust5_E105, cells= clust5_E105_WT) <- "WT_E105"
Idents(object = clust5_E105, cells= clust5_E105_KO) <- "KO_E105"
clust5_E105_findmark <- FindMarkers(clust5_E105, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust5_E105_findmark,file="clust5_E105_findmark.csv")



## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/E115")
Tbx1_CNC_E115 <- Tbx1_CNC

#0
clust0_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 0)
clust0_E115_WT <- OldWhichCells(clust0_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust0_E115_KO <- OldWhichCells(clust0_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust0_E115, cells= clust0_E115_WT) <- "WT_E115"
Idents(object = clust0_E115, cells= clust0_E115_KO) <- "KO_E115"
clust0_E115_findmark <- FindMarkers(clust0_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust0_E115_findmark, file="clust0_E115_findmark.csv")

#1
clust1_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 1)
clust1_E115_WT <- OldWhichCells(clust1_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust1_E115_KO <- OldWhichCells(clust1_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust1_E115, cells= clust1_E115_WT) <- "WT_E115"
Idents(object = clust1_E115, cells= clust1_E115_KO) <- "KO_E115"
clust1_E115_findmark <- FindMarkers(clust1_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust1_E115_findmark,file="clust1_E115_findmark.csv")

#2
clust2_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 2)
clust2_E115_WT <- OldWhichCells(clust2_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust2_E115_KO <- OldWhichCells(clust2_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust2_E115, cells= clust2_E115_WT) <- "WT_E115"
Idents(object = clust2_E115, cells= clust2_E115_KO) <- "KO_E115"
clust2_E115_findmark <- FindMarkers(clust2_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust2_E115_findmark,file="clust2_E115_findmark.csv")

#3
clust3_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 3)
clust3_E115_WT <- OldWhichCells(clust3_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust3_E115_KO <- OldWhichCells(clust3_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust3_E115, cells= clust3_E115_WT) <- "WT_E115"
Idents(object = clust3_E115, cells= clust3_E115_KO) <- "KO_E115"
clust3_E115_findmark <- FindMarkers(clust3_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust3_E115_findmark,file="clust3_E115_findmark.csv")

#4
clust4_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 4)
clust4_E115_WT <- OldWhichCells(clust4_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust4_E115_KO <- OldWhichCells(clust4_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust4_E115, cells= clust4_E115_WT) <- "WT_E115"
Idents(object = clust4_E115, cells= clust4_E115_KO) <- "KO_E115"
clust4_E115_findmark <- FindMarkers(clust4_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust4_E115_findmark,file="clust4_E115_findmark.csv")

#5
clust5_E115 <- SubsetData(Tbx1_CNC_E115, ident.use = 5)
clust5_E115_WT <- OldWhichCells(clust5_E115, subset.name = "gem.group", accept.value = c("WT_E115"))
clust5_E115_KO <- OldWhichCells(clust5_E115, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust5_E115, cells= clust5_E115_WT) <- "WT_E115"
Idents(object = clust5_E115, cells= clust5_E115_KO) <- "KO_E115"
clust5_E115_findmark <- FindMarkers(clust5_E115, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust5_E115_findmark,file="clust5_E115_findmark.csv")

##########
clust0v2 <- FindMarkers(Tbx1_CNC, ident.1 = 0, ident.2 = 2)
