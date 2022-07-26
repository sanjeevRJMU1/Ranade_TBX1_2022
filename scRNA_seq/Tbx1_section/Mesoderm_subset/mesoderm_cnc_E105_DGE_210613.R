## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE/E105")

Tbx1_mesoderm_cnc_E105 <- Tbx1_mesoderm_cnc
table(Tbx1_mesoderm_cnc_E105@active.ident,Tbx1_mesoderm_cnc_E105$gem.group)



#0
clust0 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 0)
clust0_WT <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("WT_E105"))
clust0_KO <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust0, cells= clust0_WT) <- "WT_E105"
Idents(object = clust0, cells= clust0_KO) <- "KO_E105"
clust0_findmark <- FindMarkers(clust0, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust0_findmark, file="clust0_findmark.csv")

#1
clust1 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 1)
clust1_WT <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("WT_E105"))
clust1_KO <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust1, cells= clust1_WT) <- "WT_E105"
Idents(object = clust1, cells= clust1_KO) <- "KO_E105"
clust1_findmark <- FindMarkers(clust1, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust1_findmark,file="clust1_findmark.csv")

#2
clust2 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 2)
clust2_WT <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("WT_E105"))
clust2_KO <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust2, cells= clust2_WT) <- "WT_E105"
Idents(object = clust2, cells= clust2_KO) <- "KO_E105"
clust2_findmark <- FindMarkers(clust2, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust2_findmark,file="clust2_findmark.csv")

#3
clust3 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 3)
clust3_WT <- OldWhichCells(clust3, subset.name = "gem.group", accept.value = c("WT_E105"))
clust3_KO <- OldWhichCells(clust3, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust3, cells= clust3_WT) <- "WT_E105"
Idents(object = clust3, cells= clust3_KO) <- "KO_E105"
clust3_findmark <- FindMarkers(clust3, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust3_findmark,file="clust3_findmark.csv")

#4
clust4 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 4)
clust4_WT <- OldWhichCells(clust4, subset.name = "gem.group", accept.value = c("WT_E105"))
clust4_KO <- OldWhichCells(clust4, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust4, cells= clust4_WT) <- "WT_E105"
Idents(object = clust4, cells= clust4_KO) <- "KO_E105"
clust4_findmark <- FindMarkers(clust4, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust4_findmark,file="clust4_findmark.csv")

#5
clust5 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 5)
clust5_WT <- OldWhichCells(clust5, subset.name = "gem.group", accept.value = c("WT_E105"))
clust5_KO <- OldWhichCells(clust5, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust5, cells= clust5_WT) <- "WT_E105"
Idents(object = clust5, cells= clust5_KO) <- "KO_E105"
clust5_findmark <- FindMarkers(clust5, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust5_findmark,file="clust5_findmark.csv")


#6
clust6 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 6)
clust6_WT <- OldWhichCells(clust6, subset.name = "gem.group", accept.value = c("WT_E105"))
clust6_KO <- OldWhichCells(clust6, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust6, cells= clust6_WT) <- "WT_E105"
Idents(object = clust6, cells= clust6_KO) <- "KO_E105"
clust6_findmark <- FindMarkers(clust6, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust6_findmark,file="clust6_findmark.csv")

#7
clust7 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 7)
clust7_WT <- OldWhichCells(clust7, subset.name = "gem.group", accept.value = c("WT_E105"))
clust7_KO <- OldWhichCells(clust7, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust7, cells= clust7_WT) <- "WT_E105"
Idents(object = clust7, cells= clust7_KO) <- "KO_E105"
clust7_findmark <- FindMarkers(clust7, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust7_findmark,file="clust7_findmark.csv")

#8
clust8 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 8)
clust8_WT <- OldWhichCells(clust8, subset.name = "gem.group", accept.value = c("WT_E105"))
clust8_KO <- OldWhichCells(clust8, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust8, cells= clust8_WT) <- "WT_E105"
Idents(object = clust8, cells= clust8_KO) <- "KO_E105"
clust8_findmark <- FindMarkers(clust8, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust8_findmark,file="clust8_findmark.csv")

#9
clust9 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 9)
clust9_WT <- OldWhichCells(clust9, subset.name = "gem.group", accept.value = c("WT_E105"))
clust9_KO <- OldWhichCells(clust9, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust9, cells= clust9_WT) <- "WT_E105"
Idents(object = clust9, cells= clust9_KO) <- "KO_E105"
clust9_findmark <- FindMarkers(clust9, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust9_findmark,file="clust9_findmark.csv")

#10
clust10 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 10)
clust10_WT <- OldWhichCells(clust10, subset.name = "gem.group", accept.value = c("WT_E105"))
clust10_KO <- OldWhichCells(clust10, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust10, cells= clust10_WT) <- "WT_E105"
Idents(object = clust10, cells= clust10_KO) <- "KO_E105"
clust10_findmark <- FindMarkers(clust10, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust10_findmark,file="clust10_findmark.csv")

#11
clust11 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 11)
clust11_WT <- OldWhichCells(clust11, subset.name = "gem.group", accept.value = c("WT_E105"))
clust11_KO <- OldWhichCells(clust11, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust11, cells= clust11_WT) <- "WT_E105"
Idents(object = clust11, cells= clust11_KO) <- "KO_E105"
clust11_findmark <- FindMarkers(clust11, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust11_findmark,file="clust11_findmark.csv")

#12
clust12 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 12)
clust12_WT <- OldWhichCells(clust12, subset.name = "gem.group", accept.value = c("WT_E105"))
clust12_KO <- OldWhichCells(clust12, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust12, cells= clust12_WT) <- "WT_E105"
Idents(object = clust12, cells= clust12_KO) <- "KO_E105"
clust12_findmark <- FindMarkers(clust12, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust12_findmark,file="clust12_findmark.csv")
# 
# #13
# clust13 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 13)
# clust13_WT <- OldWhichCells(clust13, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust13_KO <- OldWhichCells(clust13, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust13, cells= clust13_WT) <- "WT_E105"
# Idents(object = clust13, cells= clust13_KO) <- "KO_E105"
# clust13_findmark <- FindMarkers(clust13, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust13_findmark,file="clust13_findmark.csv")
# 
# #14
# clust14 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 14)
# clust14_WT <- OldWhichCells(clust14, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust14_KO <- OldWhichCells(clust14, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust14, cells= clust14_WT) <- "WT_E105"
# Idents(object = clust14, cells= clust14_KO) <- "KO_E105"
# clust14_findmark <- FindMarkers(clust14, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust14_findmark,file="clust14_findmark.csv")
# 
# #15
# clust15 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 15)
# clust15_WT <- OldWhichCells(clust15, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust15_KO <- OldWhichCells(clust15, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust15, cells= clust15_WT) <- "WT_E105"
# Idents(object = clust15, cells= clust15_KO) <- "KO_E105"
# clust15_findmark <- FindMarkers(clust15, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust15_findmark,file="clust15_findmark.csv")
# 
# #16
# clust16 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 16)
# clust16_WT <- OldWhichCells(clust16, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust16_KO <- OldWhichCells(clust16, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust16, cells= clust16_WT) <- "WT_E105"
# Idents(object = clust16, cells= clust16_KO) <- "KO_E105"
# clust16_findmark <- FindMarkers(clust16, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust16_findmark,file="clust16_findmark.csv")

# #17
# clust17 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 17)
# clust17_WT <- OldWhichCells(clust17, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust17_KO <- OldWhichCells(clust17, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust17, cells= clust17_WT) <- "WT_E105"
# Idents(object = clust17, cells= clust17_KO) <- "KO_E105"
# clust17_findmark <- FindMarkers(clust17, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust17_findmark,file="clust17_findmark.csv")
# 
# #18
# clust18 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 18)
# clust18_WT <- OldWhichCells(clust18, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust18_KO <- OldWhichCells(clust18, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust18, cells= clust18_WT) <- "WT_E105"
# Idents(object = clust18, cells= clust18_KO) <- "KO_E105"
# clust18_findmark <- FindMarkers(clust18, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust18_findmark,file="clust18_findmark.csv")
# 
# #19
# clust19 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 19)
# clust19_WT <- OldWhichCells(clust19, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust19_KO <- OldWhichCells(clust19, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust19, cells= clust19_WT) <- "WT_E105"
# Idents(object = clust19, cells= clust19_KO) <- "KO_E105"
# clust19_findmark <- FindMarkers(clust19, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust19_findmark,file="clust19_findmark.csv")
# 
# #20
# clust20 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 20)
# clust20_WT <- OldWhichCells(clust20, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust20_KO <- OldWhichCells(clust20, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust20, cells= clust20_WT) <- "WT_E105"
# Idents(object = clust20, cells= clust20_KO) <- "KO_E105"
# clust20_findmark <- FindMarkers(clust20, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust20_findmark,file="clust20_findmark.csv")## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE/E105")

Tbx1_mesoderm_cnc_E105 <- Tbx1_mesoderm_cnc
table(Tbx1_mesoderm_cnc_E105@active.ident,Tbx1_mesoderm_cnc_E105$gem.group)



#0
clust0 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 0)
clust0_WT <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("WT_E105"))
clust0_KO <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust0, cells= clust0_WT) <- "WT_E105"
Idents(object = clust0, cells= clust0_KO) <- "KO_E105"
clust0_findmark <- FindMarkers(clust0, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust0_findmark, file="clust0_findmark.csv")

#1
clust1 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 1)
clust1_WT <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("WT_E105"))
clust1_KO <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust1, cells= clust1_WT) <- "WT_E105"
Idents(object = clust1, cells= clust1_KO) <- "KO_E105"
clust1_findmark <- FindMarkers(clust1, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust1_findmark,file="clust1_findmark.csv")

#2
clust2 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 2)
clust2_WT <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("WT_E105"))
clust2_KO <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust2, cells= clust2_WT) <- "WT_E105"
Idents(object = clust2, cells= clust2_KO) <- "KO_E105"
clust2_findmark <- FindMarkers(clust2, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust2_findmark,file="clust2_findmark.csv")

#3
clust3 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 3)
clust3_WT <- OldWhichCells(clust3, subset.name = "gem.group", accept.value = c("WT_E105"))
clust3_KO <- OldWhichCells(clust3, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust3, cells= clust3_WT) <- "WT_E105"
Idents(object = clust3, cells= clust3_KO) <- "KO_E105"
clust3_findmark <- FindMarkers(clust3, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust3_findmark,file="clust3_findmark.csv")

#4
clust4 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 4)
clust4_WT <- OldWhichCells(clust4, subset.name = "gem.group", accept.value = c("WT_E105"))
clust4_KO <- OldWhichCells(clust4, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust4, cells= clust4_WT) <- "WT_E105"
Idents(object = clust4, cells= clust4_KO) <- "KO_E105"
clust4_findmark <- FindMarkers(clust4, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust4_findmark,file="clust4_findmark.csv")

#5
clust5 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 5)
clust5_WT <- OldWhichCells(clust5, subset.name = "gem.group", accept.value = c("WT_E105"))
clust5_KO <- OldWhichCells(clust5, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust5, cells= clust5_WT) <- "WT_E105"
Idents(object = clust5, cells= clust5_KO) <- "KO_E105"
clust5_findmark <- FindMarkers(clust5, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust5_findmark,file="clust5_findmark.csv")


#6
clust6 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 6)
clust6_WT <- OldWhichCells(clust6, subset.name = "gem.group", accept.value = c("WT_E105"))
clust6_KO <- OldWhichCells(clust6, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust6, cells= clust6_WT) <- "WT_E105"
Idents(object = clust6, cells= clust6_KO) <- "KO_E105"
clust6_findmark <- FindMarkers(clust6, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust6_findmark,file="clust6_findmark.csv")

#7
clust7 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 7)
clust7_WT <- OldWhichCells(clust7, subset.name = "gem.group", accept.value = c("WT_E105"))
clust7_KO <- OldWhichCells(clust7, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust7, cells= clust7_WT) <- "WT_E105"
Idents(object = clust7, cells= clust7_KO) <- "KO_E105"
clust7_findmark <- FindMarkers(clust7, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust7_findmark,file="clust7_findmark.csv")

#8
clust8 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 8)
clust8_WT <- OldWhichCells(clust8, subset.name = "gem.group", accept.value = c("WT_E105"))
clust8_KO <- OldWhichCells(clust8, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust8, cells= clust8_WT) <- "WT_E105"
Idents(object = clust8, cells= clust8_KO) <- "KO_E105"
clust8_findmark <- FindMarkers(clust8, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust8_findmark,file="clust8_findmark.csv")

#9
clust9 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 9)
clust9_WT <- OldWhichCells(clust9, subset.name = "gem.group", accept.value = c("WT_E105"))
clust9_KO <- OldWhichCells(clust9, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust9, cells= clust9_WT) <- "WT_E105"
Idents(object = clust9, cells= clust9_KO) <- "KO_E105"
clust9_findmark <- FindMarkers(clust9, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust9_findmark,file="clust9_findmark.csv")

#10
clust10 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 10)
clust10_WT <- OldWhichCells(clust10, subset.name = "gem.group", accept.value = c("WT_E105"))
clust10_KO <- OldWhichCells(clust10, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust10, cells= clust10_WT) <- "WT_E105"
Idents(object = clust10, cells= clust10_KO) <- "KO_E105"
clust10_findmark <- FindMarkers(clust10, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust10_findmark,file="clust10_findmark.csv")

#11
clust11 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 11)
clust11_WT <- OldWhichCells(clust11, subset.name = "gem.group", accept.value = c("WT_E105"))
clust11_KO <- OldWhichCells(clust11, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust11, cells= clust11_WT) <- "WT_E105"
Idents(object = clust11, cells= clust11_KO) <- "KO_E105"
clust11_findmark <- FindMarkers(clust11, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust11_findmark,file="clust11_findmark.csv")

#12
clust12 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 12)
clust12_WT <- OldWhichCells(clust12, subset.name = "gem.group", accept.value = c("WT_E105"))
clust12_KO <- OldWhichCells(clust12, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust12, cells= clust12_WT) <- "WT_E105"
Idents(object = clust12, cells= clust12_KO) <- "KO_E105"
clust12_findmark <- FindMarkers(clust12, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust12_findmark,file="clust12_findmark.csv")
# 
# #13
clust13 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 13)
clust13_WT <- OldWhichCells(clust13, subset.name = "gem.group", accept.value = c("WT_E105"))
clust13_KO <- OldWhichCells(clust13, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust13, cells= clust13_WT) <- "WT_E105"
Idents(object = clust13, cells= clust13_KO) <- "KO_E105"
clust13_findmark <- FindMarkers(clust13, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust13_findmark,file="clust13_findmark.csv")

#14
clust14 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 14)
clust14_WT <- OldWhichCells(clust14, subset.name = "gem.group", accept.value = c("WT_E105"))
clust14_KO <- OldWhichCells(clust14, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust14, cells= clust14_WT) <- "WT_E105"
Idents(object = clust14, cells= clust14_KO) <- "KO_E105"
clust14_findmark <- FindMarkers(clust14, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust14_findmark,file="clust14_findmark.csv")

#15
clust15 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 15)
clust15_WT <- OldWhichCells(clust15, subset.name = "gem.group", accept.value = c("WT_E105"))
clust15_KO <- OldWhichCells(clust15, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust15, cells= clust15_WT) <- "WT_E105"
Idents(object = clust15, cells= clust15_KO) <- "KO_E105"
clust15_findmark <- FindMarkers(clust15, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust15_findmark,file="clust15_findmark.csv")

#16
clust16 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 16)
clust16_WT <- OldWhichCells(clust16, subset.name = "gem.group", accept.value = c("WT_E105"))
clust16_KO <- OldWhichCells(clust16, subset.name = "gem.group", accept.value = c("KO_E105"))
Idents(object = clust16, cells= clust16_WT) <- "WT_E105"
Idents(object = clust16, cells= clust16_KO) <- "KO_E105"
clust16_findmark <- FindMarkers(clust16, ident.1 = "WT_E105", ident.2 = "KO_E105")
write.csv(clust16_findmark,file="clust16_findmark.csv")

# #17
# clust17 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 17)
# clust17_WT <- OldWhichCells(clust17, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust17_KO <- OldWhichCells(clust17, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust17, cells= clust17_WT) <- "WT_E105"
# Idents(object = clust17, cells= clust17_KO) <- "KO_E105"
# clust17_findmark <- FindMarkers(clust17, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust17_findmark,file="clust17_findmark.csv")
# 
# #18
# clust18 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 18)
# clust18_WT <- OldWhichCells(clust18, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust18_KO <- OldWhichCells(clust18, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust18, cells= clust18_WT) <- "WT_E105"
# Idents(object = clust18, cells= clust18_KO) <- "KO_E105"
# clust18_findmark <- FindMarkers(clust18, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust18_findmark,file="clust18_findmark.csv")
# 
# #19
# clust19 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 19)
# clust19_WT <- OldWhichCells(clust19, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust19_KO <- OldWhichCells(clust19, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust19, cells= clust19_WT) <- "WT_E105"
# Idents(object = clust19, cells= clust19_KO) <- "KO_E105"
# clust19_findmark <- FindMarkers(clust19, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust19_findmark,file="clust19_findmark.csv")
# 
# #20
# clust20 <- SubsetData(Tbx1_mesoderm_cnc_E105, ident.use = 20)
# clust20_WT <- OldWhichCells(clust20, subset.name = "gem.group", accept.value = c("WT_E105"))
# clust20_KO <- OldWhichCells(clust20, subset.name = "gem.group", accept.value = c("KO_E105"))
# Idents(object = clust20, cells= clust20_WT) <- "WT_E105"
# Idents(object = clust20, cells= clust20_KO) <- "KO_E105"
# clust20_findmark <- FindMarkers(clust20, ident.1 = "WT_E105", ident.2 = "KO_E105")
# write.csv(clust20_findmark,file="clust20_findmark.csv")

## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE/E115")

Tbx1_mesoderm_cnc_E115 <- Tbx1_mesoderm_cnc
table(Tbx1_mesoderm_cnc_E115@active.ident,Tbx1_mesoderm_cnc_E115$gem.group)



#0
clust0 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 0)
clust0_WT <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("WT_E115"))
clust0_KO <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust0, cells= clust0_WT) <- "WT_E115"
Idents(object = clust0, cells= clust0_KO) <- "KO_E115"
clust0_findmark <- FindMarkers(clust0, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust0_findmark, file="clust0_findmark.csv")

#1
clust1 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 1)
clust1_WT <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("WT_E115"))
clust1_KO <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust1, cells= clust1_WT) <- "WT_E115"
Idents(object = clust1, cells= clust1_KO) <- "KO_E115"
clust1_findmark <- FindMarkers(clust1, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust1_findmark,file="clust1_findmark.csv")

#2
clust2 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 2)
clust2_WT <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("WT_E115"))
clust2_KO <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust2, cells= clust2_WT) <- "WT_E115"
Idents(object = clust2, cells= clust2_KO) <- "KO_E115"
clust2_findmark <- FindMarkers(clust2, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust2_findmark,file="clust2_findmark.csv")

#3
clust3 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 3)
clust3_WT <- OldWhichCells(clust3, subset.name = "gem.group", accept.value = c("WT_E115"))
clust3_KO <- OldWhichCells(clust3, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust3, cells= clust3_WT) <- "WT_E115"
Idents(object = clust3, cells= clust3_KO) <- "KO_E115"
clust3_findmark <- FindMarkers(clust3, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust3_findmark,file="clust3_findmark.csv")

#4
clust4 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 4)
clust4_WT <- OldWhichCells(clust4, subset.name = "gem.group", accept.value = c("WT_E115"))
clust4_KO <- OldWhichCells(clust4, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust4, cells= clust4_WT) <- "WT_E115"
Idents(object = clust4, cells= clust4_KO) <- "KO_E115"
clust4_findmark <- FindMarkers(clust4, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust4_findmark,file="clust4_findmark.csv")

#5
clust5 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 5)
clust5_WT <- OldWhichCells(clust5, subset.name = "gem.group", accept.value = c("WT_E115"))
clust5_KO <- OldWhichCells(clust5, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust5, cells= clust5_WT) <- "WT_E115"
Idents(object = clust5, cells= clust5_KO) <- "KO_E115"
clust5_findmark <- FindMarkers(clust5, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust5_findmark,file="clust5_findmark.csv")


#6
clust6 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 6)
clust6_WT <- OldWhichCells(clust6, subset.name = "gem.group", accept.value = c("WT_E115"))
clust6_KO <- OldWhichCells(clust6, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust6, cells= clust6_WT) <- "WT_E115"
Idents(object = clust6, cells= clust6_KO) <- "KO_E115"
clust6_findmark <- FindMarkers(clust6, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust6_findmark,file="clust6_findmark.csv")

#7
clust7 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 7)
clust7_WT <- OldWhichCells(clust7, subset.name = "gem.group", accept.value = c("WT_E115"))
clust7_KO <- OldWhichCells(clust7, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust7, cells= clust7_WT) <- "WT_E115"
Idents(object = clust7, cells= clust7_KO) <- "KO_E115"
clust7_findmark <- FindMarkers(clust7, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust7_findmark,file="clust7_findmark.csv")

#8
clust8 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 8)
clust8_WT <- OldWhichCells(clust8, subset.name = "gem.group", accept.value = c("WT_E115"))
clust8_KO <- OldWhichCells(clust8, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust8, cells= clust8_WT) <- "WT_E115"
Idents(object = clust8, cells= clust8_KO) <- "KO_E115"
clust8_findmark <- FindMarkers(clust8, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust8_findmark,file="clust8_findmark.csv")

#9
clust9 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 9)
clust9_WT <- OldWhichCells(clust9, subset.name = "gem.group", accept.value = c("WT_E115"))
clust9_KO <- OldWhichCells(clust9, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust9, cells= clust9_WT) <- "WT_E115"
Idents(object = clust9, cells= clust9_KO) <- "KO_E115"
clust9_findmark <- FindMarkers(clust9, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust9_findmark,file="clust9_findmark.csv")

#10
clust10 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 10)
clust10_WT <- OldWhichCells(clust10, subset.name = "gem.group", accept.value = c("WT_E115"))
clust10_KO <- OldWhichCells(clust10, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust10, cells= clust10_WT) <- "WT_E115"
Idents(object = clust10, cells= clust10_KO) <- "KO_E115"
clust10_findmark <- FindMarkers(clust10, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust10_findmark,file="clust10_findmark.csv")

#11
clust11 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 11)
clust11_WT <- OldWhichCells(clust11, subset.name = "gem.group", accept.value = c("WT_E115"))
clust11_KO <- OldWhichCells(clust11, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust11, cells= clust11_WT) <- "WT_E115"
Idents(object = clust11, cells= clust11_KO) <- "KO_E115"
clust11_findmark <- FindMarkers(clust11, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust11_findmark,file="clust11_findmark.csv")

#12
clust12 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 12)
clust12_WT <- OldWhichCells(clust12, subset.name = "gem.group", accept.value = c("WT_E115"))
clust12_KO <- OldWhichCells(clust12, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust12, cells= clust12_WT) <- "WT_E115"
Idents(object = clust12, cells= clust12_KO) <- "KO_E115"
clust12_findmark <- FindMarkers(clust12, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust12_findmark,file="clust12_findmark.csv")
# 
# #13
# clust13 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 13)
# clust13_WT <- OldWhichCells(clust13, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust13_KO <- OldWhichCells(clust13, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust13, cells= clust13_WT) <- "WT_E115"
# Idents(object = clust13, cells= clust13_KO) <- "KO_E115"
# clust13_findmark <- FindMarkers(clust13, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust13_findmark,file="clust13_findmark.csv")
# 
# #14
# clust14 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 14)
# clust14_WT <- OldWhichCells(clust14, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust14_KO <- OldWhichCells(clust14, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust14, cells= clust14_WT) <- "WT_E115"
# Idents(object = clust14, cells= clust14_KO) <- "KO_E115"
# clust14_findmark <- FindMarkers(clust14, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust14_findmark,file="clust14_findmark.csv")
# 
# #15
# clust15 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 15)
# clust15_WT <- OldWhichCells(clust15, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust15_KO <- OldWhichCells(clust15, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust15, cells= clust15_WT) <- "WT_E115"
# Idents(object = clust15, cells= clust15_KO) <- "KO_E115"
# clust15_findmark <- FindMarkers(clust15, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust15_findmark,file="clust15_findmark.csv")
# 
# #16
# clust16 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 16)
# clust16_WT <- OldWhichCells(clust16, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust16_KO <- OldWhichCells(clust16, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust16, cells= clust16_WT) <- "WT_E115"
# Idents(object = clust16, cells= clust16_KO) <- "KO_E115"
# clust16_findmark <- FindMarkers(clust16, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust16_findmark,file="clust16_findmark.csv")

# #17
# clust17 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 17)
# clust17_WT <- OldWhichCells(clust17, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust17_KO <- OldWhichCells(clust17, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust17, cells= clust17_WT) <- "WT_E115"
# Idents(object = clust17, cells= clust17_KO) <- "KO_E115"
# clust17_findmark <- FindMarkers(clust17, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust17_findmark,file="clust17_findmark.csv")
# 
# #18
# clust18 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 18)
# clust18_WT <- OldWhichCells(clust18, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust18_KO <- OldWhichCells(clust18, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust18, cells= clust18_WT) <- "WT_E115"
# Idents(object = clust18, cells= clust18_KO) <- "KO_E115"
# clust18_findmark <- FindMarkers(clust18, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust18_findmark,file="clust18_findmark.csv")
# 
# #19
# clust19 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 19)
# clust19_WT <- OldWhichCells(clust19, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust19_KO <- OldWhichCells(clust19, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust19, cells= clust19_WT) <- "WT_E115"
# Idents(object = clust19, cells= clust19_KO) <- "KO_E115"
# clust19_findmark <- FindMarkers(clust19, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust19_findmark,file="clust19_findmark.csv")
# 
# #20
# clust20 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 20)
# clust20_WT <- OldWhichCells(clust20, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust20_KO <- OldWhichCells(clust20, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust20, cells= clust20_WT) <- "WT_E115"
# Idents(object = clust20, cells= clust20_KO) <- "KO_E115"
# clust20_findmark <- FindMarkers(clust20, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust20_findmark,file="clust20_findmark.csv")## DGE per cluster
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/DGE/E115")

Tbx1_mesoderm_cnc_E115 <- Tbx1_mesoderm_cnc
table(Tbx1_mesoderm_cnc_E115@active.ident,Tbx1_mesoderm_cnc_E115$gem.group)



#0
clust0 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 0)
clust0_WT <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("WT_E115"))
clust0_KO <- OldWhichCells(clust0, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust0, cells= clust0_WT) <- "WT_E115"
Idents(object = clust0, cells= clust0_KO) <- "KO_E115"
clust0_findmark <- FindMarkers(clust0, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust0_findmark, file="clust0_findmark.csv")

#1
clust1 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 1)
clust1_WT <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("WT_E115"))
clust1_KO <- OldWhichCells(clust1, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust1, cells= clust1_WT) <- "WT_E115"
Idents(object = clust1, cells= clust1_KO) <- "KO_E115"
clust1_findmark <- FindMarkers(clust1, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust1_findmark,file="clust1_findmark.csv")

#2
clust2 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 2)
clust2_WT <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("WT_E115"))
clust2_KO <- OldWhichCells(clust2, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust2, cells= clust2_WT) <- "WT_E115"
Idents(object = clust2, cells= clust2_KO) <- "KO_E115"
clust2_findmark <- FindMarkers(clust2, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust2_findmark,file="clust2_findmark.csv")

#3
clust3 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 3)
clust3_WT <- OldWhichCells(clust3, subset.name = "gem.group", accept.value = c("WT_E115"))
clust3_KO <- OldWhichCells(clust3, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust3, cells= clust3_WT) <- "WT_E115"
Idents(object = clust3, cells= clust3_KO) <- "KO_E115"
clust3_findmark <- FindMarkers(clust3, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust3_findmark,file="clust3_findmark.csv")

#4
clust4 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 4)
clust4_WT <- OldWhichCells(clust4, subset.name = "gem.group", accept.value = c("WT_E115"))
clust4_KO <- OldWhichCells(clust4, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust4, cells= clust4_WT) <- "WT_E115"
Idents(object = clust4, cells= clust4_KO) <- "KO_E115"
clust4_findmark <- FindMarkers(clust4, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust4_findmark,file="clust4_findmark.csv")

#5
clust5 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 5)
clust5_WT <- OldWhichCells(clust5, subset.name = "gem.group", accept.value = c("WT_E115"))
clust5_KO <- OldWhichCells(clust5, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust5, cells= clust5_WT) <- "WT_E115"
Idents(object = clust5, cells= clust5_KO) <- "KO_E115"
clust5_findmark <- FindMarkers(clust5, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust5_findmark,file="clust5_findmark.csv")


#6
clust6 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 6)
clust6_WT <- OldWhichCells(clust6, subset.name = "gem.group", accept.value = c("WT_E115"))
clust6_KO <- OldWhichCells(clust6, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust6, cells= clust6_WT) <- "WT_E115"
Idents(object = clust6, cells= clust6_KO) <- "KO_E115"
clust6_findmark <- FindMarkers(clust6, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust6_findmark,file="clust6_findmark.csv")

#7
clust7 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 7)
clust7_WT <- OldWhichCells(clust7, subset.name = "gem.group", accept.value = c("WT_E115"))
clust7_KO <- OldWhichCells(clust7, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust7, cells= clust7_WT) <- "WT_E115"
Idents(object = clust7, cells= clust7_KO) <- "KO_E115"
clust7_findmark <- FindMarkers(clust7, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust7_findmark,file="clust7_findmark.csv")

#8
clust8 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 8)
clust8_WT <- OldWhichCells(clust8, subset.name = "gem.group", accept.value = c("WT_E115"))
clust8_KO <- OldWhichCells(clust8, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust8, cells= clust8_WT) <- "WT_E115"
Idents(object = clust8, cells= clust8_KO) <- "KO_E115"
clust8_findmark <- FindMarkers(clust8, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust8_findmark,file="clust8_findmark.csv")

#9
clust9 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 9)
clust9_WT <- OldWhichCells(clust9, subset.name = "gem.group", accept.value = c("WT_E115"))
clust9_KO <- OldWhichCells(clust9, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust9, cells= clust9_WT) <- "WT_E115"
Idents(object = clust9, cells= clust9_KO) <- "KO_E115"
clust9_findmark <- FindMarkers(clust9, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust9_findmark,file="clust9_findmark.csv")

#10
clust10 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 10)
clust10_WT <- OldWhichCells(clust10, subset.name = "gem.group", accept.value = c("WT_E115"))
clust10_KO <- OldWhichCells(clust10, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust10, cells= clust10_WT) <- "WT_E115"
Idents(object = clust10, cells= clust10_KO) <- "KO_E115"
clust10_findmark <- FindMarkers(clust10, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust10_findmark,file="clust10_findmark.csv")

#11
clust11 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 11)
clust11_WT <- OldWhichCells(clust11, subset.name = "gem.group", accept.value = c("WT_E115"))
clust11_KO <- OldWhichCells(clust11, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust11, cells= clust11_WT) <- "WT_E115"
Idents(object = clust11, cells= clust11_KO) <- "KO_E115"
clust11_findmark <- FindMarkers(clust11, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust11_findmark,file="clust11_findmark.csv")

#12
clust12 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 12)
clust12_WT <- OldWhichCells(clust12, subset.name = "gem.group", accept.value = c("WT_E115"))
clust12_KO <- OldWhichCells(clust12, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust12, cells= clust12_WT) <- "WT_E115"
Idents(object = clust12, cells= clust12_KO) <- "KO_E115"
clust12_findmark <- FindMarkers(clust12, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust12_findmark,file="clust12_findmark.csv")
# 
# #13
clust13 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 13)
clust13_WT <- OldWhichCells(clust13, subset.name = "gem.group", accept.value = c("WT_E115"))
clust13_KO <- OldWhichCells(clust13, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust13, cells= clust13_WT) <- "WT_E115"
Idents(object = clust13, cells= clust13_KO) <- "KO_E115"
clust13_findmark <- FindMarkers(clust13, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust13_findmark,file="clust13_findmark.csv")

#14
clust14 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 14)
clust14_WT <- OldWhichCells(clust14, subset.name = "gem.group", accept.value = c("WT_E115"))
clust14_KO <- OldWhichCells(clust14, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust14, cells= clust14_WT) <- "WT_E115"
Idents(object = clust14, cells= clust14_KO) <- "KO_E115"
clust14_findmark <- FindMarkers(clust14, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust14_findmark,file="clust14_findmark.csv")

#15
clust15 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 15)
clust15_WT <- OldWhichCells(clust15, subset.name = "gem.group", accept.value = c("WT_E115"))
clust15_KO <- OldWhichCells(clust15, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust15, cells= clust15_WT) <- "WT_E115"
Idents(object = clust15, cells= clust15_KO) <- "KO_E115"
clust15_findmark <- FindMarkers(clust15, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust15_findmark,file="clust15_findmark.csv")

#16
clust16 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 16)
clust16_WT <- OldWhichCells(clust16, subset.name = "gem.group", accept.value = c("WT_E115"))
clust16_KO <- OldWhichCells(clust16, subset.name = "gem.group", accept.value = c("KO_E115"))
Idents(object = clust16, cells= clust16_WT) <- "WT_E115"
Idents(object = clust16, cells= clust16_KO) <- "KO_E115"
clust16_findmark <- FindMarkers(clust16, ident.1 = "WT_E115", ident.2 = "KO_E115")
write.csv(clust16_findmark,file="clust16_findmark.csv")

# #17
# clust17 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 17)
# clust17_WT <- OldWhichCells(clust17, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust17_KO <- OldWhichCells(clust17, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust17, cells= clust17_WT) <- "WT_E115"
# Idents(object = clust17, cells= clust17_KO) <- "KO_E115"
# clust17_findmark <- FindMarkers(clust17, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust17_findmark,file="clust17_findmark.csv")
# 
# #18
# clust18 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 18)
# clust18_WT <- OldWhichCells(clust18, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust18_KO <- OldWhichCells(clust18, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust18, cells= clust18_WT) <- "WT_E115"
# Idents(object = clust18, cells= clust18_KO) <- "KO_E115"
# clust18_findmark <- FindMarkers(clust18, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust18_findmark,file="clust18_findmark.csv")
# 
# #19
# clust19 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 19)
# clust19_WT <- OldWhichCells(clust19, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust19_KO <- OldWhichCells(clust19, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust19, cells= clust19_WT) <- "WT_E115"
# Idents(object = clust19, cells= clust19_KO) <- "KO_E115"
# clust19_findmark <- FindMarkers(clust19, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust19_findmark,file="clust19_findmark.csv")
# 
# #20
# clust20 <- SubsetData(Tbx1_mesoderm_cnc_E115, ident.use = 20)
# clust20_WT <- OldWhichCells(clust20, subset.name = "gem.group", accept.value = c("WT_E115"))
# clust20_KO <- OldWhichCells(clust20, subset.name = "gem.group", accept.value = c("KO_E115"))
# Idents(object = clust20, cells= clust20_WT) <- "WT_E115"
# Idents(object = clust20, cells= clust20_KO) <- "KO_E115"
# clust20_findmark <- FindMarkers(clust20, ident.1 = "WT_E115", ident.2 = "KO_E115")
# write.csv(clust20_findmark,file="clust20_findmark.csv")