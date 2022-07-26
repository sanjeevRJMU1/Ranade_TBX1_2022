## label each time point as a sample and plot umap with 3 time points

# Break up each cluster into timepoint and add as a metadata column as Clusters3
table(CNC_subset$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 
# 322         178          76         302        1110        1529        1522        1426         516 
# j_E115_WT_1 k_E115_WT_2 l_E115_KO_2 m_E115_KO_3 
# 1796        2509         299         678 

condition.vec.01 <- substr(CNC_subset$Sample,3,9)

unique(condition.vec.01)
#[1] "E105_WT" "E105_KO" "E115_KO" "E115_WT" "E925_WT" "E925_KO"
table(condition.vec.01)

# E105_KO E105_WT E115_KO E115_WT E925_KO E925_WT 
# 2639    2948    1493    4305     480     398 

CNC_subset$Clusters8 <- as.character(condition.vec.01)
# spot check
table(CNC_subset$Clusters8)

p17 <- plotEmbedding(ArchRProj = wt_atlas, colorBy = "cellColData", name = "Clusters3", embedding = "UMAPHarmony")


p17


#########
### plot separate UMAPs based on time point
table(CNC_subset$Sample,CNC_subset$Clusters6)
## E925_WT
# 1. clone each project 
E925_WT_subset <- CNC_subset
table(E925_WT_subset$Clusters8)
# 2. subset
idxSample_E925_WT <- BiocGenerics::which(CNC_subset$Clusters8 %in% "E925_WT")
cellsSample_E925_WT <- CNC_subset$cellNames[idxSample_E925_WT]
E925_WT_subset<-E925_WT_subset[cellsSample_E925_WT, ]
E925_WT_subset
# numberOfCells(1): 398
# medianTSS(1): 9.394
# medianFrags(1): 51559.5

# confirm it looks good before saving
col_CNC <- c("#E60015","#008939")
plotEmbedding(ArchRProj = E925_WT_subset, colorBy = "cellColData",pal = col_CNC,name = "Clusters6", embedding = "UMAPHarmony")


## E925_KO
# 1. clone each project 
E925_KO_subset <- CNC_subset
# 2. subset
idxSample_E925_KO <- BiocGenerics::which(CNC_subset$Clusters8 %in% "E925_KO")
cellsSample_E925_KO <- CNC_subset$cellNames[idxSample_E925_KO]
E925_KO_subset<-E925_KO_subset[cellsSample_E925_KO, ]
E925_KO_subset
# numberOfCells(1): 480
# medianTSS(1): 9.495
# medianFrags(1): 54934.5

# confirm it looks good before saving
table(E925_KO_subset$Clusters8,E925_KO_subset$Clusters6)

col_CNC <- c("#E60015","#008939","#921990")
plotEmbedding(ArchRProj = E925_KO_subset, colorBy = "cellColData",pal = col_CNC,name = "Clusters6", embedding = "UMAPHarmony")


## E105_WT
# 1. clone each project 
E105_WT_subset <- CNC_subset
# 2. subset
idxSample_E105_WT <- BiocGenerics::which(CNC_subset$Clusters8 %in% "E105_WT")
cellsSample_E105_WT <- CNC_subset$cellNames[idxSample_E105_WT]
E105_WT_subset<-E105_WT_subset[cellsSample_E105_WT, ]
E105_WT_subset
# numberOfCells(1): 2948
# medianTSS(1): 12.8515
# medianFrags(1): 68124

# confirm it looks good before saving
table(E105_WT_subset$Clusters8,E105_WT_subset$Clusters6)

col_CNC <- c("#E60015","#262D6E","#008939","#921990")
plotEmbedding(ArchRProj = E105_WT_subset, colorBy = "cellColData",pal = col_CNC,name = "Clusters6", embedding = "UMAPHarmony")


## E105_KO
# 1. clone each project 
E105_KO_subset <- CNC_subset
# 2. subset
idxSample_E105_KO <- BiocGenerics::which(CNC_subset$Clusters8 %in% "E105_KO")
cellsSample_E105_KO <- CNC_subset$cellNames[idxSample_E105_KO]
E105_KO_subset<-E105_KO_subset[cellsSample_E105_KO, ]
E105_KO_subset
# numberOfCells(1): 2639
# medianTSS(1): 12.874
# medianFrags(1): 64221

# confirm it looks good before saving
table(E105_KO_subset$Clusters8,E105_KO_subset$Clusters6)

col_CNC <- c("#E60015","#262D6E","#008939","#921990")
plotEmbedding(ArchRProj = E105_KO_subset, colorBy = "cellColData",pal = col_CNC,name = "Clusters6", embedding = "UMAPHarmony")


## E115_WT
# 1. clone each project 
E115_WT_subset <- CNC_subset
# 2. subset
idxSample_E115_WT <- BiocGenerics::which(CNC_subset$Clusters8 %in% "E115_WT")
cellsSample_E115_WT <- CNC_subset$cellNames[idxSample_E115_WT]
E115_WT_subset<-E115_WT_subset[cellsSample_E115_WT, ]
E115_WT_subset
# numberOfCells(1): 4305
# medianTSS(1): 11.509
# medianFrags(1): 45771

# confirm it looks good before saving
table(E115_WT_subset$Clusters8,E115_WT_subset$Clusters6)

col_CNC <- c("#E60015","#262D6E","#008939","#921990")
plotEmbedding(ArchRProj = E115_WT_subset, colorBy = "cellColData",pal = col_CNC,name = "Clusters6", embedding = "UMAPHarmony")


## E115_KO
# 1. clone each project 
E115_KO_subset <- CNC_subset
# 2. subset
idxSample_E115_KO <- BiocGenerics::which(CNC_subset$Clusters8 %in% "E115_KO")
cellsSample_E115_KO <- CNC_subset$cellNames[idxSample_E115_KO]
E115_KO_subset<-E115_KO_subset[cellsSample_E115_KO, ]
E115_KO_subset
# numberOfCells(1): 1493
# medianTSS(1): 12.288
# medianFrags(1): 42117

# confirm it looks good before saving
table(E115_KO_subset$Clusters8,E115_KO_subset$Clusters6)

col_CNC <- c("#E60015","#262D6E","#008939","#921990")
plotEmbedding(ArchRProj = E115_KO_subset, colorBy = "cellColData",pal = col_CNC,name = "Clusters6", embedding = "UMAPHarmony")




# save
#p1<-plotEmbedding(ArchRProj = E775_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony", )
p2<-plotEmbedding(ArchRProj = E825_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")
p3<-plotEmbedding(ArchRProj = E925_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")
p4<-plotEmbedding(ArchRProj = E105_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")
p5<-plotEmbedding(ArchRProj = E115_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")

#plotPDF(p1,name = "Plot-UMAP-E775_220610.pdf",ArchRProj = wt_atlas,addDOC = FALSE, width = 5, height = 5)
plotPDF(p2,name = "Plot-UMAP-E825_220610.pdf",ArchRProj = wt_atlas,addDOC = FALSE, width = 5, height = 5)
plotPDF(p3,name = "Plot-UMAP-E925_220610.pdf",ArchRProj = wt_atlas,addDOC = FALSE, width = 5, height = 5)
plotPDF(p4,name = "Plot-UMAP-E105_220610.pdf",ArchRProj = wt_atlas,addDOC = FALSE, width = 5, height = 5)
plotPDF(p5,name = "Plot-UMAP-E115_220610.pdf",ArchRProj = wt_atlas,addDOC = FALSE, width = 5, height = 5)

## problem is that E7.75 has no EndMT cells so it throws off the color palette. set manually for this timepoint.


table(E775_subset$Clusters2)
col_e775 <- c("Blood" = "#D51F26",
              "Cardiomyocyte" ="#272E6A",
              "Ectoderm"="#208A42",
              "Endoderm" = "#F47D2B",
              "Endothelium" ="#FEE500",
              "Epicardium"="#8A9FD1",
              "LPM"="#C06CAB",
              "Neural_Crest"="#D8A767",
              "Paraxial_Mesoderm"="#90D5E4",
              "SHF_Progenitor"="#89C75F")
p1<-plotEmbedding(ArchRProj = E775_subset, colorBy = "cellColData",name = "Clusters2",pal = col_e775, embedding = "UMAPHarmony")
plotPDF(p1,name = "Plot-UMAP-E775_220610.pdf",ArchRProj = wt_atlas,addDOC = FALSE, width = 5, height = 5)










plotPDF(p1,name = "Plot-UMAP-Manual-Annotations.pdf",ArchRProj = wt_atlas,addDOC = FALSE, width = 5, height = 5)
