## label each time point as a sample and plot umap with 5 time points

# Break up each cluster into timepoint and add as a metadata column as Clusters3
table(wt_atlas$Sample)
condition.vec.01 <- substr(wt_atlas$Sample,3,6)

unique(condition.vec.01)
#[1] "E105" "E115" "E925" "E825" "E775"
table(condition.vec.01)
# condition.vec.01
# E105  E115  E775  E825  E925 
# 17406 14017  4002  9496 20035 

wt_atlas$Clusters3 <- as.character(condition.vec.01)
# spot check
table(wt_atlas$Clusters3)

p17 <- plotEmbedding(ArchRProj = wt_atlas, colorBy = "cellColData", name = "Clusters3", embedding = "UMAPHarmony")


p17


#########
### plot separate UMAPs based on time point

## E775
# 1. clone each project 
E775_subset <- wt_atlas
table(E775_subset$Clusters3)
# 2. subset
idxSample_E775 <- BiocGenerics::which(wt_atlas$Clusters3 %in% "E775")
cellsSample_E775 <- wt_atlas$cellNames[idxSample_E775]
E775_subset<-E775_subset[cellsSample_E775, ]
E775_subset
# numberOfCells(1): 4002
# medianTSS(1): 13.965
# medianFrags(1): 46608

# confirm it looks good before saving
plotEmbedding(ArchRProj = E775_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")

## E825
# 1. clone each project 
E825_subset <- wt_atlas
table(E825_subset$Clusters3)
# 2. subset
idxSample_E825 <- BiocGenerics::which(wt_atlas$Clusters3 %in% "E825")
cellsSample_E825 <- wt_atlas$cellNames[idxSample_E825]
E825_subset<-E825_subset[cellsSample_E825, ]
E825_subset
# numberOfCells(1): 9496
# medianTSS(1): 11.461
# medianFrags(1): 58265.5

# confirm it looks good before saving
plotEmbedding(ArchRProj = E825_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")


## E925
# 1. clone each project 
E925_subset <- wt_atlas
table(E925_subset$Clusters3)
# 2. subset
idxSample_E925 <- BiocGenerics::which(wt_atlas$Clusters3 %in% "E925")
cellsSample_E925 <- wt_atlas$cellNames[idxSample_E925]
E925_subset<-E925_subset[cellsSample_E925, ]
E925_subset
# numberOfCells(1): 20035
# medianTSS(1): 11.45
# medianFrags(1): 43609

# confirm it looks good before saving
plotEmbedding(ArchRProj = E925_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")


## E105
# 1. clone each project 
E105_subset <- wt_atlas
table(E105_subset$Clusters3)
# 2. subset
idxSample_E105 <- BiocGenerics::which(wt_atlas$Clusters3 %in% "E105")
cellsSample_E105 <- wt_atlas$cellNames[idxSample_E105]
E105_subset<-E105_subset[cellsSample_E105, ]
E105_subset
# numberOfCells(1): 17406
# medianTSS(1): 13.3325
# medianFrags(1): 56436

# confirm it looks good before saving
plotEmbedding(ArchRProj = E105_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")


## E115
# 1. clone each project 
E115_subset <- wt_atlas
table(E115_subset$Clusters3)
# 2. subset
idxSample_E115 <- BiocGenerics::which(wt_atlas$Clusters3 %in% "E115")
cellsSample_E115 <- wt_atlas$cellNames[idxSample_E115]
E115_subset<-E115_subset[cellsSample_E115, ]
E115_subset
# numberOfCells(1): 14017
# medianTSS(1): 12.033
# medianFrags(1): 40109

# confirm it looks good before saving
plotEmbedding(ArchRProj = E115_subset, colorBy = "cellColData",name = "Clusters2", embedding = "UMAPHarmony")

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
