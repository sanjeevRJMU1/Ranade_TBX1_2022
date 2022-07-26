### Part 2: E115
### run cellchat individually on both wt and ko objects

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)
library(CellChat)
library(patchwork)

# make cell chat objects from part 1
# seurat objects:
wt_cellchat
table(wt_cellchat@active.ident)
ko_cellchat
table(ko_cellchat@active.ident)

# convert both to cell chat
data.input_wt <- GetAssayData(wt_cellchat, assay = "SCT", slot = "data") # normalized data matrix
labels_wt <- Idents(wt_cellchat)
meta_wt <- data.frame(group = labels_wt, row.names = names(labels_wt)) # create a dataframe of the cell labels

data.input_ko <- GetAssayData(ko_cellchat, assay = "SCT", slot = "data") # normalized data matrix
labels_ko <- Idents(ko_cellchat)
meta_ko <- data.frame(group = labels_ko, row.names = names(labels_ko)) # create a dataframe of the cell labels

cellchat_wt <- createCellChat(object = data.input_wt, meta = meta_wt, group.by = "group")
cellchat_wt <- addMeta(cellchat_wt, meta = meta_wt, meta.name = "group")
cellchat_wt <- setIdent(cellchat_wt, ident.use = "group") # set "labels" as default cell identity
levels(cellchat_wt@idents) # show factor levels of the cell labels

cellchat_ko <- createCellChat(object = data.input_ko, meta = meta_ko, group.by = "group")
cellchat_ko <- addMeta(cellchat_ko, meta = meta_ko, meta.name = "group")
cellchat_ko <- setIdent(cellchat_ko, ident.use = "group") # set "labels" as default cell identity
levels(cellchat_ko@idents) # show factor levels of the cell labels

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

###################################################################################################################
# wt object
# set the used database in the object
cellchat_wt@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_wt <- subsetData(cellchat_wt) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 1) # do parallel
cellchat_wt <- identifyOverExpressedGenes(cellchat_wt)
cellchat_wt <- identifyOverExpressedInteractions(cellchat_wt)
# project gene expression data onto PPI network (optional)
cellchat_wt <- projectData(cellchat_wt, PPI.mouse)
### Part II: Inference of cell-cell communication network
cellchat_wt <- computeCommunProb(cellchat_wt, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_wt <- filterCommunication(cellchat_wt, min.cells = 10)

### if you want, subset out dataframes of the communication networks, esp useful if you want to look a specific subset!!!
df.net_wt <- subsetCommunication(cellchat_wt, sources.use = c("cluster1_wt","cluster2_wt","cluster3_wt",
                                                              "cluster5_wt","cluster10_wt","cluster11_wt"),
                                              targets.use = c("wt_Migratory_NeuralCrest","wt_Craniofacial_NeuralCrest",
                                                              "wt_PA3_Cardiac_NeuralCrest"))
write.csv(df.net_wt, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/cellchat_results/unannotated/df_net_wt_unannotated_3NCsubsets_E115_wt.csv")



#Infer the cell-cell communication at a signaling pathway level
cellchat_wt <- computeCommunProbPathway(cellchat_wt)

#Calculate the aggregated cell-cell communication network
cellchat_wt <- aggregateNet(cellchat_wt)

groupSize_wt <- as.numeric(table(cellchat_wt@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_wt@net$count, vertex.weight = groupSize_wt, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_wt@net$weight, vertex.weight = groupSize_wt, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat_wt@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_wt, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

####Part III: Visualization of cell-cell communication network
cellchat_wt@netP$pathways # 49 total pathways
pathways <- as.vector(cellchat_wt@netP$pathways)
write.csv(pathways, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/cellchat_results/unannotated/pathways_unannotated_3NCsubsets_E115_wt.csv")

levels(cellchat_wt@idents) # show factor levels of the cell labels

vertex.receiver = c(1,2,3) # a numeric vector, based on the order of idents
pathways.show <- c("NOTCH") # go thru this fully!


# 1. Hierarchy plot
netVisual_aggregate(cellchat_wt, signaling = pathways.show,  vertex.receiver = vertex.receiver, top = 0.2)

# 2. Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_wt, signaling = pathways.show, layout = "circle", top = 0.2)

# 3. Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_wt, signaling = pathways.show, layout = "chord")

# 4. Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_wt, signaling = pathways.show, color.heatmap = "Reds")

## 
netAnalysis_contribution(cellchat_wt, signaling = pathways.show)

##
plotGeneExpression(cellchat_wt, signaling = "FGF")


####### Part IV: Systems analysis of cell-cell communication network
# Compute the network centrality scores
pathways.show <- c("EDN") # go thru this fully!
cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_wt, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1_wt <- netAnalysis_signalingRole_scatter(cellchat_wt)
# Signaling role analysis on the cell-cell communication networks of interest
gg2_wt <- netAnalysis_signalingRole_scatter(cellchat_wt, signaling = c("FGF"))
#gg1_wt + gg2_wt

gg1_wt
gg2_wt

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1_wt <- netAnalysis_signalingRole_heatmap(cellchat_wt, pattern = "outgoing")
ht2_wt <- netAnalysis_signalingRole_heatmap(cellchat_wt, pattern = "incoming")
#ht1_wt + ht2_wt

ht1_wt
ht2_wt

###################################################################################################################
# ko object
# set the used database in the object
cellchat_ko@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_ko <- subsetData(cellchat_ko) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 1) # do parallel
cellchat_ko <- identifyOverExpressedGenes(cellchat_ko)
cellchat_ko <- identifyOverExpressedInteractions(cellchat_ko)
# project gene expression data onto PPI network (optional)
cellchat_ko <- projectData(cellchat_ko, PPI.mouse)
### Part II: Inference of cell-cell communication network
cellchat_ko <- computeCommunProb(cellchat_ko, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_ko <- filterCommunication(cellchat_ko, min.cells = 10)

### if you want, subset out dataframes of the communication networks, esp useful if you want to look a specific subset!!!
df.net_ko <- subsetCommunication(cellchat_ko, sources.use = c("cluster1_ko","cluster2_ko","cluster3_ko",
                                                              "cluster5_ko","cluster10_ko","cluster11_ko"),
                                 targets.use = c("ko_Migratory_NeuralCrest","ko_Craniofacial_NeuralCrest",
                                                 "ko_PA3_Cardiac_NeuralCrest"))

#Infer the cell-cell communication at a signaling pathway level
cellchat_ko <- computeCommunProbPathway(cellchat_ko)

#Calculate the aggregated cell-cell communication network
cellchat_ko <- aggregateNet(cellchat_ko)

groupSize_ko <- as.numeric(table(cellchat_ko@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_ko@net$count, vertex.weight = groupSize_ko, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_ko@net$weight, vertex.weight = groupSize_ko, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat_ko@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_ko, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

####Part III: Visualization of cell-cell communication network
cellchat_ko@netP$pathways # 53 total pathways
pathways <- as.vector(cellchat_ko@netP$pathways)
write.csv(pathways, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/cellchat_results/unannotated/pathways_unannotated_3NCclusters_E115_ko.csv")

levels(cellchat_ko@idents) # show factor levels of the cell labels

vertex.receiver = c(1,2,3) # a numeric vector, based on the order of idents




# 1. Hierarchy plot
pathways.show <- c("SEMA6") # go thru this fully!
netVisual_aggregate(cellchat_ko, signaling = pathways.show,  vertex.receiver = vertex.receiver, top = 0.2)

# 2. Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_ko, signaling = pathways.show, layout = "circle", top = 0.2)

# 3. Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_ko, signaling = pathways.show, layout = "chord")

# 4. Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_ko, signaling = pathways.show, color.heatmap = "Reds")

## 
netAnalysis_contribution(cellchat_ko, signaling = pathways.show)

##
plotGeneExpression(cellchat_ko, signaling = "FGF")


####### Part IV: Systems analysis of cell-cell communication network
# Compute the network centrality scores
cellchat_ko <- netAnalysis_computeCentrality(cellchat_ko, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_ko, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1_ko <- netAnalysis_signalingRole_scatter(cellchat_ko)
# Signaling role analysis on the cell-cell communication networks of interest
gg2_ko <- netAnalysis_signalingRole_scatter(cellchat_ko, signaling = c("FGF"))
#gg1_ko + gg2_ko

gg1_ko
gg2_ko

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1_ko <- netAnalysis_signalingRole_heatmap(cellchat_ko, pattern = "outgoing")
ht2_ko <- netAnalysis_signalingRole_heatmap(cellchat_ko, pattern = "incoming")
#ht1_ko + ht2_ko

ht1_ko
ht2_ko
