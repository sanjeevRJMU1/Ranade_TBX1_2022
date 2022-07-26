### Part 2
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
# An object of class Seurat 
# 54328 features across 1858 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(wt_cellchat@active.ident)
# Anterior_Second_Heart_Field_wt   Cardiopharyngeal_Mesoderm_wt                Neural_Crest_wt 
# 109                            149                            262 
# Cardiopharyngeal_Mesenchyme_wt           Paraxial_Mesoderm_wt          Cranial_Mesenchyme_wt 
# 185                            300                             56 
# Cardiomyocyte_wt                  Epicardium_wt               Smooth_Muscle_wt 
# 270                            265                             87 
# pSHF_wt 
# 175 

KO_cellchat
# An object of class Seurat 
# 54328 features across 2500 samples within 2 assays 
# Active assay: SCT (26701 features, 0 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: mnn, umap

table(KO_cellchat@active.ident)
# Anterior_Second_Heart_Field_KO   Cardiopharyngeal_Mesoderm_KO                Neural_Crest_KO 
# 146                            169                            482 
# Cardiopharyngeal_Mesenchyme_KO           Paraxial_Mesoderm_KO          Cranial_Mesenchyme_KO 
# 236                            284                             83 
# Cardiomyocyte_KO                  Epicardium_KO               Smooth_Muscle_KO 
# 294                            406                            101 
# pSHF_KO 
# 299 

# convert both to cell chat
data.input_wt <- GetAssayData(wt_cellchat, assay = "SCT", slot = "data") # normalized data matrix
labels_wt <- Idents(wt_cellchat)
meta_wt <- data.frame(group = labels_wt, row.names = names(labels_wt)) # create a dataframe of the cell labels

data.input_KO <- GetAssayData(KO_cellchat, assay = "SCT", slot = "data") # normalized data matrix
labels_KO <- Idents(KO_cellchat)
meta_KO <- data.frame(group = labels_KO, row.names = names(labels_KO)) # create a dataframe of the cell labels

cellchat_wt <- createCellChat(object = data.input_wt, meta = meta_wt, group.by = "group")
cellchat_wt <- addMeta(cellchat_wt, meta = meta_wt, meta.name = "group")
cellchat_wt <- setIdent(cellchat_wt, ident.use = "group") # set "labels" as default cell identity
levels(cellchat_wt@idents) # show factor levels of the cell labels

cellchat_KO <- createCellChat(object = data.input_KO, meta = meta_KO, group.by = "group")
cellchat_KO <- addMeta(cellchat_KO, meta = meta_KO, meta.name = "group")
cellchat_KO <- setIdent(cellchat_KO, ident.use = "group") # set "labels" as default cell identity
levels(cellchat_KO@idents) # show factor levels of the cell labels

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
df.net_wt <- subsetCommunication(cellchat_wt, sources.use = c("Anterior_Second_Heart_Field_wt","Cardiopharyngeal_Mesoderm_wt"), targets.use = c("Neural_Crest_wt"))

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
cellchat_wt@netP$pathways # 50 total pathways
pathways <- as.vector(cellchat_wt@netP$pathways)
write.csv(pathways, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/cellchat_results/E95_cellchat/pathways_wt.csv")

levels(cellchat_wt@idents) # show factor levels of the cell labels

vertex.receiver = c(1,2,3,4) # a numeric vector, based on the order of idents
pathways.show <- c("VEGF") # go thru this fully!


# 1. Hierarchy plot
pathways.show <- c("VEGF")
netVisual_aggregate(cellchat_wt, signaling = pathways.show,  vertex.receiver = vertex.receiver, top = 0.2)

# 2. Circle plot
# pathways.show <- c("SEMA3")
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_wt, signaling = pathways.show, layout = "circle", top = 0.2)

# 3. Chord diagram
# pathways.show <- c("SEMA3")
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat_wt, signaling = pathways.show, layout = "chord")

# 4. Heatmap
pathways.show <- c("TGFb")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_wt, signaling = pathways.show, color.heatmap = "Reds")

## 
netAnalysis_contribution(cellchat_wt, signaling = pathways.show)

##
plotGeneExpression(cellchat_wt, signaling = "FGF")


####### Part IV: Systems analysis of cell-cell communication network
# Compute the network centrality scores
cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_wt, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1_wt <- netAnalysis_signalingRole_scatter(cellchat_wt)
# Signaling role analysis on the cell-cell communication networks of interest
gg2_wt <- netAnalysis_signalingRole_scatter(cellchat_wt, signaling = c("BMP"))
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
cellchat_KO@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat_KO <- subsetData(cellchat_KO) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 1) # do parallel
cellchat_KO <- identifyOverExpressedGenes(cellchat_KO)
cellchat_KO <- identifyOverExpressedInteractions(cellchat_KO)
# project gene expression data onto PPI network (optional)
cellchat_KO <- projectData(cellchat_KO, PPI.mouse)
### Part II: Inference of cell-cell communication network
cellchat_KO <- computeCommunProb(cellchat_KO, raw.use = T, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_KO <- filterCommunication(cellchat_KO, min.cells = 10)

### if you want, subset out dataframes of the communication networks, esp useful if you want to look a specific subset!!!
df.net_KO <- subsetCommunication(cellchat_KO, sources.use = c("Anterior_Second_Heart_Field_KO","Cardiopharyngeal_Mesoderm_KO"), targets.use = c("Neural_Crest_KO"))

#Infer the cell-cell communication at a signaling pathway level
cellchat_KO <- computeCommunProbPathway(cellchat_KO)

#Calculate the aggregated cell-cell communication network
cellchat_KO <- aggregateNet(cellchat_KO)

groupSize_KO <- as.numeric(table(cellchat_KO@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_KO@net$count, vertex.weight = groupSize_KO, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_KO@net$weight, vertex.weight = groupSize_KO, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat_KO@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_KO, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

####Part III: Visualization of cell-cell communication network
cellchat_KO@netP$pathways # 51 total pathways
pathways <- as.vector(cellchat_KO@netP$pathways)
write.csv(pathways, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/cellchat_results/E95_cellchat/pathways_KO.csv")

levels(cellchat_KO@idents) # show factor levels of the cell labels

vertex.receiver = c(1,2,3,4) # a numeric vector, based on the order of idents
pathways.show <- c("FGF") # go thru this fully!


# 1. Hierarchy plot
netVisual_aggregate(cellchat_KO, signaling = pathways.show,  vertex.receiver = vertex.receiver, top = 0.2)

# 2. Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_KO, signaling = pathways.show, layout = "circle", top = 0.2)

# 3. Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_KO, signaling = pathways.show, layout = "chord")

# 4. Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_KO, signaling = pathways.show, color.heatmap = "Reds")

## 
netAnalysis_contribution(cellchat_KO, signaling = pathways.show)

##
plotGeneExpression(cellchat_KO, signaling = "FGF")


####### Part IV: Systems analysis of cell-cell communication network
# Compute the network centrality scores
cellchat_KO <- netAnalysis_computeCentrality(cellchat_KO, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_KO, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1_KO <- netAnalysis_signalingRole_scatter(cellchat_KO)
# Signaling role analysis on the cell-cell communication networks of interest
gg2_KO <- netAnalysis_signalingRole_scatter(cellchat_KO, signaling = c("FGF"))
#gg1_KO + gg2_KO

gg1_KO
gg2_KO

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1_KO <- netAnalysis_signalingRole_heatmap(cellchat_KO, pattern = "outgoing")
ht2_KO <- netAnalysis_signalingRole_heatmap(cellchat_KO, pattern = "incoming")
#ht1_KO + ht2_KO

ht1_KO
ht2_KO
