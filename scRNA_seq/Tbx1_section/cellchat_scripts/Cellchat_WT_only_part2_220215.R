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
library(ggalluvial)

# make cell chat objects from part 1
# seurat objects:
Tbx1_mesoderm_cnc_wt_only
table(Tbx1_mesoderm_cnc_wt_only@active.ident)

# convert both to cell chat
data.input_wt <- GetAssayData(Tbx1_mesoderm_cnc_wt_only, assay = "SCT", slot = "data") # normalized data matrix
labels_wt <- Idents(Tbx1_mesoderm_cnc_wt_only)
meta_wt <- data.frame(group = labels_wt, row.names = names(labels_wt)) # create a dataframe of the cell labels

cellchat_wt <- createCellChat(object = data.input_wt, meta = meta_wt, group.by = "group")
cellchat_wt <- addMeta(cellchat_wt, meta = meta_wt, meta.name = "group")
cellchat_wt <- setIdent(cellchat_wt, ident.use = "group") # set "labels" as default cell identity
levels(cellchat_wt@idents) # show factor levels of the cell labels
# [1] "wt_Migratory_NeuralCrest"     "wt_Craniofacial_NeuralCrest"  "wt_PA3_Cardiac_NeuralCrest"  
# [4] "wt_Cardiac_NeuralCrest"       "Anterior_Second_Heart_Field"  "Cardiopharyngeal_Mesoderm"   
# [7] "Cardiopharyngeal_Mesenchyme"  "Cranial_Mesenchyme"           "Smooth_Muscle"               
# [10] "Epicardium"                   "Paraxial_Mesoderm"            "Posterior_Second_Heart_Field"
# [13] "Cardiomyocyte" 

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
df.net_wt <- subsetCommunication(cellchat_wt, sources.use = c("Anterior_Second_Heart_Field","Cardiopharyngeal_Mesoderm","Cardiopharyngeal_Mesenchyme","Cranial_Mesenchyme","Smooth_Muscle","Epicardium","Paraxial_Mesoderm","Posterior_Second_Heart_Field","Cardiomyocyte"), targets.use = c("wt_Migratory_NeuralCrest","wt_Craniofacial_NeuralCrest","wt_PA3_Cardiac_NeuralCrest","wt_Cardiac_NeuralCrest"))
write.csv(df.net_wt, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/cellchat_results/df_net_wt_only_220215.csv")

unique_pathways_communic<-unique(df.net_wt$pathway_name)

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
pathways_wt <- as.vector(cellchat_wt@netP$pathways)
write.csv(pathways_wt, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/cellchat_results/pathways_wt_only_220215.csv")

levels(cellchat_wt@idents) # show factor levels of the cell labels

vertex.receiver = c(1,2,3,4) # a numeric vector, based on the order of idents
pathways.show <- c("TGFb") # go thru this fully!


# 1. Hierarchy plot
pathways.show <- c("BMP")
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
ht1_wt + ht2_wt

ht1_wt
ht2_wt

###
library(NMF)
library(ggalluvial)
library(ggalluvial)

selectK(cellchat_wt, pattern = "outgoing")

nPatterns = 3
cellchat_wt <- identifyCommunicationPatterns(cellchat_wt, pattern = "outgoing", k = 3)
netAnalysis_river(cellchat_wt, pattern = c("outgoing"), targets.use = c(1:4), cutoff = 0.4)
cellchat_wt <- identifyCommunicationPatterns(cellchat_wt, pattern = "incoming", k = 5)
net_in<- netAnalysis_river(cellchat_wt, pattern = c("incoming"),targets.use = c(1:4), cutoff = 0.4)
#> Please make sure you have load `library(ggalluvial)` when running this function

netAnalysis_dot(cellchat_wt, pattern = "incoming")
netAnalysis_dot(cellchat_wt, pattern = "outgoing")

###################################################################################################################
# optimize conditions for incoming patterns 
# varying patterns (k), cutoff at 0.4
## pathways used as from df.net (meaning they pass some statistical value for communicating to any one of neural crest cells)

cellchat_wt <- identifyCommunicationPatterns(cellchat_wt, pattern = "incoming", k = 4)
netAnalysis_river(cellchat_wt, pattern = c("incoming"),targets.use = c(1:4), cutoff = 0.4, signaling = unique_pathways_communic)


### final result is that 4 patterns is good. pathways used = those from df_net run earlier
