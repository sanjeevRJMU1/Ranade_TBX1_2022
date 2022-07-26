#### plots to keep
# p# 1. Hierarchy plot
pathways.show <- c("FGF")
netVisual_aggregate(cellchat_wt, signaling = pathways.show,  vertex.receiver = vertex.receiver, top = 0.2)
netVisual_aggregate(cellchat_KO, signaling = pathways.show,  vertex.receiver = vertex.receiver, top = 0.2)


# 4. Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_wt, signaling = pathways.show, color.heatmap = "Reds")
netVisual_heatmap(cellchat_KO, signaling = pathways.show, color.heatmap = "Reds")

## great one
# Compute the network centrality scores
pathways.show <- c("FGF")
cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netAnalysis_signalingRole_network(cellchat_wt, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
gg1_wt <- netAnalysis_signalingRole_scatter(cellchat_wt)
gg2_wt <- netAnalysis_signalingRole_scatter(cellchat_wt, signaling = c("TGFb"))
#gg1_wt + gg2_wt

gg1_wt
gg2_wt

#########
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1_wt <- netAnalysis_signalingRole_heatmap(cellchat_wt, pattern = "outgoing")
ht2_wt <- netAnalysis_signalingRole_heatmap(cellchat_wt, pattern = "incoming")
#ht1_wt + ht2_wt

ht1_wt
ht2_wt


###############################################################
# wt and ko analyzed separately
# p# 1. Hierarchy plot
pathways.show <- c("FGF")
netVisual_aggregate(cellchat_wt, signaling = pathways.show,  vertex.receiver = vertex.receiver, top = 0.2)
netVisual_aggregate(cellchat_KO, signaling = pathways.show,  vertex.receiver = vertex.receiver, top = 0.2)


### DGE


# Compute the network centrality scores
pathways.show <- c("SEMA3")
cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netAnalysis_signalingRole_network(cellchat_wt, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

cellchat_KO <- netAnalysis_computeCentrality(cellchat_KO, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netAnalysis_signalingRole_network(cellchat_KO, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

cellchat_merged@meta$datasets = factor(cellchat_merged@meta$datasets, levels = c("WT", "KO")) # set factor level
plotGeneExpression(cellchat_merged, signaling = "EPHA", split.by = "datasets", colors.ggplot = T)
