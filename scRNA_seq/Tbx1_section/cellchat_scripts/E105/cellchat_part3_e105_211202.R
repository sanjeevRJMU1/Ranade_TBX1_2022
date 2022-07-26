### Part 3
### run cellchat on wt v ko

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)
library(CellChat)
library(patchwork)

### run wt v ko analysis

# merge cellchat objects
object.list <- list(WT = cellchat_wt, KO = cellchat_ko)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
cellchat_merged
# An object of class CellChat created from a merged object with multiple datasets 
# 972 signaling genes.
# 15291 cells.

## Part I: Predict general principles of cell-cell communication
gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = T)
netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat_merged)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
gg[[1]]
gg[[2]]

###
# gg1 <- netAnalysis_signalingChanges_scatter(cellchat_merged, idents.use = "NeuralCrest_Cardiac_2", signaling.exclude = "NeuralCrest_Craniofacial")
# #> Visualizing differential outgoing and incoming signaling changes from NL to LS
# #> The following `from` values were not present in `x`: 0, -1, 1
# #> The following `from` values were not present in `x`: 0, -1
# gg2 <- netAnalysis_signalingChanges_scatter(cellchat_merged, idents.use = "NeuralCrest_Cardiac_2", signaling.exclude = c("NeuralCrest_Craniofacial"))
# #> Visualizing differential outgoing and incoming signaling changes from NL to LS
# #> The following `from` values were not present in `x`: 0, -1, 1
# #> The following `from` values were not present in `x`: 0, -1
# patchwork::wrap_plots(plots = list(gg1,gg2))


####Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
table(cellchat_merged@idents)
netVisual_bubble(cellchat_merged, sources.use = c(2,3,5,7), targets.use = c(1,4,6,8), comparison = c(1,2), angle.x = 45)
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(2,3,5,7), targets.use = c(1,4,6,8),comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(2,3,5,7), targets.use = c(1,4,6,8),comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


gg1 <- rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merged, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

