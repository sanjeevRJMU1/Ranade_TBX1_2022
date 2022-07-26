####### bubble plots of each important pathway disrupted in Tbx1 KO
## main point: ID which genes in pathway (Lig or Receptor) driving the differences

# idents: 
# 1 = migratory prog
# 2 = craniofacial
# 3 = cardiopharyngeal
# 4 = cardiac crest #2

# TGFb
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "TGFb")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "TGFb")
#> Comparing communications on a merged object
gg1 + gg2

# BMP
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "BMP")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "BMP")
#> Comparing communications on a merged object
gg1 + gg2

# GDF
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "GDF")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "GDF")
#> Comparing communications on a merged object
gg1 + gg2

# WNT
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "WNT")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "WNT")
#> Comparing communications on a merged object
gg1 + gg2

# ncWNT
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "ncWNT")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "ncWNT")
#> Comparing communications on a merged object
gg1 + gg2

# EGF
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "EGF")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "EGF")
#> Comparing communications on a merged object
gg1 + gg2

# NRG
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "NRG")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "NRG")
#> Comparing communications on a merged object
gg1 + gg2

# FGF
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "FGF")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "FGF")
#> Comparing communications on a merged object
gg1 + gg2

# PDGF
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "PDGF")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "PDGF")
#> Comparing communications on a merged object
gg1 + gg2

# IGF
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "IGF")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "IGF")
#> Comparing communications on a merged object
gg1 + gg2

# CXCL
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "CXCL")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "CXCL")
#> Comparing communications on a merged object
gg1 + gg2

# MIF
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "MIF")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "MIF")
#> Comparing communications on a merged object
gg1 + gg2


# EDN
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "EDN")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "EDN")
#> Comparing communications on a merged object
gg1 + gg2

# SEMA3
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "SEMA3")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "SEMA3")
#> Comparing communications on a merged object
gg1 + gg2

# EPHA
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "EPHA")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "EPHA")
#> Comparing communications on a merged object
gg1 + gg2

# NOTCH
gg1 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T, signaling = "NOTCH")
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_merged, sources.use = c(5:9), targets.use = c(3),comparison = c(2,1), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T,signaling = "NOTCH")
#> Comparing communications on a merged object
gg1 + gg2






