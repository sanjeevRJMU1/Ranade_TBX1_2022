####### plot UMAPs split by time point, label = F, nolegend

# mesoderm subset
table(Tbx1_mesoderm_cnc$gem.group)
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap",split.by = "gem.group", label = F, pt.size = 0.2, ncol = 3) + NoLegend()

# CNC subset
table(Tbx1_CNC$gem.group)
DimPlot(Tbx1_CNC, reduction = "umap",split.by = "gem.group", label = F, pt.size = 0.2, ncol = 3) + NoLegend()
