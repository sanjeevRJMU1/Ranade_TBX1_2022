Tbx1_CNC <- readRDS("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_CNC_annotated.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()

FeaturePlot(Tbx1_CNC, c("Sema3c"))
FeaturePlot(Tbx1_CNC, c("Nrp1"))
FeaturePlot(Tbx1_CNC, c("Nrp2"))
FeaturePlot(Tbx1_CNC, c("Plxna2"))
FeaturePlot(Tbx1_CNC, c("Plxna4"))
FeaturePlot(Tbx1_CNC, c("Sema3c"))
VlnPlot(Tbx1_CNC, c("Plxna2","Sema3c"), pt.size = 0.1)

features=c("Fgfr2","Fgfr1")
DotPlot(Tbx1_CNC, features = features, cols = c("grey","red")) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev)
