#####
## re-start 08/29
Tbx1_mesoderm_cnc <- readRDS(file = "Tbx1_mesoderm_cnc_annotated.RDS")
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
Tbx1_mesoderm_cnc.markers <- FindAllMarkers(Tbx1_mesoderm_cnc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VlnPlot(Tbx1_mesoderm_cnc, c("nCount_RNA"), pt.size = 0.1, ncol = 1)

Tbx1_mesoderm_cnc <- RunUMAP(Tbx1_mesoderm_cnc, dims = 1:20, reduction = "mnn", verbose = TRUE)
Tbx1_mesoderm_cnc <- FindNeighbors(Tbx1_mesoderm_cnc,reduction = "mnn", dims = 1:20, verbose = TRUE, force.recalc = T)
Tbx1_mesoderm_cnc <- FindClusters(Tbx1_mesoderm_cnc, verbose = TRUE, resolution = 0.5)
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = T)
saveRDS(Tbx1_mesoderm_cnc, file = "Tbx1_mesoderm_cnc_unannotated_210829.RDS")

new.cluster.ids <- c("Posterior_Second_Heart_Field",
                     "Paraxial_Mesoderm",
                     "Cardiopharyngeal_Mesoderm",
                     "Cardiopharyngeal_Mesenchyme",
                     "Neural_Crest",
                     "Cranial_Mesenchyme",
                     "Neural_Crest",
                     "Smooth_Muscle",
                     "Neural_Crest",
                     "Epicardium",
                     "Cardiopharyngeal_Mesenchyme",
                     "Anterior_Second_Heart_Field",
                     "Cardiomyocyte",
                     "Neural_Crest",
                     "Posterior_Second_Heart_Field",
                     "Epicardium",
                     "Epicardium")
names(new.cluster.ids) <- levels(Tbx1_mesoderm_cnc)
Tbx1_mesoderm_cnc <- RenameIdents(Tbx1_mesoderm_cnc, new.cluster.ids)
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()

features <- c(
  "Foxf1","Osr1",
  "Foxd1","Meox2","Zic1",
  "Tbx1","Foxc1","Foxc2",
  "Lix1","Msx2","Alx1",
  "Dlx2","Dlx5",
  "Bdnf","Ptx3",
  "Rgs5","Isl1","Hand2",
  "Tbx18","Wt1",
  "Myog","Myf5","Eya4","Pax3",
  "Actc1","Tnnt2"
)
DotPlot(Tbx1_mesoderm_cnc, features = features, cols = c("grey","red")) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev)

#
Tbx1_mesoderm_cnc.markers <- FindAllMarkers(Tbx1_mesoderm_cnc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_mesoderm_cnc.markers, file = "Tbx1_mesoderm_cnc_annotated-markers_210830.csv")
saveRDS(Tbx1_mesoderm_cnc, file = "Tbx1_mesoderm_cnc_annotated_210830.RDS")



DoHeatmap(subset(Tbx1_mesoderm_cnc, downsample=100), features = features,slot = "scale.data",assay = "SCT", label = F) + scale_fill_gradientn(colors = c("cornflowerblue", "white", "red")) 

