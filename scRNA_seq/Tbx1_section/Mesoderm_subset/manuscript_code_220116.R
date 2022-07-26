
#####################
# re-starting on 220103

setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

## re-start 08/29
Tbx1_mesoderm_cnc <- readRDS(file = "Tbx1_mesoderm_cnc_annotated_210830.RDS")
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
table(Tbx1_mesoderm_cnc@active.ident,Tbx1_mesoderm_cnc$gem.group)





cran_vs_cpm_mes <- FindMarkers(Tbx1_mesoderm_cnc, ident.1 ="Cranial_Mesenchyme",ident.2 = "Cardiopharyngeal_Mesenchyme")




card_vs_parax <- FindMarkers(Tbx1_mesoderm_cnc, ident.1 ="Cardiopharyngeal_Mesoderm",ident.2 = "Paraxial_Mesoderm")

#######
Tbx1_mesoderm_cnc <- readRDS(file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_mesoderm_cnc_annotated_210830.RDS")
DimPlot(Tbx1_mesoderm_cnc, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
table(Tbx1_mesoderm_cnc@active.ident, Tbx1_mesoderm_cnc$gem.group)

my_order <- c("Cardiomyocyte",
              "Posterior_Second_Heart_Field",
              "Anterior_Second_Heart_Field",
              "Cardiopharyngeal_Mesoderm",
              "Paraxial_Mesoderm",
              "Cardiopharyngeal_Mesenchyme",
              "Cranial_Mesenchyme",
              "Epicardium",
              "Smooth_Muscle",
              "Neural_Crest"
)
Tbx1_mesoderm_cnc@active.ident <- factor(x = Tbx1_mesoderm_cnc@active.ident, levels = my_order)


features <- c(
  "Actc1","Tnnt2",
  "Foxf1","Osr1",
  "Mybpc1","Tbx1","Fgf8","Fgf10","Six2","Six1","Eya4",
  "Foxd1","Meox2","Zic1","Foxc1","Foxc2","Irx2","Hoxd4",
  "Lix1","Msx2","Alx1","Alx4",
  "Bdnf","Ptx3",
  "Tbx18","Wt1",
  "Rgs5","Isl1","Hand2",
  "Dlx2","Dlx5"
)



DotPlot(Tbx1_mesoderm_cnc, features = features, cols = c("grey","red")) + RotatedAxis() + coord_flip()

DotPlot(Tbx1_mesoderm_cnc, features = features, cols = c("grey","red")) + RotatedAxis() 
DoHeatmap(subset(Tbx1_mesoderm_cnc, downsample=100), features = features,slot = "scale.data",assay = "SCT", label = T) + scale_fill_gradientn(colors = c("cornflowerblue", "white", "red")) + NoLegend() 


