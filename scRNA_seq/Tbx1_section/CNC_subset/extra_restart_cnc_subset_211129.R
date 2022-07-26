## CNC subset re-start


setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)


##
###
Tbx1_CNC <-readRDS(file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_CNC_annotated.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()

DimPlot(Tbx1_CNC, reduction = "umap",group.by = "gem.group", label = TRUE, pt.size = 0.2) + NoLegend()



# re-run markers and re-make dotplot for figure with more genes 
Tbx1_CNC.markers <- FindAllMarkers(Tbx1_CNC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


#### dot plot for figure
table(Tbx1_CNC@active.ident)
my_order <-c("Migratory_NeuralCrest","Craniofacial_NeuralCrest","Cardiac_NeuralCrest","PA3_Cardiac_NeuralCrest")
Tbx1_CNC@active.ident <- factor(x = Tbx1_CNC@active.ident, levels = my_order)

dotfigures <- c("Sox10","Zeb2","Foxd3","Cdh19","Pax3","Ednrb","Sox8","Sox2",
               "Pantr1","Gm29260","Pou3f3","Smoc1","Emx2","Ebf1","Dlx3","Barx1",
               "Hand1","Hand2","Dlk1","Isl1","Tbx2","Tbx3","Rgs5","Foxf1","Gsc","Gata3","Gata6",
               "Hoxa3","Hoxb3","Hoxd3","Hoxd4","Six2","Shox2","Bdnf")
#pdf(file = "/Users/sranade/Dropbox (Gladstone)/NewManuscriptFigures_211128/Figure_3/clusterdotplot_noleg.pdf", height = 6, width = 4)
DotPlot(Tbx1_CNC, features = dotfigures, cols = c("grey","red"), col.min = -1, col.max = 1.5, dot.scale = 5) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev) + NoLegend()
#dev.off( )

pdf(file = "/Users/sranade/Dropbox (Gladstone)/NewManuscriptFigures_211128/Figure_3/clusterdotplot.pdf", height = 6, width = 4)
DotPlot(Tbx1_CNC, features = dotfigures, cols = c("grey","red"), col.min = -1, col.max = 1.5, dot.scale = 5) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev)
dev.off( )

DotPlot(Tbx1_CNC, features = dotfigures, cols = c("grey","red"), col.min = -1, col.max = 1.5, dot.scale = 5) + RotatedAxis() + scale_x_discrete(limits = rev) + NoLegend()


dotfigures <- c("Twist1","Snai1","Sox10","Foxd3",
                "Dlx2","Dlx5","Sox9",
                "Hand2","Dlk1","Isl1",
                "Hoxa3","Hoxb3","Hoxd4")
DotPlot(Tbx1_CNC, features = dotfigures, cols = c("grey","red"), col.min = -1, col.max = 1.5, dot.scale = 5) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev)


pdf(file = "/Users/sranade/Desktop/Weinstein_Seminar_2022/CNC_cluster_dot_trim.pdf", height = 6, width = 4)
DotPlot(Tbx1_CNC, features = dotfigures, cols = c("grey","red"), col.min = -1, col.max = 1.5, dot.scale = 5) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev) + NoLegend()
dev.off( )


dot_cardiac <- c("Hand1","Hand2","Tbx2","Tbx3","Gata6","Dlk1",
                "Hoxa3","Hoxb3","Hoxd3","Hoxd4")
DotPlot(Tbx1_CNC, features = dot_cardiac, cols = c("grey","red"), col.min = -1, col.max = 1.5, dot.scale = 5) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev) + NoLegend()
DotPlot(Tbx1_CNC, features = dot_cardiac, cols = c("grey","red"), col.min = -1, col.max = 1.5, dot.scale = 5) + RotatedAxis() + coord_flip() + scale_x_discrete(limits = rev)



## check expression of genes from Ghandi/Bronner papers
VlnPlot(Tbx1_CNC, c("Hmga1"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Ets1"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Tgfbr2"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Foxd3"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Sox8"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Pax7"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Pax3"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Tfap2a"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Hoxa1"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Mafb"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Arid3b"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Pantr1","Pou3f3"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Zic1"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Zic2"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Zic3"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Wnt1"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Nkx6-3"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Pdgfra"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Gata3"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Six1"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Six2"), pt.size = 0.1)
VlnPlot(Tbx1_CNC, c("Hoxa2"), pt.size = 0.1)




#
FeaturePlot(Tbx1_CNC, c("Hoxa1"))
FeaturePlot(Tbx1_CNC, c("Hoxb5"))
FeaturePlot(Tbx1_CNC, c("Foxd3"))
FeaturePlot(Tbx1_CNC, c("Pax7"))
FeaturePlot(Tbx1_CNC, c("Tgif1"))
FeaturePlot(Tbx1_CNC, c("Lhx5"))
FeaturePlot(Tbx1_CNC, c("Hoxa2"))
FeaturePlot(Tbx1_CNC, c("Nr6a1"))
FeaturePlot(Tbx1_CNC, c("Phox2a"))
FeaturePlot(Tbx1_CNC, c("Phox2b"))
FeaturePlot(Tbx1_CNC, c("Mafb"))

########### UMAP by time point (all samples collapsed)
table(Tbx1_CNC$gem.group)
condition.vec.01 <- substr(Tbx1_CNC$gem.group,4,7)
table(condition.vec.01)
# condition.vec.01
# E105 E115 E925 
# 4317 5121  793 

Tbx1_CNC$TimePoint <- as.character(condition.vec.01)
table(Tbx1_CNC$TimePoint)
# E10   E11   E77   E85   E92 
# 14228 10021 10483  8537 18261 

DimPlot(Tbx1_CNC, reduction = "umap",group.by = "TimePoint", label = F) + NoLegend()
DimPlot(Tbx1_CNC, reduction = "umap",group.by = "TimePoint", label = F)


########### UMAP by time point x genotype (all samples collapsed)
table(Tbx1_CNC$gem.group)
condition.vec.02 <- substr(Tbx1_CNC$gem.group,1,7)
table(condition.vec.02)
# condition.vec.02
# KO_E105 KO_E115 KO_E925 WT_E105 WT_E115 WT_E925 
# 2208    1418     507    2109    3703     286 


length(unique(condition.vec.02))
#6

Tbx1_CNC$GenotypeTimpoint <- as.character(condition.vec.02)
table(Tbx1_CNC$GenotypeTimpoint)
# E10   E11   E77   E85   E92 
# 14228 10021 10483  8537 18261 

DimPlot(Tbx1_CNC, reduction = "umap",group.by = "GenotypeTimpoint", label = F) + NoLegend()
DimPlot(Tbx1_CNC, reduction = "umap",group.by = "GenotypeTimpoint", label = F)





