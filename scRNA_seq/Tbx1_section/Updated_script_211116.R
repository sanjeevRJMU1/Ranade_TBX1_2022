# on 11/16 I decided to revist the annotations. Mostly focused on what discriminates the "mesoderm/mesenchyme" clusters
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)

# load Mesoderm-cnc unannotated object
Tbx1_mesoderm_cnc_unannotated_210829 <- readRDS("~/scRNA-seq/2020_Tbx1/all_aggr_210429/seurat_objects/Tbx1_mesoderm_cnc_unannotated_210829.RDS")
DimPlot(Tbx1_mesoderm_cnc_unannotated_210829, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()
Tbx1_mesoderm_cnc_unannotated_210829.markers <- FindAllMarkers(Tbx1_mesoderm_cnc_unannotated_210829, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tbx1_mesoderm_cnc_unannotated_210829.markers, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/findallmarkers_outputs/Tbx1_mesoderm_cnc_unannotated_210829.markers_211116.csv")

table(Tbx1_mesoderm_cnc_unannotated_210829@active.ident,Tbx1_mesoderm_cnc_unannotated_210829$gem.group)
#     KO_E105 KO_E115 KO_E925 WT_E105 WT_E115 WT_E925
# 0      735    2293     299     283     870     175
# 1      728    1549     284     563     459     300
# 2      536    1561     169     409     445     149
# 3      385    1627     183     291     608     146
# 4      692     186     187     661    1197      97
# 5      408    1298      83     300     872      56
# 6      507     337       0     600     910       1
# 7      279     799     101     192     701      87
# 8      525     156      87     478     868      32
# 9      433     603     294     310     243     183
# 10     140     803      53     106     404      39
# 11     126     346     146     240     543     109
# 12     164     562     294      64      90     270
# 13     298     405     208     198     177     132
# 14     149     570       0      72     162       0
# 15     278     331      68     171      23      42
# 16     148     213      44     110      20      40


## marker test
clust5v10 <- FindMarkers(Tbx1_mesoderm_cnc_unannotated_210829, ident.1 = 5, ident.2 = 10)
clust1v2 <- FindMarkers(Tbx1_mesoderm_cnc_unannotated_210829, ident.1 = 1, ident.2 = 2)


