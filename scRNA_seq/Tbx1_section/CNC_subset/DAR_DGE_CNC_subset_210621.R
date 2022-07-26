###### DAR/DGE
setwd("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)


Tbx1_CNC <-readRDS(file = "Tbx1_CNC_annotated.RDS")
DimPlot(Tbx1_CNC, reduction = "umap", label = TRUE, pt.size = 0.2) + NoLegend()

#### after running homer annotate on motif1 (Tbox motif)
CNC_motif1 <- read.csv(file="/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/homer_annotate/Clusters_WTvKO/PA3_motifs/CNC_motif1_outputfile.txt", sep = '\t')
CNC_motif2 <- read.csv(file="/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/homer_annotate/Clusters_WTvKO/PA3_motifs/CNC_motif2_outputfile.txt", sep = '\t')
CNC_motif3 <- read.csv(file="/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/homer_annotate/Clusters_WTvKO/PA3_motifs/CNC_motif3_outputfile.txt", sep = '\t')
CNC_motif3 <- read.csv(file="/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/homer_annotate/Clusters_WTvKO/PA3_motifs/CNC_motif3_outputfile.txt", sep = '\t')



# filter rows without a motif
CNC_motif1_filtered <- CNC_motif1[CNC_motif1$X1.NGTAATTA.BestGuess.DLX5.Homeobox..BasalGanglia.Dlx5.ChIP.seq.GSE124936..Homer.0.987..Distance.From.Peak.sequence.strand.conservation.!="",]
write.csv(CNC_motif1_filtered, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/homer_annotate/Clusters_WTvKO/PA3_motifs/CNC_motif1_filtered.csv", row.names = F)
CNC_motif2_filtered <- CNC_motif2[CNC_motif2$X3.CTGKTTTTTBTC.BestGuess.Hand2.bHLH..Mesoderm.Hand2.ChIP.Seq.GSE61475..Homer.0.929..Distance.From.Peak.sequence.strand.conservation.!="",]
write.csv(CNC_motif2_filtered, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/homer_annotate/Clusters_WTvKO/PA3_motifs/CNC_motif2_filtered.csv", row.names = F)
CNC_motif3_filtered <- CNC_motif3[CNC_motif3$X5.TTGTTTAC.BestGuess.FOXO3.MA0157.2.Jaspar.0.978..Distance.From.Peak.sequence.strand.conservation.!="",]
write.csv(CNC_motif3_filtered, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/homer_annotate/Clusters_WTvKO/PA3_motifs/CNC_motif3_filtered.csv", row.names = F)



####################################################
### Intersecting DAR with DGE
library(dplyr)
library(data.table)
library(purrr)


# Read in your list of genes that are differentially expressed per cluster
CNC_DGE_list <- read.csv("/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/annotated/E115/PA3_DGE_list_down.txt", header = F)
CNC_DGE_list<-dplyr::pull(CNC_DGE_list)


# Read in your Tbx1 motif containing peaks from homer annotate
#CNC_motif1 <- read.csv(file="/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports/CNC/homer_annotate/CNC_motif1_outputfile.txt", sep = '\t')
# filter rows without a motif
#CNC_motif1_filtered <- CNC_motif1[CNC_motif1$X1.GTGWTGAMAGAT.BestGuess.Tbx20.T.box..Heart.Tbx20.ChIP.Seq.GSE29636..Homer.0.856..Distance.From.Peak.sequence.strand.conservation.!="",]

# subset where "Genes" is the column in the df that contains the gene names

CNC_Dlx_DAR_DGE_match <- subset(CNC_motif1_filtered, CNC_motif1_filtered$Gene.Name %in% CNC_DGE_list)
CNC_Dlx_DAR_DGE_match <- CNC_Dlx_DAR_DGE_match[,c(2,3,4,6,5,9,10,16)]

CNC_Hand_DAR_DGE_match <- subset(CNC_motif2_filtered, CNC_motif2_filtered$Gene.Name %in% CNC_DGE_list)
CNC_Hand_DAR_DGE_match <- CNC_Hand_DAR_DGE_match[,c(2,3,4,6,5,9,10,16)]

CNC_Fox_DAR_DGE_match <- subset(CNC_motif3_filtered, CNC_motif3_filtered$Gene.Name %in% CNC_DGE_list)
CNC_Fox_DAR_DGE_match <- CNC_Fox_DAR_DGE_match[,c(2,3,4,6,5,9,10,16)]


# sort by genes with multiple peaks
CNC_Dlx_DAR_DGE_match_sorted <- as.data.frame(sort(table(CNC_Dlx_DAR_DGE_match$Gene.Name)))
CNC_Dlx_DAR_DGE_match_sorted<-CNC_Dlx_DAR_DGE_match_sorted[!(CNC_Dlx_DAR_DGE_match_sorted$Freq==0),]
CNC_Dlx_DAR_DGE_match_sorted <- map_dfr(CNC_Dlx_DAR_DGE_match_sorted, rev)


write.csv(CNC_Dlx_DAR_DGE_match_sorted, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/annotated/E115/CNC_Dlx_sorted.csv", quote = F, row.names = F, sep = '\t')

# subset where "Genes" is the column in the df that contains the gene names

CNC_Hand_DAR_DGE_match <- subset(CNC_motif2_filtered, CNC_motif2_filtered$Gene.Name %in% CNC_DGE_list)
CNC_Hand_DAR_DGE_match <- CNC_Hand_DAR_DGE_match[,c(2,3,4,6,5,9,10,16)]


# sort by genes with multiple peaks
CNC_Hand_DAR_DGE_match_sorted <- as.data.frame(sort(table(CNC_Hand_DAR_DGE_match$Gene.Name)))
CNC_Hand_DAR_DGE_match_sorted<-CNC_Hand_DAR_DGE_match_sorted[!(CNC_Hand_DAR_DGE_match_sorted$Freq==0),]
CNC_Hand_DAR_DGE_match_sorted <- map_dfr(CNC_Hand_DAR_DGE_match_sorted, rev)


write.csv(CNC_Hand_DAR_DGE_match_sorted, file = "/Users/sranade/scRNA-seq/2020_Tbx1/all_aggr_210429/CNC_DGE/annotated/E115/CNC_Hand_sorted.csv", quote = F, row.names = F)


##################
CNC_subset <- SubsetData(Tbx1_mesoderm_cnc, ident.use = c("CNC"))
CNC_subset_WT <- OldWhichCells(CNC_subset, subset.name = "gem.group", accept.value = c("WT_E925"))
CNC_subset_KO <- OldWhichCells(CNC_subset, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = CNC_subset, cells= CNC_subset_WT) <- "WT_E925"
Idents(object = CNC_subset, cells= CNC_subset_KO) <- "KO_E925"

CNC_subset <- SubsetData(CNC_subset, ident.use = c("KO_E925","WT_E925"))
table(CNC_subset@active.ident)


my_levels<-c("WT_E925","KO_E925","CNC")
CNC_subset@active.ident <- factor(x = CNC_subset@active.ident, levels = my_levels)
c7_dot <- c("Sema3c","Sema6d","Sema6a","Sema3e","Fgf8","Fgf10","Fgf1","Fgf13","Wnt5a","Tgfb1","Aldh1a2","Crabp1","Efnb2","Plxna4","Plxna2","Slit1","Lef1","Hoxb1","Nkx2-6","Isl1","Gata6","Ebf1","Ebf2","Pitx2","Nr2f2","Six1")
DotPlot(CNC_subset,cols = c("grey","red") ,features = c7_dot) + coord_flip()

no_tf_dot <- c("Sema3c","Sema6d","Sema6a","Sema3e","Fgf8","Fgf10","Fgf1","Fgf13","Wnt5a","Tgfb1","Aldh1a2","Crabp1","Efnb2","Plxna4","Plxna2","Slit1")
DotPlot(CNC_subset,cols = c("grey","red") ,features = no_tf_dot) + coord_flip()




