###### DAR/DGE of Tbx1 motifs in AHF 
setwd("/Users/sranade/scATAC-seq/Mesoderm_subset3_210517")

library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(data.table)


######## 
# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517")
mesoderm_subset
######## 

#### after running homer annotate on motif1 (Tbox motif)
AHF_motif1 <- read.csv(file="/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_220110/AHF/homer_annotate/AHF_E95_down_motif1_outputfile.txt", sep = '\t')
# filter rows without a motif
AHF_motif1_filtered <- AHF_motif1[AHF_motif1$X1.GTGHTGACAGAT.BestGuess.Tbx20.T.box..Heart.Tbx20.ChIP.Seq.GSE29636..Homer.0.874..Distance.From.Peak.sequence.strand.conservation.!="",]
write.csv(AHF_motif1_filtered, file = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_220110/AHF/DAR_DGE/AHF_motif1_filtered.csv", row.names = F)

####################################################
### Intersecting DAR with DGE
library(dplyr)
library(data.table)
library(purrr)


# Read in your list of genes that are differentially expressed per cluster
AHF_DGE_list <- read.csv("/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_220110/AHF/homer_annotate/AHF_E95_DGE_list.txt", header = F)
AHF_DGE_list<-dplyr::pull(AHF_DGE_list)


#### after running homer annotate on motif1 (Tbox motif)
AHF_motif1 <- read.csv(file="/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_220110/AHF/homer_annotate/AHF_E95_down_motif1_outputfile.txt", sep = '\t')
# filter rows without a motif
AHF_motif1_filtered <- AHF_motif1[AHF_motif1$X1.GTGHTGACAGAT.BestGuess.Tbx20.T.box..Heart.Tbx20.ChIP.Seq.GSE29636..Homer.0.874..Distance.From.Peak.sequence.strand.conservation.!="",]

# subset where "Genes" is the column in the df that contains the gene names

AHF_Tbx1_DAR_DGE_match <- subset(AHF_motif1_filtered, AHF_motif1_filtered$Gene.Name %in% AHF_DGE_list)
AHF_Tbx1_DAR_DGE_match <- AHF_Tbx1_DAR_DGE_match[,c(2,3,4,6,5,9,10,16)]

# sort by genes with multiple peaks
AHF_Tbx1_DAR_DGE_match_sorted <- as.data.frame(sort(table(AHF_Tbx1_DAR_DGE_match$Gene.Name)))
AHF_Tbx1_DAR_DGE_match_sorted<-AHF_Tbx1_DAR_DGE_match_sorted[!(AHF_Tbx1_DAR_DGE_match_sorted$Freq==0),]
AHF_Tbx1_DAR_DGE_match_sorted <- map_dfr(AHF_Tbx1_DAR_DGE_match_sorted, rev)


## write
write.csv(AHF_Tbx1_DAR_DGE_match,file = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_220110/AHF/DAR_DGE/AHF_DARDGE_unfiltered.csv")

# ## read back list where only intronic and intergenic peaks are kept
# AHF_filtered_list <- read.csv("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/output_results/AHF_DARDGE_filtered_list.csv", header = F)
# length(table(AHF_filtered_list))
#
# AHF_unfiltered_list <- read.csv("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/output_results/AHF_DARDGE_unfiltered_list.csv", header = F)
# length(table(AHF_unfiltered_list))

# trim the columns
# 
# AHF_Tbx1_DAR_DGE_match <-AHF_Tbx1_DAR_DGE_match[AHF_Tbx1_DAR_DGE_match$Detailed.Annotation == "Intergenic",]




# sort by genes with multiple peaks
AHF_Tbx1_DAR_DGE_match_sorted <- as.data.frame(sort(table(AHF_Tbx1_DAR_DGE_match$Gene.Name)))
AHF_Tbx1_DAR_DGE_match_sorted<-AHF_Tbx1_DAR_DGE_match_sorted[!(AHF_Tbx1_DAR_DGE_match_sorted$Freq==0),]
AHF_Tbx1_DAR_DGE_match_sorted <- map_dfr(AHF_Tbx1_DAR_DGE_match_sorted, rev)



AHF_Tbx1_DAR_DGE_match_sorted_trim<-AHF_Tbx1_DAR_DGE_match_sorted[,1]
write.csv(AHF_Tbx1_DAR_DGE_match_sorted_trim, file = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_220110/AHF/DAR_DGE/AHF_TF_sorted.csv", quote = F, row.names = F, sep = '\t')


##################


AHF_subset <- SubsetData(Tbx1_mesoderm_cnc, ident.use = c("AHF"))
AHF_subset_WT <- OldWhichCells(AHF_subset, subset.name = "gem.group", accept.value = c("WT_E925"))
AHF_subset_KO <- OldWhichCells(AHF_subset, subset.name = "gem.group", accept.value = c("KO_E925"))
Idents(object = AHF_subset, cells= AHF_subset_WT) <- "WT_E925"
Idents(object = AHF_subset, cells= AHF_subset_KO) <- "KO_E925"

AHF_subset <- SubsetData(AHF_subset, ident.use = c("KO_E925","WT_E925"))
AHF_subset@active.ident<-droplevels(AHF_subset@active.ident)
table(AHF_subset@active.ident)


my_levels<-c("WT_E925","KO_E925")
AHF_subset@active.ident <- factor(x = AHF_subset@active.ident, levels = my_levels)
c7_dot <- c("Sema3c","Sema6d","Sema6a","Sema3e","Fgf8","Fgf10","Fgf1","Fgf13","Wnt5a","Tgfb1","Aldh1a2","Crabp1","Efnb2","Plxna4","Plxna2","Slit1","Lef1","Hoxb1","Nkx2-6","Isl1","Gata6","Ebf1","Ebf2","Pitx2","Nr2f2","Six1")
DotPlot(AHF_subset,cols = c("grey","red") ,features = c7_dot) + coord_flip()

no_tf_dot <- c("Sema3c","Sema6d","Sema6a","Sema3e","Efnb2","Slit1","Slit2","Robo2","Plxna4","Plxna2","Fgf8","Fgf10","Fgf1","Fgf13","Wnt5a","Tgfb1","Aldh1a2","Crabp1")
DotPlot(AHF_subset,cols = c("grey","red") ,features = no_tf_dot) + coord_flip()




