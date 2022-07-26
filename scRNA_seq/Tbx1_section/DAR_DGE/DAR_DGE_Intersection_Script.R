### Intersecting DAR with DGE
library(dplyr)
library(data.table)
library(purrr)


# Read in your list of genes that are differentially expressed per cluster
AHF_DGE_list <- read.csv("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/AHF_DGE_FC0-3.csv", header = F)
AHF_DGE_list<-dplyr::pull(AHF_DGE_list)


# Read in your Tbx1 motif containing peaks from homer annotate
AHF_Tbx1_DAR <- read.csv(file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/motif_peaks/AHF/AHF_motif1_filtered.csv")

# subset where "Genes" is the column in the df that contains the gene names

AHF_Tbx1_DAR_DGE_match <- subset(AHF_Tbx1_DAR, AHF_Tbx1_DAR$Gene.Name %in% AHF_DGE_list)
AHF_Tbx1_DAR_DGE_match <- AHF_Tbx1_DAR_DGE_match[,c(2,3,4,6,5,9,10,16)]

## write
write.csv(AHF_Tbx1_DAR_DGE_match,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/output_results/AHF_DARDGE_unfiltered.csv")

## read back list where only intronic and intergenic peaks are kept
AHF_filtered_list <- read.csv("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/output_results/AHF_DARDGE_filtered_list.csv", header = F)
length(table(AHF_filtered_list))

AHF_unfiltered_list <- read.csv("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/output_results/AHF_DARDGE_unfiltered_list.csv", header = F)
length(table(AHF_unfiltered_list))

# trim the columns

AHF_Tbx1_DAR_DGE_match <-AHF_Tbx1_DAR_DGE_match[AHF_Tbx1_DAR_DGE_match$Detailed.Annotation == "Intergenic",] 




# sort by genes with multiple peaks
AHF_Tbx1_DAR_DGE_match_sorted <- as.data.frame(sort(table(AHF_Tbx1_DAR_DGE_match$Gene.Name)))
AHF_Tbx1_DAR_DGE_match_sorted<-AHF_Tbx1_DAR_DGE_match_sorted[!(AHF_Tbx1_DAR_DGE_match_sorted$Freq==0),]
AHF_Tbx1_DAR_DGE_match_sorted <- map_dfr(AHF_Tbx1_DAR_DGE_match_sorted, rev)




write.table(AHF_TF_sorted, file = "AHF_TF_sorted.csv", quote = F, row.names = F, sep = '\t')



####################################################
### Intersecting DAR with DGE
library(dplyr)
library(data.table)
library(purrr)


# Read in your list of genes that are differentially expressed per cluster
PharyngealMesoderm_DGE_list <- read.csv("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/PharyngealMesoderm_DGE.csv", header = F)
PharyngealMesoderm_DGE_list<-dplyr::pull(PharyngealMesoderm_DGE_list)


# Read in your Tbx1 motif containing peaks from homer annotate
PharyngealMesoderm_Tbx1_DAR <- read.csv(file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/motif_peaks/PharyngealMesoderm/PharyngealMesoderm_motif1_filtered.csv")

# subset where "Genes" is the column in the df that contains the gene names

PharyngealMesoderm_Tbx1_DAR_DGE_match <- subset(PharyngealMesoderm_Tbx1_DAR, PharyngealMesoderm_Tbx1_DAR$Gene.Name %in% PharyngealMesoderm_DGE_list)
PharyngealMesoderm_Tbx1_DAR_DGE_match <- PharyngealMesoderm_Tbx1_DAR_DGE_match[,c(2,3,4,6,5,9,10,16)]

## write
write.csv(PharyngealMesoderm_Tbx1_DAR_DGE_match,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/output_results/PharyngealMesoderm_DARDGE_unfiltered.csv")

## read back list where only intronic and intergenic peaks are kept
PharyngealMesoderm_filtered_list <- read.csv("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/output_results/PharyngealMesoderm_DARDGE_filtered_list.csv", header = F)
length(table(PharyngealMesoderm_filtered_list))

PharyngealMesoderm_unfiltered_list <- read.csv("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_DGE_Intersection/output_results/PharyngealMesoderm_DARDGE_unfiltered_list.csv", header = F)
PharyngealMesoderm_unfiltered_list<-dplyr::pull(PharyngealMesoderm_unfiltered_list)
table(PharyngealMesoderm_unfiltered_list)
length(table(PharyngealMesoderm_unfiltered_list))

# trim the columns

PharyngealMesoderm_Tbx1_DAR_DGE_match <-PharyngealMesoderm_Tbx1_DAR_DGE_match[PharyngealMesoderm_Tbx1_DAR_DGE_match$Detailed.Annotation == "Intergenic",] 




# sort by genes with multiple peaks
PharyngealMesoderm_Tbx1_DAR_DGE_match_sorted <- as.data.frame(sort(table(PharyngealMesoderm_Tbx1_DAR_DGE_match$Gene.Name)))
PharyngealMesoderm_Tbx1_DAR_DGE_match_sorted<-PharyngealMesoderm_Tbx1_DAR_DGE_match_sorted[!(PharyngealMesoderm_Tbx1_DAR_DGE_match_sorted$Freq==0),]
PharyngealMesoderm_Tbx1_DAR_DGE_match_sorted <- map_dfr(PharyngealMesoderm_Tbx1_DAR_DGE_match_sorted, rev)




write.table(PharyngealMesoderm_TF_sorted, file = "PharyngealMesoderm_TF_sorted.csv", quote = F, row.names = F, sep = '\t')









