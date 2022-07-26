#### export peaks as bed file
###### 1. all peaks in WT atlas

# look at the peak set and then save it as a separate object
getPeakSet(wt_atlas)
#GRanges object with 548312 ranges and 13 metadata columns:
merged_peak_setGR <- getPeakSet(wt_atlas)



## Export all non-overlapping 548312 peaks
setwd("/Users/sranade/scATAC-seq/new_wt_section_211111/Exports")
merged_peak_df <- as.data.frame(merged_peak_setGR@ranges)

merged_peak_df <- data.frame(seqnames=seqnames(merged_peak_setGR),
                 starts=start(merged_peak_setGR)-1,
                 ends=end(merged_peak_setGR),
                 names=c(rep(".", length(merged_peak_setGR))),
                 scores=merged_peak_setGR@elementMetadata@listData[["score"]],
                 strands=merged_peak_setGR@strand@values)

# format in same manner as you would use for homer 
merged_peak_df <- merged_peak_df %>% mutate(id = row_number())
merged_peak_df$blankVar <- NA
merged_peak_df_homerformat <- merged_peak_df[,c(1,2,3,7,8,6)]
write.table(merged_peak_df_homerformat, file="allPeaks_211215.bed", quote=F, sep="\t", row.names=F, col.names=F)


###################################################
# 2. all cluster enriched peaks
# Export cluster specific peaks for 11 clusters:
markersPeaks <- getMarkerFeatures(
  ArchRProj = wt_atlas, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
# class: SummarizedExperiment 
# dim: 548312 11 
# metadata(2): MatchInfo Params
# assays(7): Log2FC Mean ... AUC MeanBGD
# rownames(548312): 1 2 ... 548311 548312
# rowData names(4): seqnames idx start end
# colnames(11): Blood Cardiomyocyte ... Paraxial_Mesoderm SHF_Progenitor
# colData names(0):


markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1.5", returnGR = TRUE)
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.5",
  transpose = F,
  nLabel = 1,
  nPrint = 1)
# Identified 205037 markers!

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#################
## export marker peaks and sort by meanDiff and format for homer
library(dplyr)
EndMT_df <- as.data.frame(markerList@listData[["EndMT"]])
Endocardium_df <- as.data.frame(markerList@listData[["Endothelium"]])
LateralPlateMesoderm_df <- as.data.frame(markerList@listData[["LPM"]])
Cardiomyocyte_df <- as.data.frame(markerList@listData[["Cardiomyocyte"]])
Epicardium_df <- as.data.frame(markerList@listData[["Epicardium"]])
Blood_df <- as.data.frame(markerList@listData[["Blood"]])
NeuralCrest_df <- as.data.frame(markerList@listData[["Neural_Crest"]])
SHF_Progenitor_df <- as.data.frame(markerList@listData[["SHF_Progenitor"]])
Ectoderm_df <- as.data.frame(markerList@listData[["Ectoderm"]])
Endoderm_df <- as.data.frame(markerList@listData[["Endoderm"]])
ParaxialMesoderm_df <- as.data.frame(markerList@listData[["Paraxial_Mesoderm"]])

# add unique row id
EndMT_df <- EndMT_df %>% mutate(id = row_number())
Endocardium_df<- Endocardium_df %>% mutate(id = row_number())
LateralPlateMesoderm_df<- LateralPlateMesoderm_df %>% mutate(id = row_number())
Cardiomyocyte_df<- Cardiomyocyte_df %>% mutate(id = row_number())
Epicardium_df<- Epicardium_df %>% mutate(id = row_number())
Blood_df<- Blood_df %>% mutate(id = row_number())
NeuralCrest_df<- NeuralCrest_df %>% mutate(id = row_number())
SHF_Progenitor_df<- SHF_Progenitor_df %>% mutate(id = row_number())
Ectoderm_df<- Ectoderm_df %>% mutate(id = row_number())
Endoderm_df<- Endoderm_df %>% mutate(id = row_number())
ParaxialMesoderm_df<- ParaxialMesoderm_df %>% mutate(id = row_number())



# add blank col
EndMT_df$blankVar <- NA
Endocardium_df$blankVar <- NA
LateralPlateMesoderm_df$blankVar <- NA
Cardiomyocyte_df$blankVar <- NA
Epicardium_df$blankVar <- NA
Blood_df$blankVar <- NA
NeuralCrest_df$blankVar <- NA
SHF_Progenitor_df$blankVar <- NA
Ectoderm_df$blankVar <- NA
Endoderm_df$blankVar <- NA
ParaxialMesoderm_df$blankVar <- NA


# format for homer
EndMT_df <- EndMT_df[,c(1,2,3,9,10,5)]
Endocardium_df<- Endocardium_df[,c(1,2,3,9,10,5)]
LateralPlateMesoderm_df<- LateralPlateMesoderm_df[,c(1,2,3,9,10,5)]
Cardiomyocyte_df<- Cardiomyocyte_df[,c(1,2,3,9,10,5)]
Epicardium_df<- Epicardium_df[,c(1,2,3,9,10,5)]
Blood_df<- Blood_df[,c(1,2,3,9,10,5)]
NeuralCrest_df<- NeuralCrest_df[,c(1,2,3,9,10,5)]
SHF_Progenitor_df<- SHF_Progenitor_df[,c(1,2,3,9,10,5)]
Ectoderm_df<- Ectoderm_df[,c(1,2,3,9,10,5)]
Endoderm_df<- Endoderm_df[,c(1,2,3,9,10,5)]
ParaxialMesoderm_df<- ParaxialMesoderm_df[,c(1,2,3,9,10,5)]



# save
setwd("/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/cluster_enriched_peaks_211215")
write.table(EndMT_df, file="EndMT_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(Endocardium_df, file="Endocardium_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(LateralPlateMesoderm_df, file="LateralPlateMesoderm_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(Cardiomyocyte_df, file="Cardiomyocyte_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(Epicardium_df, file="Epicardium_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(Blood_df, file="Blood_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(NeuralCrest_df, file="NeuralCrest_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(SHF_Progenitor_df, file="SHF_Progenitor_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(Ectoderm_df, file="Ectoderm_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(Endoderm_df, file="Endoderm_df.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(ParaxialMesoderm_df, file="ParaxialMesoderm_df.bed", quote=F, sep="\t", row.names=F, col.names=F)




