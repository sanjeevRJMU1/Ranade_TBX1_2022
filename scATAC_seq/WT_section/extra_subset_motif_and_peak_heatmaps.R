### subsetting by full rowname match
#### look thru enrichMotifs and find the ones that are the most enriched per cluster. choose 3 per group and then put in that
#### full name for subseting

# used homer for motif enrichment
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


wt_atlas <- addMotifAnnotations(ArchRProj = wt_atlas, motifSet = "cisbp", name = "Motif", force = T)
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = wt_atlas,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)
enrichMotifs
# class: SummarizedExperiment 
# dim: 332 11 
# metadata(0):
#   assays(10): mlog10Padj mlog10p ... CompareFrequency feature
# rownames(332): AP.1.bZIP_1 AP.2gamma.AP2_2 ... ZNF711.Zf_331 ZSCAN22.Zf_332
# rowData names(0):
#   colnames(11): Blood Cardiomyocyte ... Paraxial_Mesoderm SHF_Progenitor
# colData names(0):

#motif_df <- as.data.frame(enrichMotifs@assays$data$mlog10Padj)
motif_df <- as.data.frame(enrichMotifs@assays@data@listData[["mlog10Padj"]])

write.table(motif_df, file="/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/enrichmotif_logP_cisbp.txt", quote=F, sep="\t", row.names=T, col.names=T)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

############ Homer Version
### TODO: add support for partial rownames
object2subset <- enrichMotifs
## create vector with FULL name matches to rows of interest
rownames2match <- c("GATA.SCL.Zf.bHLH_111",
                    "EKLF.Zf_65",
                    "Hnf1.Homeobox_125",
                    "Bach1.bZIP_16",
                    "Mef2d.MADS_164",
                    "Tbx20.T.box_284",
                    "Gata4.Zf_108",
                    "Nkx2.1.Homeobox_182",
                    "Sox2.HMG_261",
                    "OCT4.SOX2.TCF.NANOG.POU.Homeobox.HMG_198",
                    "Brn1.POU.Homeobox_27",
                    "HOXA2.Homeobox_128",
                    "c.Jun.CRE.bZIP_142",
                    "Fosl2.bZIP_86",
                    "Atf2.bZIP_11",
                    "Etv2.ETS_82",
                    "GRHL2.CP2_119",
                    "AP.2gamma.AP2_2",
                    "FOXA1.Forkhead_87",
                    "Six2.Homeobox_255",
                    "ERG.ETS_72",
                    "EWS.FLI1.fusion.ETS_84",
                    "ETV1.ETS_81",
                    "Fli1.ETS_85",
                    "Gata2.Zf_102",
                    "TEAD.TEA_294",
                    "NF1.CTF_176",
                    "Tlx..NR_297",
                    "TEAD2.TEA_292",
                    "Gata6.Zf_109",
                    "Tcf3.HMG_288",
                    "Foxo3.Forkhead_96",
                    "Pitx1.Ebox.Homeobox.bHLH_224",
                    "ZBTB18.Zf_307",
                    "NeuroD1.bHLH_173",
                    "Olig2.bHLH_202",
                    "AP.2alpha.AP2_3",
                    "Sox9.HMG_265",
                    "Lhx2.Homeobox_151",
                    "RXR.NR..DR1_251"
                    )
## pull indices where rownames match - full string matches only
indices4subset <- which(rownames(object2subset) %in% rownames2match, arr.ind=TRUE)
## subset SummarizedExperiment by row index
SE_subset <- object2subset[indices4subset, 1:ncol(object2subset)]

### sanity checks
object2subset
SE_subset

heatmapEM <- plotEnrichHeatmap(SE_subset, transpose = F)
order_rows <- c("Cardiomyocyte","Epicardium","LPM","EndMT","Endothelium","Blood","Ectoderm","Endoderm","Neural_Crest","Paraxial_Mesoderm","SHF_Progenitor")
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = wt_atlas, addDOC = FALSE)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot", row_order= order_rows)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = wt_atlas, addDOC = FALSE)



## hack: export all peaks from peak heatmap
peak_matrix_heatmap<-heatmapPeaks@ht_list[["Row Z-Scores
209065 features
PeakMatrix"]]@matrix
peak_matrix_heatmap<-as.data.frame(peak_matrix_heatmap)

peak_matrix_blood <- peak_matrix_heatmap[peak_matrix_heatmap$Blood==1,]
peak_matrix_CM <- peak_matrix_heatmap[peak_matrix_heatmap$Cardiomyocyte==1,]
peak_matrix_Ectoderm <- peak_matrix_heatmap[peak_matrix_heatmap$Ectoderm==1,]
peak_matrix_EndMT <- peak_matrix_heatmap[peak_matrix_heatmap$EndMT==1,]
peak_matrix_Endoderm <- peak_matrix_heatmap[peak_matrix_heatmap$Endoderm==1,]
peak_matrix_Endothelium <- peak_matrix_heatmap[peak_matrix_heatmap$Endothelium==1,]
peak_matrix_Epicardium <- peak_matrix_heatmap[peak_matrix_heatmap$Epicardium==1,]
peak_matrix_ParaxialMesoderm <- peak_matrix_heatmap[peak_matrix_heatmap$ParaxialMesoderm==1,]
peak_matrix_LPM <- peak_matrix_heatmap[peak_matrix_heatmap$LPM==1,]
peak_matrix_SHF_Progenitors <- peak_matrix_heatmap[peak_matrix_heatmap$SHF_Progenitors==1,]
peak_matrix_CNC <- peak_matrix_heatmap[peak_matrix_heatmap$NeuralCrest==1,]

## Format for homer, repeating for each cluster (I really need to learn how to write loop scripts!!)
library(stringr)
result <- data.frame("")
data <- peak_matrix_CNC

data$col2split <- rownames(data)
result <- data.frame(result, # function 1
                     do.call(rbind, # function 2
                             str_split(gsub(data$col2split, # function 3, referred to column by index in original, but it broke. Ha!
                                            pattern = ":",
                                            replacement = "-"),
                                       pattern = "-")
                     )
)
# remove empty column
result[,1] <- NULL
colnames(result) <- c("chr", "start", "stop")
# add unique row id
library(dplyr)
setwd("/Users/sranade/scATAC-seq/wt_atlas_210209/ExportsforFigures/updated_210820/markerpeak_heatmap_annotations/cluster_bed_files")
result$ID <- seq.int(nrow(result))
result$blankVar <- NA
result$Strand <- "*"
write.table(result, file="CNC_Progenitors_peakmatrix_homer.bed", quote=F, sep="\t", row.names=F, col.names=F)
rm(data)
rm(result)



