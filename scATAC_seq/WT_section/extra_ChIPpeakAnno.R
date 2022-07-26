#### Peak annotation of all peaks
library(ChIPpeakAnno)
## loading packages
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(rtracklayer)
library(GenomicRanges)
library(regioneR)

## export peaks as a bed file
merged_peak_setGR <- getPeakSet(wt_atlas)
merged_peak_df <- as.data.frame(merged_peak_setGR@ranges)
merged_peak_df <- data.frame(seqnames=seqnames(merged_peak_setGR),
                             starts=start(merged_peak_setGR)-1,
                             ends=end(merged_peak_setGR),
                             names=c(rep(".", length(merged_peak_setGR))),
                             scores=merged_peak_setGR@elementMetadata@listData[["score"]],
                             strands=merged_peak_setGR@strand@values)
write.table(merged_peak_df, file="/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/allPeaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

## you can re-load the saved peaks or simply use the df you have already created
#1.
gr_obj =  read.table("/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/allPeaks.bed", header = F, sep = '\t')

#2.
#gr_obj <- merged_peak_df
gr_obj <- toGRanges(gr_obj)

peakAnno <- annotatePeak(peak = gr_obj, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
plotAnnoPie(peakAnno)

plotdist<- plotDistToTSS(peakAnno,
              title="Distribution of Peaks to TSS")

plotdist

dist_df <- as.data.frame(plotdist[["data"]])
write.table(dist_df, file="/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/PlotDist.txt", quote=F, sep="\t", row.names=F, col.names=T)

annotation_df<-as.data.frame(peakAnno@anno@elementMetadata@listData[["annotation"]])


write.table(annotation_df, file="/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/annotation_df.txt", quote=F, sep="\t", row.names=F, col.names=T)
