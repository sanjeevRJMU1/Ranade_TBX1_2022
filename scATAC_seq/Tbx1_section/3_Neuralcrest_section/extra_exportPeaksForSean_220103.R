#### exports for sean, redo 220103

# look at the peak set and then save it as a separate object
getPeakSet(CNC_subset)
#GRanges object with 303570 ranges and 13 metadata columns:
merged_peak_setGR <- getPeakSet(CNC_subset)



## Export all non-overlapping 303570 peaks
setwd("/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/forsean")
merged_peak_df <- as.data.frame(merged_peak_setGR@ranges)

merged_peak_df <- data.frame(seqnames=seqnames(merged_peak_setGR),
                             starts=start(merged_peak_setGR)-1,
                             ends=end(merged_peak_setGR),
                             names=c(rep(".", length(merged_peak_setGR))),
                             scores=merged_peak_setGR@elementMetadata@listData[["score"]],
                             strands=merged_peak_setGR@strand@values)

# format in same manner as you would use for homer 
merged_peak_df <- merged_peak_df %>% mutate(id = row_number())
merged_peak_df$blankVar <- ''
merged_peak_df_homerformat <- merged_peak_df[,c(1,2,3,7,8,6)]
write.table(merged_peak_df_homerformat, file="allPeaks_220103.bed", quote=F, sep="\t", row.names=F, col.names=F)
