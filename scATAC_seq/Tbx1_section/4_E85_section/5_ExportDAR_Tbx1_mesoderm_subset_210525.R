setwd("/Users/sranade/scATAC-seq/Tbx1_E825_Mesoderm_subset_210524/Exports/mesoderm_subset_DAR")
# 
markerList_C5 <- getMarkers(markerTest_C5, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_C5_df <- as.data.frame(markerList_C5@listData[["C5_W"]])
markerList_C5_df_homer <- markerList_C5_df[order(markerList_C5_df$Log2FC, decreasing = TRUE),]
markerList_C5_df_homer$rank <- seq_len(nrow(markerList_C5_df_homer))
markerList_C5_df_homer$blank <- ""
markerList_C5_df_homer_trim <- markerList_C5_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C5_df_homer_trim, file = "markerList_C5_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_C5_up <- getMarkers(markerTest_C5, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_C5_up_df <- as.data.frame(markerList_C5_up@listData[["C5_W"]])
markerList_C5_up_df_homer <- markerList_C5_up_df[order(markerList_C5_up_df$Log2FC, decreasing = F),]
markerList_C5_up_df_homer$rank <- seq_len(nrow(markerList_C5_up_df_homer))
markerList_C5_up_df_homer$blank <-""
markerList_C5_up_df_homer_trim <- markerList_C5_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C5_up_df_homer_trim, file = "markerList_C5_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_C10 <- getMarkers(markerTest_C10, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_C10_df <- as.data.frame(markerList_C10@listData[["C10_W"]])
markerList_C10_df_homer <- markerList_C10_df[order(markerList_C10_df$Log2FC, decreasing = TRUE),]
markerList_C10_df_homer$rank <- seq_len(nrow(markerList_C10_df_homer))
markerList_C10_df_homer$blank <- ""
markerList_C10_df_homer_trim <- markerList_C10_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C10_df_homer_trim, file = "markerList_C10_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_C10_up <- getMarkers(markerTest_C10, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_C10_up_df <- as.data.frame(markerList_C10_up@listData[["C10_W"]])
markerList_C10_up_df_homer <- markerList_C10_up_df[order(markerList_C10_up_df$Log2FC, decreasing = F),]
markerList_C10_up_df_homer$rank <- seq_len(nrow(markerList_C10_up_df_homer))
markerList_C10_up_df_homer$blank <-""
markerList_C10_up_df_homer_trim <- markerList_C10_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C10_up_df_homer_trim, file = "markerList_C10_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_C11 <- getMarkers(markerTest_C11, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_C11_df <- as.data.frame(markerList_C11@listData[["C11_W"]])
markerList_C11_df_homer <- markerList_C11_df[order(markerList_C11_df$Log2FC, decreasing = TRUE),]
markerList_C11_df_homer$rank <- seq_len(nrow(markerList_C11_df_homer))
markerList_C11_df_homer$blank <- ""
markerList_C11_df_homer_trim <- markerList_C11_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C11_df_homer_trim, file = "markerList_C11_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_C11_up <- getMarkers(markerTest_C11, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_C11_up_df <- as.data.frame(markerList_C11_up@listData[["C11_W"]])
markerList_C11_up_df_homer <- markerList_C11_up_df[order(markerList_C11_up_df$Log2FC, decreasing = F),]
markerList_C11_up_df_homer$rank <- seq_len(nrow(markerList_C11_up_df_homer))
markerList_C11_up_df_homer$blank <-""
markerList_C11_up_df_homer_trim <- markerList_C11_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C11_up_df_homer_trim, file = "markerList_C11_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_C12 <- getMarkers(markerTest_C12, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_C12_df <- as.data.frame(markerList_C12@listData[["C12_W"]])
markerList_C12_df_homer <- markerList_C12_df[order(markerList_C12_df$Log2FC, decreasing = TRUE),]
markerList_C12_df_homer$rank <- seq_len(nrow(markerList_C12_df_homer))
markerList_C12_df_homer$blank <- ""
markerList_C12_df_homer_trim <- markerList_C12_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C12_df_homer_trim, file = "markerList_C12_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_C12_up <- getMarkers(markerTest_C12, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_C12_up_df <- as.data.frame(markerList_C12_up@listData[["C12_W"]])
markerList_C12_up_df_homer <- markerList_C12_up_df[order(markerList_C12_up_df$Log2FC, decreasing = F),]
markerList_C12_up_df_homer$rank <- seq_len(nrow(markerList_C12_up_df_homer))
markerList_C12_up_df_homer$blank <-""
markerList_C12_up_df_homer_trim <- markerList_C12_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_C12_up_df_homer_trim, file = "markerList_C12_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)