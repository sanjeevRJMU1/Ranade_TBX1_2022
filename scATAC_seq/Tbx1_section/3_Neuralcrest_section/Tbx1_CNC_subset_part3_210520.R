# 
markerList_PA3_Cardiac_NeuralCrest_E115 <- getMarkers(markerTest_PA3_Cardiac_NeuralCrest_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_PA3_Cardiac_NeuralCrest_E115_df <- as.data.frame(markerList_PA3_Cardiac_NeuralCrest_E115@listData[["PA3_Cardiac_NeuralCrest_15_W"]])
markerList_PA3_Cardiac_NeuralCrest_E115_df_homer <- markerList_PA3_Cardiac_NeuralCrest_E115_df[order(markerList_PA3_Cardiac_NeuralCrest_E115_df$Log2FC, decreasing = TRUE),]
markerList_PA3_Cardiac_NeuralCrest_E115_df_homer$rank <- seq_len(nrow(markerList_PA3_Cardiac_NeuralCrest_E115_df_homer))
markerList_PA3_Cardiac_NeuralCrest_E115_df_homer$blank <- ""
markerList_PA3_Cardiac_NeuralCrest_E115_df_homer_trim <- markerList_PA3_Cardiac_NeuralCrest_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_PA3_Cardiac_NeuralCrest_E115_df_homer_trim, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/markerList_PA3_Cardiac_NeuralCrest_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_PA3_Cardiac_NeuralCrest_E115_up <- getMarkers(markerTest_PA3_Cardiac_NeuralCrest_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_PA3_Cardiac_NeuralCrest_E115_up_df <- as.data.frame(markerList_PA3_Cardiac_NeuralCrest_E115_up@listData[["PA3_Cardiac_NeuralCrest_15_W"]])
markerList_PA3_Cardiac_NeuralCrest_E115_up_df_homer <- markerList_PA3_Cardiac_NeuralCrest_E115_up_df[order(markerList_PA3_Cardiac_NeuralCrest_E115_up_df$Log2FC, decreasing = F),]
markerList_PA3_Cardiac_NeuralCrest_E115_up_df_homer$rank <- seq_len(nrow(markerList_PA3_Cardiac_NeuralCrest_E115_up_df_homer))
markerList_PA3_Cardiac_NeuralCrest_E115_up_df_homer$blank <-""
markerList_PA3_Cardiac_NeuralCrest_E115_up_df_homer_trim <- markerList_PA3_Cardiac_NeuralCrest_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_PA3_Cardiac_NeuralCrest_E115_up_df_homer_trim, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/markerList_PA3_Cardiac_NeuralCrest_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_Craniofacial_NeuralCrest_E115 <- getMarkers(markerTest_Craniofacial_NeuralCrest_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Craniofacial_NeuralCrest_E115_df <- as.data.frame(markerList_Craniofacial_NeuralCrest_E115@listData[["Craniofacial_NeuralCrest_15_W"]])
markerList_Craniofacial_NeuralCrest_E115_df_homer <- markerList_Craniofacial_NeuralCrest_E115_df[order(markerList_Craniofacial_NeuralCrest_E115_df$Log2FC, decreasing = TRUE),]
markerList_Craniofacial_NeuralCrest_E115_df_homer$rank <- seq_len(nrow(markerList_Craniofacial_NeuralCrest_E115_df_homer))
markerList_Craniofacial_NeuralCrest_E115_df_homer$blank <- ""
markerList_Craniofacial_NeuralCrest_E115_df_homer_trim <- markerList_Craniofacial_NeuralCrest_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Craniofacial_NeuralCrest_E115_df_homer_trim, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/markerList_Craniofacial_NeuralCrest_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Craniofacial_NeuralCrest_E115_up <- getMarkers(markerTest_Craniofacial_NeuralCrest_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Craniofacial_NeuralCrest_E115_up_df <- as.data.frame(markerList_Craniofacial_NeuralCrest_E115_up@listData[["Craniofacial_NeuralCrest_15_W"]])
markerList_Craniofacial_NeuralCrest_E115_up_df_homer <- markerList_Craniofacial_NeuralCrest_E115_up_df[order(markerList_Craniofacial_NeuralCrest_E115_up_df$Log2FC, decreasing = F),]
markerList_Craniofacial_NeuralCrest_E115_up_df_homer$rank <- seq_len(nrow(markerList_Craniofacial_NeuralCrest_E115_up_df_homer))
markerList_Craniofacial_NeuralCrest_E115_up_df_homer$blank <-""
markerList_Craniofacial_NeuralCrest_E115_up_df_homer_trim <- markerList_Craniofacial_NeuralCrest_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Craniofacial_NeuralCrest_E115_up_df_homer_trim, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/markerList_Craniofacial_NeuralCrest_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Craniofacial_NeuralCrest_E105 <- getMarkers(markerTest_Craniofacial_NeuralCrest_E105, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Craniofacial_NeuralCrest_E105_df <- as.data.frame(markerList_Craniofacial_NeuralCrest_E105@listData[["Craniofacial_NeuralCrest_05_W"]])
markerList_Craniofacial_NeuralCrest_E105_df_homer <- markerList_Craniofacial_NeuralCrest_E105_df[order(markerList_Craniofacial_NeuralCrest_E105_df$Log2FC, decreasing = TRUE),]
markerList_Craniofacial_NeuralCrest_E105_df_homer$rank <- seq_len(nrow(markerList_Craniofacial_NeuralCrest_E105_df_homer))
markerList_Craniofacial_NeuralCrest_E105_df_homer$blank <- ""
markerList_Craniofacial_NeuralCrest_E105_df_homer_trim <- markerList_Craniofacial_NeuralCrest_E105_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Craniofacial_NeuralCrest_E105_df_homer_trim, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/markerList_Craniofacial_NeuralCrest_E105_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Craniofacial_NeuralCrest_E105_up <- getMarkers(markerTest_Craniofacial_NeuralCrest_E105, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Craniofacial_NeuralCrest_E105_up_df <- as.data.frame(markerList_Craniofacial_NeuralCrest_E105_up@listData[["Craniofacial_NeuralCrest_05_W"]])
markerList_Craniofacial_NeuralCrest_E105_up_df_homer <- markerList_Craniofacial_NeuralCrest_E105_up_df[order(markerList_Craniofacial_NeuralCrest_E105_up_df$Log2FC, decreasing = F),]
markerList_Craniofacial_NeuralCrest_E105_up_df_homer$rank <- seq_len(nrow(markerList_Craniofacial_NeuralCrest_E105_up_df_homer))
markerList_Craniofacial_NeuralCrest_E105_up_df_homer$blank <-""
markerList_Craniofacial_NeuralCrest_E105_up_df_homer_trim <- markerList_Craniofacial_NeuralCrest_E105_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Craniofacial_NeuralCrest_E105_up_df_homer_trim, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/markerList_Craniofacial_NeuralCrest_E105_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Cardiac_NeuralCrest_E115 <- getMarkers(markerTest_Cardiac_NeuralCrest_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Cardiac_NeuralCrest_E115_df <- as.data.frame(markerList_Cardiac_NeuralCrest_E115@listData[["Cardiac_NeuralCrest_15_W"]])
markerList_Cardiac_NeuralCrest_E115_df_homer <- markerList_Cardiac_NeuralCrest_E115_df[order(markerList_Cardiac_NeuralCrest_E115_df$Log2FC, decreasing = TRUE),]
markerList_Cardiac_NeuralCrest_E115_df_homer$rank <- seq_len(nrow(markerList_Cardiac_NeuralCrest_E115_df_homer))
markerList_Cardiac_NeuralCrest_E115_df_homer$blank <- ""
markerList_Cardiac_NeuralCrest_E115_df_homer_trim <- markerList_Cardiac_NeuralCrest_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiac_NeuralCrest_E115_df_homer_trim, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/markerList_Cardiac_NeuralCrest_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Cardiac_NeuralCrest_E115_up <- getMarkers(markerTest_Cardiac_NeuralCrest_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Cardiac_NeuralCrest_E115_up_df <- as.data.frame(markerList_Cardiac_NeuralCrest_E115_up@listData[["Cardiac_NeuralCrest_15_W"]])
markerList_Cardiac_NeuralCrest_E115_up_df_homer <- markerList_Cardiac_NeuralCrest_E115_up_df[order(markerList_Cardiac_NeuralCrest_E115_up_df$Log2FC, decreasing = F),]
markerList_Cardiac_NeuralCrest_E115_up_df_homer$rank <- seq_len(nrow(markerList_Cardiac_NeuralCrest_E115_up_df_homer))
markerList_Cardiac_NeuralCrest_E115_up_df_homer$blank <-""
markerList_Cardiac_NeuralCrest_E115_up_df_homer_trim <- markerList_Cardiac_NeuralCrest_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiac_NeuralCrest_E115_up_df_homer_trim, file = "/Users/sranade/scATAC-seq/CNC_subset_210607/Exports/markerList_Cardiac_NeuralCrest_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)
