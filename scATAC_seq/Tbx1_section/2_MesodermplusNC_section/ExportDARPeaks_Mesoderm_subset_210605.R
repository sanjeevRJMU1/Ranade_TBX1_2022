setwd("/Users/sranade/scATAC-seq/Mesoderm_subset3_210517/Exports/DAR_exports_210110")
# 
markerList_Anterior_Second_Heart_Field_E925 <- getMarkers(markerTest_Anterior_Second_Heart_Field_E925, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Anterior_Second_Heart_Field_E925_df <- as.data.frame(markerList_Anterior_Second_Heart_Field_E925@listData[["Anterior_Second_Heart_Field_25_W"]])
markerList_Anterior_Second_Heart_Field_E925_df_homer <- markerList_Anterior_Second_Heart_Field_E925_df[order(markerList_Anterior_Second_Heart_Field_E925_df$Log2FC, decreasing = TRUE),]
markerList_Anterior_Second_Heart_Field_E925_df_homer$rank <- seq_len(nrow(markerList_Anterior_Second_Heart_Field_E925_df_homer))
markerList_Anterior_Second_Heart_Field_E925_df_homer$blank <- ""
markerList_Anterior_Second_Heart_Field_E925_df_homer_trim <- markerList_Anterior_Second_Heart_Field_E925_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Anterior_Second_Heart_Field_E925_df_homer_trim, file = "markerList_Anterior_Second_Heart_Field_E925_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Anterior_Second_Heart_Field_E925_up <- getMarkers(markerTest_Anterior_Second_Heart_Field_E925, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Anterior_Second_Heart_Field_E925_up_df <- as.data.frame(markerList_Anterior_Second_Heart_Field_E925_up@listData[["Anterior_Second_Heart_Field_25_W"]])
markerList_Anterior_Second_Heart_Field_E925_up_df_homer <- markerList_Anterior_Second_Heart_Field_E925_up_df[order(markerList_Anterior_Second_Heart_Field_E925_up_df$Log2FC, decreasing = F),]
markerList_Anterior_Second_Heart_Field_E925_up_df_homer$rank <- seq_len(nrow(markerList_Anterior_Second_Heart_Field_E925_up_df_homer))
markerList_Anterior_Second_Heart_Field_E925_up_df_homer$blank <-""
markerList_Anterior_Second_Heart_Field_E925_up_df_homer_trim <- markerList_Anterior_Second_Heart_Field_E925_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Anterior_Second_Heart_Field_E925_up_df_homer_trim, file = "markerList_Anterior_Second_Heart_Field_E925_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_Anterior_Second_Heart_Field_E105 <- getMarkers(markerTest_Anterior_Second_Heart_Field_E105, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Anterior_Second_Heart_Field_E105_df <- as.data.frame(markerList_Anterior_Second_Heart_Field_E105@listData[["Anterior_Second_Heart_Field_05_W"]])
markerList_Anterior_Second_Heart_Field_E105_df_homer <- markerList_Anterior_Second_Heart_Field_E105_df[order(markerList_Anterior_Second_Heart_Field_E105_df$Log2FC, decreasing = TRUE),]
markerList_Anterior_Second_Heart_Field_E105_df_homer$rank <- seq_len(nrow(markerList_Anterior_Second_Heart_Field_E105_df_homer))
markerList_Anterior_Second_Heart_Field_E105_df_homer$blank <- ""
markerList_Anterior_Second_Heart_Field_E105_df_homer_trim <- markerList_Anterior_Second_Heart_Field_E105_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Anterior_Second_Heart_Field_E105_df_homer_trim, file = "markerList_Anterior_Second_Heart_Field_E105_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Anterior_Second_Heart_Field_E105_up <- getMarkers(markerTest_Anterior_Second_Heart_Field_E105, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Anterior_Second_Heart_Field_E105_up_df <- as.data.frame(markerList_Anterior_Second_Heart_Field_E105_up@listData[["Anterior_Second_Heart_Field_05_W"]])
markerList_Anterior_Second_Heart_Field_E105_up_df_homer <- markerList_Anterior_Second_Heart_Field_E105_up_df[order(markerList_Anterior_Second_Heart_Field_E105_up_df$Log2FC, decreasing = F),]
markerList_Anterior_Second_Heart_Field_E105_up_df_homer$rank <- seq_len(nrow(markerList_Anterior_Second_Heart_Field_E105_up_df_homer))
markerList_Anterior_Second_Heart_Field_E105_up_df_homer$blank <-""
markerList_Anterior_Second_Heart_Field_E105_up_df_homer_trim <- markerList_Anterior_Second_Heart_Field_E105_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Anterior_Second_Heart_Field_E105_up_df_homer_trim, file = "markerList_Anterior_Second_Heart_Field_E105_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Anterior_Second_Heart_Field_E115 <- getMarkers(markerTest_Anterior_Second_Heart_Field_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Anterior_Second_Heart_Field_E115_df <- as.data.frame(markerList_Anterior_Second_Heart_Field_E115@listData[["Anterior_Second_Heart_Field_15_W"]])
markerList_Anterior_Second_Heart_Field_E115_df_homer <- markerList_Anterior_Second_Heart_Field_E115_df[order(markerList_Anterior_Second_Heart_Field_E115_df$Log2FC, decreasing = TRUE),]
markerList_Anterior_Second_Heart_Field_E115_df_homer$rank <- seq_len(nrow(markerList_Anterior_Second_Heart_Field_E115_df_homer))
markerList_Anterior_Second_Heart_Field_E115_df_homer$blank <- ""
markerList_Anterior_Second_Heart_Field_E115_df_homer_trim <- markerList_Anterior_Second_Heart_Field_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Anterior_Second_Heart_Field_E115_df_homer_trim, file = "markerList_Anterior_Second_Heart_Field_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Anterior_Second_Heart_Field_E115_up <- getMarkers(markerTest_Anterior_Second_Heart_Field_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Anterior_Second_Heart_Field_E115_up_df <- as.data.frame(markerList_Anterior_Second_Heart_Field_E115_up@listData[["Anterior_Second_Heart_Field_15_W"]])
markerList_Anterior_Second_Heart_Field_E115_up_df_homer <- markerList_Anterior_Second_Heart_Field_E115_up_df[order(markerList_Anterior_Second_Heart_Field_E115_up_df$Log2FC, decreasing = F),]
markerList_Anterior_Second_Heart_Field_E115_up_df_homer$rank <- seq_len(nrow(markerList_Anterior_Second_Heart_Field_E115_up_df_homer))
markerList_Anterior_Second_Heart_Field_E115_up_df_homer$blank <-""
markerList_Anterior_Second_Heart_Field_E115_up_df_homer_trim <- markerList_Anterior_Second_Heart_Field_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Anterior_Second_Heart_Field_E115_up_df_homer_trim, file = "markerList_Anterior_Second_Heart_Field_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

#############
markerList_Cardiopharyngeal_Mesoderm_E925 <- getMarkers(markerTest_Cardiopharyngeal_Mesoderm_E925, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Cardiopharyngeal_Mesoderm_E925_df <- as.data.frame(markerList_Cardiopharyngeal_Mesoderm_E925@listData[["Cardiopharyngeal_Mesoderm_25_W"]])
markerList_Cardiopharyngeal_Mesoderm_E925_df_homer <- markerList_Cardiopharyngeal_Mesoderm_E925_df[order(markerList_Cardiopharyngeal_Mesoderm_E925_df$Log2FC, decreasing = TRUE),]
markerList_Cardiopharyngeal_Mesoderm_E925_df_homer$rank <- seq_len(nrow(markerList_Cardiopharyngeal_Mesoderm_E925_df_homer))
markerList_Cardiopharyngeal_Mesoderm_E925_df_homer$blank <- ""
markerList_Cardiopharyngeal_Mesoderm_E925_df_homer_trim <- markerList_Cardiopharyngeal_Mesoderm_E925_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiopharyngeal_Mesoderm_E925_df_homer_trim, file = "markerList_Cardiopharyngeal_Mesoderm_E925_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Cardiopharyngeal_Mesoderm_E925_up <- getMarkers(markerTest_Cardiopharyngeal_Mesoderm_E925, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Cardiopharyngeal_Mesoderm_E925_up_df <- as.data.frame(markerList_Cardiopharyngeal_Mesoderm_E925_up@listData[["Cardiopharyngeal_Mesoderm_25_W"]])
markerList_Cardiopharyngeal_Mesoderm_E925_up_df_homer <- markerList_Cardiopharyngeal_Mesoderm_E925_up_df[order(markerList_Cardiopharyngeal_Mesoderm_E925_up_df$Log2FC, decreasing = F),]
markerList_Cardiopharyngeal_Mesoderm_E925_up_df_homer$rank <- seq_len(nrow(markerList_Cardiopharyngeal_Mesoderm_E925_up_df_homer))
markerList_Cardiopharyngeal_Mesoderm_E925_up_df_homer$blank <-""
markerList_Cardiopharyngeal_Mesoderm_E925_up_df_homer_trim <- markerList_Cardiopharyngeal_Mesoderm_E925_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiopharyngeal_Mesoderm_E925_up_df_homer_trim, file = "markerList_Cardiopharyngeal_Mesoderm_E925_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

# 
markerList_Cardiopharyngeal_Mesoderm_E105 <- getMarkers(markerTest_Cardiopharyngeal_Mesoderm_E105, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Cardiopharyngeal_Mesoderm_E105_df <- as.data.frame(markerList_Cardiopharyngeal_Mesoderm_E105@listData[["Cardiopharyngeal_Mesoderm_05_W"]])
markerList_Cardiopharyngeal_Mesoderm_E105_df_homer <- markerList_Cardiopharyngeal_Mesoderm_E105_df[order(markerList_Cardiopharyngeal_Mesoderm_E105_df$Log2FC, decreasing = TRUE),]
markerList_Cardiopharyngeal_Mesoderm_E105_df_homer$rank <- seq_len(nrow(markerList_Cardiopharyngeal_Mesoderm_E105_df_homer))
markerList_Cardiopharyngeal_Mesoderm_E105_df_homer$blank <- ""
markerList_Cardiopharyngeal_Mesoderm_E105_df_homer_trim <- markerList_Cardiopharyngeal_Mesoderm_E105_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiopharyngeal_Mesoderm_E105_df_homer_trim, file = "markerList_Cardiopharyngeal_Mesoderm_E105_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Cardiopharyngeal_Mesoderm_E105_up <- getMarkers(markerTest_Cardiopharyngeal_Mesoderm_E105, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Cardiopharyngeal_Mesoderm_E105_up_df <- as.data.frame(markerList_Cardiopharyngeal_Mesoderm_E105_up@listData[["Cardiopharyngeal_Mesoderm_05_W"]])
markerList_Cardiopharyngeal_Mesoderm_E105_up_df_homer <- markerList_Cardiopharyngeal_Mesoderm_E105_up_df[order(markerList_Cardiopharyngeal_Mesoderm_E105_up_df$Log2FC, decreasing = F),]
markerList_Cardiopharyngeal_Mesoderm_E105_up_df_homer$rank <- seq_len(nrow(markerList_Cardiopharyngeal_Mesoderm_E105_up_df_homer))
markerList_Cardiopharyngeal_Mesoderm_E105_up_df_homer$blank <-""
markerList_Cardiopharyngeal_Mesoderm_E105_up_df_homer_trim <- markerList_Cardiopharyngeal_Mesoderm_E105_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiopharyngeal_Mesoderm_E105_up_df_homer_trim, file = "markerList_Cardiopharyngeal_Mesoderm_E105_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Cardiopharyngeal_Mesoderm_E115 <- getMarkers(markerTest_Cardiopharyngeal_Mesoderm_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Cardiopharyngeal_Mesoderm_E115_df <- as.data.frame(markerList_Cardiopharyngeal_Mesoderm_E115@listData[["Cardiopharyngeal_Mesoderm_15_W"]])
markerList_Cardiopharyngeal_Mesoderm_E115_df_homer <- markerList_Cardiopharyngeal_Mesoderm_E115_df[order(markerList_Cardiopharyngeal_Mesoderm_E115_df$Log2FC, decreasing = TRUE),]
markerList_Cardiopharyngeal_Mesoderm_E115_df_homer$rank <- seq_len(nrow(markerList_Cardiopharyngeal_Mesoderm_E115_df_homer))
markerList_Cardiopharyngeal_Mesoderm_E115_df_homer$blank <- ""
markerList_Cardiopharyngeal_Mesoderm_E115_df_homer_trim <- markerList_Cardiopharyngeal_Mesoderm_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiopharyngeal_Mesoderm_E115_df_homer_trim, file = "markerList_Cardiopharyngeal_Mesoderm_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Cardiopharyngeal_Mesoderm_E115_up <- getMarkers(markerTest_Cardiopharyngeal_Mesoderm_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Cardiopharyngeal_Mesoderm_E115_up_df <- as.data.frame(markerList_Cardiopharyngeal_Mesoderm_E115_up@listData[["Cardiopharyngeal_Mesoderm_15_W"]])
markerList_Cardiopharyngeal_Mesoderm_E115_up_df_homer <- markerList_Cardiopharyngeal_Mesoderm_E115_up_df[order(markerList_Cardiopharyngeal_Mesoderm_E115_up_df$Log2FC, decreasing = F),]
markerList_Cardiopharyngeal_Mesoderm_E115_up_df_homer$rank <- seq_len(nrow(markerList_Cardiopharyngeal_Mesoderm_E115_up_df_homer))
markerList_Cardiopharyngeal_Mesoderm_E115_up_df_homer$blank <-""
markerList_Cardiopharyngeal_Mesoderm_E115_up_df_homer_trim <- markerList_Cardiopharyngeal_Mesoderm_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiopharyngeal_Mesoderm_E115_up_df_homer_trim, file = "markerList_Cardiopharyngeal_Mesoderm_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)
##############
# 
markerList_Cardiomyocyte_E115 <- getMarkers(markerTest_Cardiomyocyte_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Cardiomyocyte_E115_df <- as.data.frame(markerList_Cardiomyocyte_E115@listData[["Cardiomyocyte_15_W"]])
markerList_Cardiomyocyte_E115_df_homer <- markerList_Cardiomyocyte_E115_df[order(markerList_Cardiomyocyte_E115_df$Log2FC, decreasing = TRUE),]
markerList_Cardiomyocyte_E115_df_homer$rank <- seq_len(nrow(markerList_Cardiomyocyte_E115_df_homer))
markerList_Cardiomyocyte_E115_df_homer$blank <- ""
markerList_Cardiomyocyte_E115_df_homer_trim <- markerList_Cardiomyocyte_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiomyocyte_E115_df_homer_trim, file = "markerList_Cardiomyocyte_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Cardiomyocyte_E115_up <- getMarkers(markerTest_Cardiomyocyte_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Cardiomyocyte_E115_up_df <- as.data.frame(markerList_Cardiomyocyte_E115_up@listData[["Cardiomyocyte_15_W"]])
markerList_Cardiomyocyte_E115_up_df_homer <- markerList_Cardiomyocyte_E115_up_df[order(markerList_Cardiomyocyte_E115_up_df$Log2FC, decreasing = F),]
markerList_Cardiomyocyte_E115_up_df_homer$rank <- seq_len(nrow(markerList_Cardiomyocyte_E115_up_df_homer))
markerList_Cardiomyocyte_E115_up_df_homer$blank <-""
markerList_Cardiomyocyte_E115_up_df_homer_trim <- markerList_Cardiomyocyte_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cardiomyocyte_E115_up_df_homer_trim, file = "markerList_Cardiomyocyte_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

##############
# 
markerList_Neural_Crest_E105 <- getMarkers(markerTest_Neural_Crest_E105, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Neural_Crest_E105_df <- as.data.frame(markerList_Neural_Crest_E105@listData[["Neural_Crest_05_W"]])
markerList_Neural_Crest_E105_df_homer <- markerList_Neural_Crest_E105_df[order(markerList_Neural_Crest_E105_df$Log2FC, decreasing = TRUE),]
markerList_Neural_Crest_E105_df_homer$rank <- seq_len(nrow(markerList_Neural_Crest_E105_df_homer))
markerList_Neural_Crest_E105_df_homer$blank <- ""
markerList_Neural_Crest_E105_df_homer_trim <- markerList_Neural_Crest_E105_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Neural_Crest_E105_df_homer_trim, file = "markerList_Neural_Crest_E105_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Neural_Crest_E105_up <- getMarkers(markerTest_Neural_Crest_E105, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Neural_Crest_E105_up_df <- as.data.frame(markerList_Neural_Crest_E105_up@listData[["Neural_Crest_05_W"]])
markerList_Neural_Crest_E105_up_df_homer <- markerList_Neural_Crest_E105_up_df[order(markerList_Neural_Crest_E105_up_df$Log2FC, decreasing = F),]
markerList_Neural_Crest_E105_up_df_homer$rank <- seq_len(nrow(markerList_Neural_Crest_E105_up_df_homer))
markerList_Neural_Crest_E105_up_df_homer$blank <-""
markerList_Neural_Crest_E105_up_df_homer_trim <- markerList_Neural_Crest_E105_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Neural_Crest_E105_up_df_homer_trim, file = "markerList_Neural_Crest_E105_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Neural_Crest_E115 <- getMarkers(markerTest_Neural_Crest_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Neural_Crest_E115_df <- as.data.frame(markerList_Neural_Crest_E115@listData[["Neural_Crest_15_W"]])
markerList_Neural_Crest_E115_df_homer <- markerList_Neural_Crest_E115_df[order(markerList_Neural_Crest_E115_df$Log2FC, decreasing = TRUE),]
markerList_Neural_Crest_E115_df_homer$rank <- seq_len(nrow(markerList_Neural_Crest_E115_df_homer))
markerList_Neural_Crest_E115_df_homer$blank <- ""
markerList_Neural_Crest_E115_df_homer_trim <- markerList_Neural_Crest_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Neural_Crest_E115_df_homer_trim, file = "markerList_Neural_Crest_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Neural_Crest_E115_up <- getMarkers(markerTest_Neural_Crest_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Neural_Crest_E115_up_df <- as.data.frame(markerList_Neural_Crest_E115_up@listData[["Neural_Crest_15_W"]])
markerList_Neural_Crest_E115_up_df_homer <- markerList_Neural_Crest_E115_up_df[order(markerList_Neural_Crest_E115_up_df$Log2FC, decreasing = F),]
markerList_Neural_Crest_E115_up_df_homer$rank <- seq_len(nrow(markerList_Neural_Crest_E115_up_df_homer))
markerList_Neural_Crest_E115_up_df_homer$blank <-""
markerList_Neural_Crest_E115_up_df_homer_trim <- markerList_Neural_Crest_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Neural_Crest_E115_up_df_homer_trim, file = "markerList_Neural_Crest_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

##############
# 
markerList_Cranial_Mesenchyme_E115 <- getMarkers(markerTest_Cranial_Mesenchyme_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Cranial_Mesenchyme_E115_df <- as.data.frame(markerList_Cranial_Mesenchyme_E115@listData[["Cranial_Mesenchyme_15_W"]])
markerList_Cranial_Mesenchyme_E115_df_homer <- markerList_Cranial_Mesenchyme_E115_df[order(markerList_Cranial_Mesenchyme_E115_df$Log2FC, decreasing = TRUE),]
markerList_Cranial_Mesenchyme_E115_df_homer$rank <- seq_len(nrow(markerList_Cranial_Mesenchyme_E115_df_homer))
markerList_Cranial_Mesenchyme_E115_df_homer$blank <- ""
markerList_Cranial_Mesenchyme_E115_df_homer_trim <- markerList_Cranial_Mesenchyme_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cranial_Mesenchyme_E115_df_homer_trim, file = "markerList_Cranial_Mesenchyme_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Cranial_Mesenchyme_E115_up <- getMarkers(markerTest_Cranial_Mesenchyme_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Cranial_Mesenchyme_E115_up_df <- as.data.frame(markerList_Cranial_Mesenchyme_E115_up@listData[["Cranial_Mesenchyme_15_W"]])
markerList_Cranial_Mesenchyme_E115_up_df_homer <- markerList_Cranial_Mesenchyme_E115_up_df[order(markerList_Cranial_Mesenchyme_E115_up_df$Log2FC, decreasing = F),]
markerList_Cranial_Mesenchyme_E115_up_df_homer$rank <- seq_len(nrow(markerList_Cranial_Mesenchyme_E115_up_df_homer))
markerList_Cranial_Mesenchyme_E115_up_df_homer$blank <-""
markerList_Cranial_Mesenchyme_E115_up_df_homer_trim <- markerList_Cranial_Mesenchyme_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Cranial_Mesenchyme_E115_up_df_homer_trim, file = "markerList_Cranial_Mesenchyme_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)

##############
# 
markerList_Paraxial_Mesoderm_E115 <- getMarkers(markerTest_Paraxial_Mesoderm_E115, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Paraxial_Mesoderm_E115_df <- as.data.frame(markerList_Paraxial_Mesoderm_E115@listData[["Paraxial_Mesoderm_15_W"]])
markerList_Paraxial_Mesoderm_E115_df_homer <- markerList_Paraxial_Mesoderm_E115_df[order(markerList_Paraxial_Mesoderm_E115_df$Log2FC, decreasing = TRUE),]
markerList_Paraxial_Mesoderm_E115_df_homer$rank <- seq_len(nrow(markerList_Paraxial_Mesoderm_E115_df_homer))
markerList_Paraxial_Mesoderm_E115_df_homer$blank <- ""
markerList_Paraxial_Mesoderm_E115_df_homer_trim <- markerList_Paraxial_Mesoderm_E115_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Paraxial_Mesoderm_E115_df_homer_trim, file = "markerList_Paraxial_Mesoderm_E115_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Paraxial_Mesoderm_E115_up <- getMarkers(markerTest_Paraxial_Mesoderm_E115, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Paraxial_Mesoderm_E115_up_df <- as.data.frame(markerList_Paraxial_Mesoderm_E115_up@listData[["Paraxial_Mesoderm_15_W"]])
markerList_Paraxial_Mesoderm_E115_up_df_homer <- markerList_Paraxial_Mesoderm_E115_up_df[order(markerList_Paraxial_Mesoderm_E115_up_df$Log2FC, decreasing = F),]
markerList_Paraxial_Mesoderm_E115_up_df_homer$rank <- seq_len(nrow(markerList_Paraxial_Mesoderm_E115_up_df_homer))
markerList_Paraxial_Mesoderm_E115_up_df_homer$blank <-""
markerList_Paraxial_Mesoderm_E115_up_df_homer_trim <- markerList_Paraxial_Mesoderm_E115_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Paraxial_Mesoderm_E115_up_df_homer_trim, file = "markerList_Paraxial_Mesoderm_E115_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)
