### Part 1: Export Peaks Enriched in WT v KO (ie, down-regulated in KO) for clusters of interest

markerList_AHF_df_homer <- markerList_AHF_df[order(markerList_AHF_df$Log2FC, decreasing = TRUE),]
markerList_AHF_df_homer$rank <- seq_len(nrow(markerList_AHF_df_homer))
markerList_AHF_df_homer_trim <- markerList_AHF_df_homer[,c(1,2,3,9,5)]
write.table(markerList_AHF_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/markerList_AHF_df_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)

markerList_PharyngealMesoderm_df_homer <- markerList_PharyngealMesoderm_df[order(markerList_PharyngealMesoderm_df$Log2FC, decreasing = TRUE),]
markerList_PharyngealMesoderm_df_homer$rank <- seq_len(nrow(markerList_PharyngealMesoderm_df_homer))
markerList_PharyngealMesoderm_df_homer_trim <- markerList_PharyngealMesoderm_df_homer[,c(1,2,3,9,5)]
write.table(markerList_PharyngealMesoderm_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/markerList_PharyngealMesoderm_df_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)

markerList_Endoderm_df_homer <- markerList_Endoderm_df[order(markerList_Endoderm_df$Log2FC, decreasing = TRUE),]
markerList_Endoderm_df_homer$rank <- seq_len(nrow(markerList_Endoderm_df_homer))
markerList_Endoderm_df_homer_trim <- markerList_Endoderm_df_homer[,c(1,2,3,9,5)]
write.table(markerList_Endoderm_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/markerList_Endoderm_df_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)

markerList_Ectoderm_df_homer <- markerList_Ectoderm_df[order(markerList_Ectoderm_df$Log2FC, decreasing = TRUE),]
markerList_Ectoderm_df_homer$rank <- seq_len(nrow(markerList_Ectoderm_df_homer))
markerList_Ectoderm_df_homer_trim <- markerList_Ectoderm_df_homer[,c(1,2,3,9,5)]
write.table(markerList_Ectoderm_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/markerList_Ectoderm_df_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)


markerList_PharyngealMesoderm_up_df_homer <- markerList_PharyngealMesoderm_up_df[order(markerList_PharyngealMesoderm_up_df$Log2FC, decreasing = TRUE),]
markerList_PharyngealMesoderm_up_df_homer$rank <- seq_len(nrow(markerList_PharyngealMesoderm_up_df_homer))
markerList_PharyngealMesoderm_up_df_homer_trim <- markerList_PharyngealMesoderm_up_df_homer[,c(1,2,3,9,5)]
write.table(markerList_PharyngealMesoderm_up_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/markerList_PharyngealMesoderm_up_df_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)

## go to the folder DAR_homer and run script 5f

### Part 2: Export Peaks Enriched in WT v KO (ie, up-regulated in KO) for clusters of interest
## run this in a new session after running 5a, selecting for peaks <= -1

markerList_AHF_df_homer <- markerList_AHF_df[order(markerList_AHF_df$Log2FC, decreasing = F),]
markerList_AHF_df_homer$rank <- seq_len(nrow(markerList_AHF_df_homer))
markerList_AHF_df_homer_trim <- markerList_AHF_df_homer[,c(1,2,3,9,5)]
write.table(markerList_AHF_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/upKO/markerList_AHF_df_upKO_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)

markerList_PharyngealMesoderm_df_homer <- markerList_PharyngealMesoderm_df[order(markerList_PharyngealMesoderm_df$Log2FC, decreasing = F),]
markerList_PharyngealMesoderm_df_homer$rank <- seq_len(nrow(markerList_PharyngealMesoderm_df_homer))
markerList_PharyngealMesoderm_df_homer_trim <- markerList_PharyngealMesoderm_df_homer[,c(1,2,3,9,5)]
write.table(markerList_PharyngealMesoderm_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/upKO/markerList_PharyngealMesoderm_df_upKO_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)

markerList_Endoderm_df_homer <- markerList_Endoderm_df[order(markerList_Endoderm_df$Log2FC, decreasing = F),]
markerList_Endoderm_df_homer$rank <- seq_len(nrow(markerList_Endoderm_df_homer))
markerList_Endoderm_df_homer_trim <- markerList_Endoderm_df_homer[,c(1,2,3,9,5)]
write.table(markerList_Endoderm_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/upKO/markerList_Endoderm_df_upKO_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)

markerList_Ectoderm_df_homer <- markerList_Ectoderm_df[order(markerList_Ectoderm_df$Log2FC, decreasing = F),]
markerList_Ectoderm_df_homer$rank <- seq_len(nrow(markerList_Ectoderm_df_homer))
markerList_Ectoderm_df_homer_trim <- markerList_Ectoderm_df_homer[,c(1,2,3,9,5)]
write.table(markerList_Ectoderm_df_homer_trim, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/upKO/markerList_Ectoderm_df_upKO_homer_trim.txt",
            quote = F, sep = '\t',col.names = F, row.names = F)

## go to the folder DAR_homer and run script 5f

