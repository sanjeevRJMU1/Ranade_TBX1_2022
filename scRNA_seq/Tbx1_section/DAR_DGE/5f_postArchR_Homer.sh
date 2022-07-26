### script for taking DAR from ArchR script 5 --> Homer
## 1. open the txt file in excel and simply add a blank column in column 5
## 2. save in excel as a "txt" but change the file name, add .bed --> technically this is still a txt file, just with a name edit

findMotifsGenome.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/markerList_AHF_df_homer_trim.bed.txt mm10 /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/AHF_motif_output/ -size given -mask
findMotifsGenome.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/markerList_Ectoderm_df_homer_trim.bed.txt mm10 /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/Ectoderm_motif_output/ -size given -mask
findMotifsGenome.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/markerList_Endoderm_df_homer_trim.bed.txt mm10 /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/Endoderm_motif_output/ -size given -mask
findMotifsGenome.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/markerList_PharyngealMesoderm_df_homer_trim.bed.txt mm10 /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/PharyngealMesoderm_motif_output/ -size given -mask


annotatePeaks.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/markerList_AHF_df_homer_trim.bed.txt mm10 > /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/annotatemotifs/AHF_annotate.txt
annotatePeaks.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/markerList_Ectoderm_df_homer_trim.bed.txt mm10 > /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/annotatemotifs/Ectoderm_annotate.txt
annotatePeaks.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/markerList_Endoderm_df_homer_trim.bed.txt mm10 > /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/annotatemotifs/Endoderm_annotate.txt
annotatePeaks.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/markerList_PharyngealMesoderm_df_homer_trim.bed.txt mm10 > /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/annotatemotifs/PharyngealMesoderm_annotate.txt



##### on up-regulated peaks
findMotifsGenome.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/upKO/markerList_AHF_df_upKO_homer_trim.txt mm10 /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/DAR_markerList/upKO/AHF_upKO_motif_output/ -size given -mask
findMotifsGenome.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/markerList_Ectoderm_df_homer_trim.bed.txt mm10 /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/Ectoderm_motif_output/ -size given -mask
findMotifsGenome.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/markerList_Endoderm_df_homer_trim.bed.txt mm10 /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/Endoderm_motif_output/ -size given -mask
findMotifsGenome.pl /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/markerList_PharyngealMesoderm_df_homer_trim.bed.txt mm10 /Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/PharyngealMesoderm_motif_output/ -size given -mask
