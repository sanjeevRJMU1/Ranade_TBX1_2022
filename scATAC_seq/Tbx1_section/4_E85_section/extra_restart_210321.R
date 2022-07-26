###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
##### Re-start project on 03/21


library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 0.9.5


######## 
# Re-load ArchR project!
Tbx1_E825 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825")
Tbx1_E825
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/Tbx1_all_aggr/Tbx1_E825 
# samples(7): c_KO_2 b_WT_1 ... f_WT_3 e_WT_2
# sampleColData names(1): ArrowFiles
# cellColData names(23): Sample TSSEnrichment ... ReadsInPeaks FRIP
# numberOfCells(1): 15178
# medianTSS(1): 11.38
# medianFrags(1): 51252
######## 

## Table of Clusters
# Clusters = simple clustering all cells
table(Tbx1_E825$Clusters)
# C1  C10  C11  C12  C13  C14  C15  C16  C17  C18  C19   C2  C20  C21  C22  C23  C24  C25   C3   C4   C5   C6   C7   C8   C9 
# 739  154  181 1715  277  187  453 1155  482  770  743  322  383  479 1146  639 1225  771  624  761   68  312  136  163 1293 

# Clusters2 = manually annotated clusters
table(Tbx1_E825$Clusters2)
# AHF          CM_Atrium             CM_OFT       CM_Ventricle                CNC           Ectoderm 
# 479                770                482               1155                299               3655 
# Endocardium           Endoderm         Epicardium                LPM PharyngealMesoderm                 PM 
# 639               2978               1126                453                771               1225 
# pSHF 
# 1146 

# Clusters3 = manually annotated clusters, divided by WT v KO
table(Tbx1_E825$Clusters3)
# AHF_K                AHF_W          CM_Atrium_K          CM_Atrium_W             CM_OFT_K 
# 182                  297                  368                  402                  127 
# CM_OFT_W       CM_Ventricle_K       CM_Ventricle_W                CNC_K                CNC_W 
# 355                  397                  758                  259                   40 
# Ectoderm_K           Ectoderm_W        Endocardium_K        Endocardium_W           Endoderm_K 
# 2289                 1366                  300                  339                 1493 
# Endoderm_W         Epicardium_K         Epicardium_W                LPM_K                LPM_W 
# 1485                  514                  612                  216                  237 
# PharyngealMesoderm_K PharyngealMesoderm_W                 PM_K                 PM_W               pSHF_K 
# 496                  275                  630                  595                  642 
# pSHF_W 
# 504


########################################################################################################################
########################################################################################################################
# What I did was re-run the entire workflow from script #2 up to #5 (DAR). I tried different versions of 
# clustering and liked the harmony + dims 2:30. The peak calls right now are for harmony etc for Clusters 4
# overall, the conclusions are that the mesoderm progenitors are fucked, the myocytes are spared and one 
# population of endoderm cells (that have GAS for Tbx1) do show slight defects...
# also, i am really liking the whole aggr object + subset'ing. In the end, that may be the winner here!
# save project and move to hub on 03/22.