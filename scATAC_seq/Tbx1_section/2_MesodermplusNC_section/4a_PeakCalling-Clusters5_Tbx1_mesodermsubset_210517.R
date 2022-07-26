###### Tbx1 WT v KO mesoderm subset1
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
#### Chapters 9 and 10 Pseudo-bulk Replicates and Peak calling
#### This version = groupcoverage on Clusters5 = Clusters1 divided into WT and KO per time point

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


# Re-load ArchR project!
mesoderm_subset <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517")
mesoderm_subset
######## 


## 9. Pseudo-bulk Replicates in ArchR
mesoderm_subset <- addGroupCoverages(ArchRProj = mesoderm_subset, groupBy = "Clusters5", minCells = 40, maxCells = 1500, force = T)
#2021-05-17 23:22:50 : Finished Creation of Coverage Files!, 67.955 mins elapsed.

## 10. Peak calling
mesoderm_subset <- addReproduciblePeakSet(
  ArchRProj = mesoderm_subset, 
  groupBy = "Clusters5", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 1500,
  maxPeaks = 150000,
  minCells = 25,
  force = T
)

# Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
# C1_05_K   C1_05_K    201        201           2   67  134   150000
# C1_05_W   C1_05_W    115        115           2   56   59   150000
# C1_15_K   C1_15_K    727        727           3  125  355   150000
# C1_15_W   C1_15_W    226        226           2   98  128   150000
# C1_25_K   C1_25_K    657        657           2  312  345   150000
# C1_25_W   C1_25_W    783        783           2  240  543   150000
# C2_05_K   C2_05_K     54         54           2   40   40    81000
# C2_05_W   C2_05_W     48         44           2   27   29    66000
# C2_15_K   C2_15_K     51         51           2   40   40    76500
# C2_15_W   C2_15_W    132        132           2   53   79   150000
# C2_25_K   C2_25_K      1          1           2    1    1     1500
# C2_25_W   C2_25_W      5          5           2    5    5     7500
# C3_05_K   C3_05_K    509        509           2  253  256   150000
# C3_05_W   C3_05_W    601        601           2  271  330   150000
# C3_15_K   C3_15_K   2121       2121           3  335 1107   150000
# C3_15_W   C3_15_W   1112       1112           2  392  720   150000
# C3_25_K   C3_25_K    506        506           2  177  329   150000
# C3_25_W   C3_25_W    527        527           2  262  265   150000
# C4_05_K   C4_05_K    500        500           2  155  345   150000
# C4_05_W   C4_05_W    511        511           2  242  269   150000
# C4_15_K   C4_15_K   1525       1525           3  312  749   150000
# C4_15_W   C4_15_W   1057       1057           2  399  658   150000
# C4_25_K   C4_25_K    289        289           2   95  194   150000
# C4_25_W   C4_25_W    356        356           2  138  218   150000
# C5_05_K   C5_05_K     63         59           2   40   40    88500
# C5_05_W   C5_05_W    124        124           2   61   63   150000
# C5_15_K   C5_15_K    356        356           3   67  192   150000
# C5_15_W   C5_15_W    584        584           2  269  315   150000
# C5_25_K   C5_25_K      4          4           2    4    4     6000
# C6_05_K   C6_05_K    426        426           2  115  311   150000
# C6_05_W   C6_05_W    364        364           2  125  239   150000
# C6_15_K   C6_15_K    161        136           2   67   69   150000
# C6_15_W   C6_15_W    462        462           2  152  310   150000
# C6_25_K   C6_25_K    197        197           2   72  125   150000
# C6_25_W   C6_25_W    138        138           2   40   98   150000
# C7_05_K   C7_05_K   1582       1582           2  709  873   150000
# C7_05_W   C7_05_W   1725       1725           2  841  884   150000
# C7_15_K   C7_15_K    481        481           3   85  210   150000
# C7_15_W   C7_15_W   2085       2085           2  905 1180   150000
# C7_25_K   C7_25_K    169        169           2   56  113   150000
# C7_25_W   C7_25_W    170        170           2   40  130   150000
# C8_05_K   C8_05_K    779        779           2  143  636   150000
# C8_05_W   C8_05_W    565        565           2  226  339   150000
# C8_15_K   C8_15_K   1028       1028           3  145  568   150000
# C8_15_W   C8_15_W    212        212           2   92  120   150000
# C8_25_K   C8_25_K    551        551           2  182  369   150000
# C8_25_W   C8_25_W    460        460           2  221  239   150000
# C9_05_K   C9_05_K     62         59           2   40   40    88500
# C9_05_W   C9_05_W     64         60           2   40   40    90000
# C9_15_K   C9_15_K    463        463           3  125  192   150000
# C9_15_W   C9_15_W    276        276           2   94  182   150000
# C9_25_K   C9_25_K     49         44           2   25   31    66000
# C9_25_W   C9_25_W     81         81           2   40   41   121500
# C10_05_K C10_05_K    568        568           2  235  333   150000
# C10_05_W C10_05_W    735        735           2  338  397   150000
# C10_15_K C10_15_K    495        495           3  122  207   150000
# C10_15_W C10_15_W   1174       1174           2  470  704   150000
# C10_25_K C10_25_K    110        110           2   47   63   150000
# C10_25_W C10_25_W     90         90           2   40   50   135000
# C11_05_K C11_05_K     73         63           2   40   40    94500
# C11_05_W C11_05_W     77         65           2   40   40    97500
# C11_15_K C11_15_K    163        133           2   54   79   150000
# C11_15_W C11_15_W    131        131           2   57   74   150000
# C11_25_K C11_25_K     57         55           2   40   40    82500
# C11_25_W C11_25_W     81         81           2   40   41   121500
# C12_05_K C12_05_K    515        515           2  181  334   150000
# C12_05_W C12_05_W    603        603           2  280  323   150000
# C12_15_K C12_15_K   2343       2343           3  443 1251   150000
# C12_15_W C12_15_W    496        496           2  123  373   150000
# C12_25_K C12_25_K    377        377           2  166  211   150000
# C12_25_W C12_25_W    555        555           2  253  302   150000
# C13_05_K C13_05_K    872        872           2  346  526   150000
# C13_05_W C13_05_W    741        741           2  315  426   150000
# C13_15_K C13_15_K   1284       1284           3  290  628   150000
# C13_15_W C13_15_W    788        788           2  240  548   150000
# C13_25_K C13_25_K    565        565           2  237  328   150000
# C13_25_W C13_25_W    565        565           2  215  350   150000
# C14_05_W C14_05_W      1          1           2    1    1     1500
# C14_15_K C14_15_K     15         15           2   10   13    22500
# C14_15_W C14_15_W    151        151           2   40  111   150000
# C14_25_K C14_25_K      5          5           2    5    5     7500
# C14_25_W C14_25_W      1          1           2    1    1     1500
# C15_05_K C15_05_K    935        935           2  333  602   150000
# C15_05_W C15_05_W    550        550           2  250  300   150000
# C15_15_K C15_15_K   3163       3124           3  702 1500   150000
# C15_15_W C15_15_W   1667       1667           2  703  964   150000
# C15_25_K C15_25_K    696        696           2  271  425   150000
# C15_25_W C15_25_W    605        605           2  291  314   150000
# C16_05_K C16_05_K    200        200           2   98  102   150000
# C16_05_W C16_05_W    296        296           2  140  156   150000
# C16_15_K C16_15_K    437        437           3  101  208   150000
# C16_15_W C16_15_W    597        597           2  258  339   150000
# C16_25_K C16_25_K    286        286           2  106  180   150000
# C16_25_W C16_25_W    389        389           2  134  255   150000
#2021-05-18 01:58:12 : Finished Creating Union Peak Set (530923)!, 155.361 mins elapsed.


# look at the peak set and then save it as a separate object
getPeakSet(mesoderm_subset)
merged_peak_setGR <- getPeakSet(mesoderm_subset)

# add peakmatrix to object
mesoderm_subset <- addPeakMatrix(mesoderm_subset)
getAvailableMatrices(mesoderm_subset)

saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517", load = TRUE)

