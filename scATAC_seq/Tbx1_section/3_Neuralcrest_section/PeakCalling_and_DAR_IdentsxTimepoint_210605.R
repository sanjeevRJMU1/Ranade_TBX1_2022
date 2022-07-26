############################ 
## Peak calling on named Idents x Genotype x Time Point = Clusters7
############################ 
table(mesoderm_subset$Clusters7)

## 9. Pseudo-bulk Replicates in ArchR
mesoderm_subset <- addGroupCoverages(ArchRProj = mesoderm_subset, groupBy = "Clusters7", minCells = 40, maxCells = 1500, force = T)
#2021-06-05 13:05:42 : Finished Creation of Coverage Files!, 66.708 mins elapsed.

## 10. Peak calling
mesoderm_subset <- addReproduciblePeakSet(
  ArchRProj = mesoderm_subset, 
  groupBy = "Clusters7", 
  pathToMacs2 = "//anaconda3/bin/macs2",
  peaksPerCell = 1500,
  maxPeaks = 150000,
  minCells = 25,
  force = T
)

# Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
# Anterior_Second_Heart_Field_05_K   Anterior_Second_Heart_Field_05_K    200        200           2   98  102   150000
# Anterior_Second_Heart_Field_05_W   Anterior_Second_Heart_Field_05_W    296        296           2  140  156   150000
# Anterior_Second_Heart_Field_15_K   Anterior_Second_Heart_Field_15_K    437        437           3  101  208   150000
# Anterior_Second_Heart_Field_15_W   Anterior_Second_Heart_Field_15_W    597        597           2  258  339   150000
# Anterior_Second_Heart_Field_25_K   Anterior_Second_Heart_Field_25_K    286        286           2  106  180   150000
# Anterior_Second_Heart_Field_25_W   Anterior_Second_Heart_Field_25_W    389        389           2  134  255   150000
# Cardiomyocyte_05_K                               Cardiomyocyte_05_K    255        255           2   95  160   150000
# Cardiomyocyte_05_W                               Cardiomyocyte_05_W    163        163           2   79   84   150000
# Cardiomyocyte_15_K                               Cardiomyocyte_15_K    778        778           3  133  380   150000
# Cardiomyocyte_15_W                               Cardiomyocyte_15_W    358        358           2  151  207   150000
# Cardiomyocyte_25_K                               Cardiomyocyte_25_K    658        658           2  313  345   150000
# Cardiomyocyte_25_W                               Cardiomyocyte_25_W    788        788           2  241  547   150000
# Cardiopharyngeal_Mesenchyme_05_K   Cardiopharyngeal_Mesenchyme_05_K    582        582           2  290  292   150000
# Cardiopharyngeal_Mesenchyme_05_W   Cardiopharyngeal_Mesenchyme_05_W    678        678           2  302  376   150000
# Cardiopharyngeal_Mesenchyme_15_K   Cardiopharyngeal_Mesenchyme_15_K   2284       2284           3  365 1186   150000
# Cardiopharyngeal_Mesenchyme_15_W   Cardiopharyngeal_Mesenchyme_15_W   1243       1243           2  449  794   150000
# Cardiopharyngeal_Mesenchyme_25_K   Cardiopharyngeal_Mesenchyme_25_K    563        563           2  202  361   150000
# Cardiopharyngeal_Mesenchyme_25_W   Cardiopharyngeal_Mesenchyme_25_W    608        608           2  303  305   150000
# Cardiopharyngeal_Mesoderm_05_K       Cardiopharyngeal_Mesoderm_05_K    872        872           2  346  526   150000
# Cardiopharyngeal_Mesoderm_05_W       Cardiopharyngeal_Mesoderm_05_W    741        741           2  315  426   150000
# Cardiopharyngeal_Mesoderm_15_K       Cardiopharyngeal_Mesoderm_15_K   1284       1284           3  290  628   150000
# Cardiopharyngeal_Mesoderm_15_W       Cardiopharyngeal_Mesoderm_15_W    788        788           2  240  548   150000
# Cardiopharyngeal_Mesoderm_25_K       Cardiopharyngeal_Mesoderm_25_K    565        565           2  237  328   150000
# Cardiopharyngeal_Mesoderm_25_W       Cardiopharyngeal_Mesoderm_25_W    565        565           2  215  350   150000
# Cranial_Mesenchyme_05_K                     Cranial_Mesenchyme_05_K    500        500           2  155  345   150000
# Cranial_Mesenchyme_05_W                     Cranial_Mesenchyme_05_W    511        511           2  242  269   150000
# Cranial_Mesenchyme_15_K                     Cranial_Mesenchyme_15_K   1525       1525           3  312  749   150000
# Cranial_Mesenchyme_15_W                     Cranial_Mesenchyme_15_W   1057       1057           2  399  658   150000
# Cranial_Mesenchyme_25_K                     Cranial_Mesenchyme_25_K    289        289           2   95  194   150000
# Cranial_Mesenchyme_25_W                     Cranial_Mesenchyme_25_W    356        356           2  138  218   150000
# Epicardium_05_K                                     Epicardium_05_K    779        779           2  143  636   150000
# Epicardium_05_W                                     Epicardium_05_W    565        565           2  226  339   150000
# Epicardium_15_K                                     Epicardium_15_K   1028       1028           3  145  568   150000
# Epicardium_15_W                                     Epicardium_15_W    212        212           2   92  120   150000
# Epicardium_25_K                                     Epicardium_25_K    551        551           2  182  369   150000
# Epicardium_25_W                                     Epicardium_25_W    460        460           2  221  239   150000
# Juxtacardiac_Field_05_K                     Juxtacardiac_Field_05_K    568        568           2  235  333   150000
# Juxtacardiac_Field_05_W                     Juxtacardiac_Field_05_W    735        735           2  338  397   150000
# Juxtacardiac_Field_15_K                     Juxtacardiac_Field_15_K    495        495           3  122  207   150000
# Juxtacardiac_Field_15_W                     Juxtacardiac_Field_15_W   1174       1174           2  470  704   150000
# Juxtacardiac_Field_25_K                     Juxtacardiac_Field_25_K    110        110           2   47   63   150000
# Juxtacardiac_Field_25_W                     Juxtacardiac_Field_25_W     90         90           2   40   50   135000
# Neural_Crest_05_K                                 Neural_Crest_05_K   2071       2071           2  875 1196   150000
# Neural_Crest_05_W                                 Neural_Crest_05_W   2213       2213           2 1029 1184   150000
# Neural_Crest_15_K                                 Neural_Crest_15_K    998        998           3  177  471   150000
# Neural_Crest_15_W                                 Neural_Crest_15_W   3131       2826           2 1326 1500   150000
# Neural_Crest_25_K                                 Neural_Crest_25_K    370        370           2  131  239   150000
# Neural_Crest_25_W                                 Neural_Crest_25_W    308        308           2   54  254   150000
# Paraxial_Mesoderm_05_K                       Paraxial_Mesoderm_05_K    515        515           2  181  334   150000
# Paraxial_Mesoderm_05_W                       Paraxial_Mesoderm_05_W    604        604           2  280  324   150000
# Paraxial_Mesoderm_15_K                       Paraxial_Mesoderm_15_K   2358       2358           3  447 1259   150000
# Paraxial_Mesoderm_15_W                       Paraxial_Mesoderm_15_W    647        647           2  262  385   150000
# Paraxial_Mesoderm_25_K                       Paraxial_Mesoderm_25_K    382        382           2  169  213   150000
# Paraxial_Mesoderm_25_W                       Paraxial_Mesoderm_25_W    556        556           2  254  302   150000
# Posterior_Second_Heart_Field_05_K Posterior_Second_Heart_Field_05_K    935        935           2  333  602   150000
# Posterior_Second_Heart_Field_05_W Posterior_Second_Heart_Field_05_W    550        550           2  250  300   150000
# Posterior_Second_Heart_Field_15_K Posterior_Second_Heart_Field_15_K   3163       3124           3  702 1500   150000
# Posterior_Second_Heart_Field_15_W Posterior_Second_Heart_Field_15_W   1667       1667           2  703  964   150000
# Posterior_Second_Heart_Field_25_K Posterior_Second_Heart_Field_25_K    696        696           2  271  425   150000
# Posterior_Second_Heart_Field_25_W Posterior_Second_Heart_Field_25_W    605        605           2  291  314   150000
# Smooth_Muscle_05_K                               Smooth_Muscle_05_K     62         59           2   40   40    88500
# Smooth_Muscle_05_W                               Smooth_Muscle_05_W     64         60           2   40   40    90000
# Smooth_Muscle_15_K                               Smooth_Muscle_15_K    463        463           3  125  192   150000
# Smooth_Muscle_15_W                               Smooth_Muscle_15_W    276        276           2   94  182   150000
# Smooth_Muscle_25_K                               Smooth_Muscle_25_K     49         44           2   25   31    66000
# Smooth_Muscle_25_W                               Smooth_Muscle_25_W     81         81           2   40   41   121500
#2021-06-05 15:28:37 : Finished Creating Union Peak Set (491076)!, 142.896 mins elapsed.

# look at the peak set and then save it as a separate object
getPeakSet(mesoderm_subset)
merged_peak_setGR <- getPeakSet(mesoderm_subset)

# add peakmatrix to object
mesoderm_subset <- addPeakMatrix(mesoderm_subset)
getAvailableMatrices(mesoderm_subset)

saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517", load = TRUE)


############################ ############################ ############################ ############################ 
# DAR
#
markerTest_Anterior_Second_Heart_Field_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Anterior_Second_Heart_Field_25_W",
  bgdGroups = "Anterior_Second_Heart_Field_25_K"
)
pma_Anterior_Second_Heart_Field_E925 <- plotMarkers(seMarker = markerTest_Anterior_Second_Heart_Field_E925, name = "Anterior_Second_Heart_Field_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Anterior_Second_Heart_Field_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Anterior_Second_Heart_Field_05_W",
  bgdGroups = "Anterior_Second_Heart_Field_05_K"
)
pma_Anterior_Second_Heart_Field_E105 <- plotMarkers(seMarker = markerTest_Anterior_Second_Heart_Field_E105, name = "Anterior_Second_Heart_Field_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Anterior_Second_Heart_Field_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Anterior_Second_Heart_Field_15_W",
  bgdGroups = "Anterior_Second_Heart_Field_15_K"
)
pma_Anterior_Second_Heart_Field_E115 <- plotMarkers(seMarker = markerTest_Anterior_Second_Heart_Field_E115, name = "Anterior_Second_Heart_Field_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Cardiomyocyte_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiomyocyte_25_W",
  bgdGroups = "Cardiomyocyte_25_K"
)
pma_Cardiomyocyte_E925 <- plotMarkers(seMarker = markerTest_Cardiomyocyte_E925, name = "Cardiomyocyte_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cardiomyocyte_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiomyocyte_05_W",
  bgdGroups = "Cardiomyocyte_05_K"
)
pma_Cardiomyocyte_E105 <- plotMarkers(seMarker = markerTest_Cardiomyocyte_E105, name = "Cardiomyocyte_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cardiomyocyte_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiomyocyte_15_W",
  bgdGroups = "Cardiomyocyte_15_K",
  maxCells = 358
)
pma_Cardiomyocyte_E115 <- plotMarkers(seMarker = markerTest_Cardiomyocyte_E115, name = "Cardiomyocyte_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Cardiopharyngeal_Mesenchyme_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiopharyngeal_Mesenchyme_25_W",
  bgdGroups = "Cardiopharyngeal_Mesenchyme_25_K"
)
pma_Cardiopharyngeal_Mesenchyme_E925 <- plotMarkers(seMarker = markerTest_Cardiopharyngeal_Mesenchyme_E925, name = "Cardiopharyngeal_Mesenchyme_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cardiopharyngeal_Mesenchyme_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiopharyngeal_Mesenchyme_05_W",
  bgdGroups = "Cardiopharyngeal_Mesenchyme_05_K"
)
pma_Cardiopharyngeal_Mesenchyme_E105 <- plotMarkers(seMarker = markerTest_Cardiopharyngeal_Mesenchyme_E105, name = "Cardiopharyngeal_Mesenchyme_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cardiopharyngeal_Mesenchyme_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiopharyngeal_Mesenchyme_15_W",
  bgdGroups = "Cardiopharyngeal_Mesenchyme_15_K"
)
pma_Cardiopharyngeal_Mesenchyme_E115 <- plotMarkers(seMarker = markerTest_Cardiopharyngeal_Mesenchyme_E115, name = "Cardiopharyngeal_Mesenchyme_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Cardiopharyngeal_Mesoderm_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiopharyngeal_Mesoderm_25_W",
  bgdGroups = "Cardiopharyngeal_Mesoderm_25_K"
)
pma_Cardiopharyngeal_Mesoderm_E925 <- plotMarkers(seMarker = markerTest_Cardiopharyngeal_Mesoderm_E925, name = "Cardiopharyngeal_Mesoderm_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cardiopharyngeal_Mesoderm_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiopharyngeal_Mesoderm_05_W",
  bgdGroups = "Cardiopharyngeal_Mesoderm_05_K"
)
pma_Cardiopharyngeal_Mesoderm_E105 <- plotMarkers(seMarker = markerTest_Cardiopharyngeal_Mesoderm_E105, name = "Cardiopharyngeal_Mesoderm_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cardiopharyngeal_Mesoderm_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cardiopharyngeal_Mesoderm_15_W",
  bgdGroups = "Cardiopharyngeal_Mesoderm_15_K"
)
pma_Cardiopharyngeal_Mesoderm_E115 <- plotMarkers(seMarker = markerTest_Cardiopharyngeal_Mesoderm_E115, name = "Cardiopharyngeal_Mesoderm_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Cranial_Mesenchyme_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cranial_Mesenchyme_25_W",
  bgdGroups = "Cranial_Mesenchyme_25_K"
)
pma_Cranial_Mesenchyme_E925 <- plotMarkers(seMarker = markerTest_Cranial_Mesenchyme_E925, name = "Cranial_Mesenchyme_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cranial_Mesenchyme_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cranial_Mesenchyme_05_W",
  bgdGroups = "Cranial_Mesenchyme_05_K"
)
pma_Cranial_Mesenchyme_E105 <- plotMarkers(seMarker = markerTest_Cranial_Mesenchyme_E105, name = "Cranial_Mesenchyme_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Cranial_Mesenchyme_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Cranial_Mesenchyme_15_W",
  bgdGroups = "Cranial_Mesenchyme_15_K"
)
pma_Cranial_Mesenchyme_E115 <- plotMarkers(seMarker = markerTest_Cranial_Mesenchyme_E115, name = "Cranial_Mesenchyme_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Epicardium_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Epicardium_25_W",
  bgdGroups = "Epicardium_25_K"
)
pma_Epicardium_E925 <- plotMarkers(seMarker = markerTest_Epicardium_E925, name = "Epicardium_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Epicardium_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Epicardium_05_W",
  bgdGroups = "Epicardium_05_K"
)
pma_Epicardium_E105 <- plotMarkers(seMarker = markerTest_Epicardium_E105, name = "Epicardium_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Epicardium_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Epicardium_15_W",
  bgdGroups = "Epicardium_15_K"
)
pma_Epicardium_E115 <- plotMarkers(seMarker = markerTest_Epicardium_E115, name = "Epicardium_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Juxtacardiac_Field_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Juxtacardiac_Field_25_W",
  bgdGroups = "Juxtacardiac_Field_25_K"
)
pma_Juxtacardiac_Field_E925 <- plotMarkers(seMarker = markerTest_Juxtacardiac_Field_E925, name = "Juxtacardiac_Field_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Juxtacardiac_Field_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Juxtacardiac_Field_05_W",
  bgdGroups = "Juxtacardiac_Field_05_K"
)
pma_Juxtacardiac_Field_E105 <- plotMarkers(seMarker = markerTest_Juxtacardiac_Field_E105, name = "Juxtacardiac_Field_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Juxtacardiac_Field_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Juxtacardiac_Field_15_W",
  bgdGroups = "Juxtacardiac_Field_15_K"
)
pma_Juxtacardiac_Field_E115 <- plotMarkers(seMarker = markerTest_Juxtacardiac_Field_E115, name = "Juxtacardiac_Field_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Neural_Crest_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Neural_Crest_25_W",
  bgdGroups = "Neural_Crest_25_K"
)
pma_Neural_Crest_E925 <- plotMarkers(seMarker = markerTest_Neural_Crest_E925, name = "Neural_Crest_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Neural_Crest_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Neural_Crest_05_W",
  bgdGroups = "Neural_Crest_05_K"
)
pma_Neural_Crest_E105 <- plotMarkers(seMarker = markerTest_Neural_Crest_E105, name = "Neural_Crest_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Neural_Crest_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Neural_Crest_15_W",
  bgdGroups = "Neural_Crest_15_K"
)
pma_Neural_Crest_E115 <- plotMarkers(seMarker = markerTest_Neural_Crest_E115, name = "Neural_Crest_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Paraxial_Mesoderm_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Paraxial_Mesoderm_25_W",
  bgdGroups = "Paraxial_Mesoderm_25_K"
)
pma_Paraxial_Mesoderm_E925 <- plotMarkers(seMarker = markerTest_Paraxial_Mesoderm_E925, name = "Paraxial_Mesoderm_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Paraxial_Mesoderm_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Paraxial_Mesoderm_05_W",
  bgdGroups = "Paraxial_Mesoderm_05_K"
)
pma_Paraxial_Mesoderm_E105 <- plotMarkers(seMarker = markerTest_Paraxial_Mesoderm_E105, name = "Paraxial_Mesoderm_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Paraxial_Mesoderm_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Paraxial_Mesoderm_15_W",
  bgdGroups = "Paraxial_Mesoderm_15_K"
)
pma_Paraxial_Mesoderm_E115 <- plotMarkers(seMarker = markerTest_Paraxial_Mesoderm_E115, name = "Paraxial_Mesoderm_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Posterior_Second_Heart_Field_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Posterior_Second_Heart_Field_25_W",
  bgdGroups = "Posterior_Second_Heart_Field_25_K"
)
pma_Posterior_Second_Heart_Field_E925 <- plotMarkers(seMarker = markerTest_Posterior_Second_Heart_Field_E925, name = "Posterior_Second_Heart_Field_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Posterior_Second_Heart_Field_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Posterior_Second_Heart_Field_05_W",
  bgdGroups = "Posterior_Second_Heart_Field_05_K"
)
pma_Posterior_Second_Heart_Field_E105 <- plotMarkers(seMarker = markerTest_Posterior_Second_Heart_Field_E105, name = "Posterior_Second_Heart_Field_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Posterior_Second_Heart_Field_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Posterior_Second_Heart_Field_15_W",
  bgdGroups = "Posterior_Second_Heart_Field_15_K"
)
pma_Posterior_Second_Heart_Field_E115 <- plotMarkers(seMarker = markerTest_Posterior_Second_Heart_Field_E115, name = "Posterior_Second_Heart_Field_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

#
markerTest_Smooth_Muscle_E925 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Smooth_Muscle_25_W",
  bgdGroups = "Smooth_Muscle_25_K"
)
pma_Smooth_Muscle_E925 <- plotMarkers(seMarker = markerTest_Smooth_Muscle_E925, name = "Smooth_Muscle_25_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Smooth_Muscle_E105 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Smooth_Muscle_05_W",
  bgdGroups = "Smooth_Muscle_05_K"
)
pma_Smooth_Muscle_E105 <- plotMarkers(seMarker = markerTest_Smooth_Muscle_E105, name = "Smooth_Muscle_05_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")

markerTest_Smooth_Muscle_E115 <- getMarkerFeatures(
  ArchRProj = mesoderm_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters7",
  binarize = F,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Smooth_Muscle_15_W",
  bgdGroups = "Smooth_Muscle_15_K"
)
pma_Smooth_Muscle_E115 <- plotMarkers(seMarker = markerTest_Smooth_Muscle_E115, name = "Smooth_Muscle_15_W", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")


