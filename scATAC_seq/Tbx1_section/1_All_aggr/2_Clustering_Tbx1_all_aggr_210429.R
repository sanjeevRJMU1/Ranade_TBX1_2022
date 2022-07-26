###### Tbx1 WT v KO all aggr
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
##### Chapters 4 - 6, Clustering with batch correction

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1


######## 
# Re-load ArchR project!
Tbx1_all_aggr <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_all_aggr")
Tbx1_all_aggr
######## 

### 4. Iterative LSI (TF-IDF and SVD)
Tbx1_all_aggr <- addIterativeLSI(
  ArchRProj = Tbx1_all_aggr,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 3, 
  clusterParams = list(
    resolution = c(0.1,0.2),
    sampleCells = 10000,
    maxClusters = 40,
    n.start = 10), 
  varFeatures = 25000, 
  sampleCellsPre = 10000,
  dimsToUse = 2:30,
  LSIMethod = 2,
  seed = 1,
  force = T
)

# # # # Batch correction with Harmony --> from now on, use "Harmony" for reducedDims and "UMAP" for embedding
Tbx1_all_aggr <- addHarmony(
  ArchRProj = Tbx1_all_aggr,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  dimsToUse = 2:30,
  groupBy = "Sample", force = T
)

#### 5. Clustering
Tbx1_all_aggr <- addClusters(
  input = Tbx1_all_aggr,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  maxClusters = 40,
  dimsToUse = 2:30,
  resolution = 0.6, force = T
)

# ## Confusion matrix of sample number vs cluster 
cM <- confusionMatrix(paste0(Tbx1_all_aggr$Clusters), paste0(Tbx1_all_aggr$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p
plotPDF(p, name = "Plot-ConfusionMatrix-Heatmap.pdf", ArchRProj = Tbx1_all_aggr, addDOC = FALSE, width = 5, height = 5)

#### 6. Single-cell Embeddings
### UMAP with Harmony
Tbx1_all_aggr <- addUMAP(
  ArchRProj = Tbx1_all_aggr, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony",
  metric = "cosine",
  force = T
)

p9 <- plotEmbedding(ArchRProj = Tbx1_all_aggr, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p10 <- plotEmbedding(ArchRProj = Tbx1_all_aggr, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p9, p10, type = "h")
plotPDF(p9,p10, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = Tbx1_all_aggr, addDOC = FALSE, width = 5, height = 5)

table(Tbx1_all_aggr$Clusters,Tbx1_all_aggr$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 j_E115_WT_1
# C1          226         185         151         210         197         153         204         149         281         188
# C10         368         286         352         492         773         330         357         239         911         533
# C11          42          31          48          34          61          48          69          49          91         376
# C12         254         128          58         241         729        1070        1013         849         237         994
# C13         218         159         252         207         169         171         207         198         328         311
# C14         495         308         382         437         438         724         592         489         741         430
# C15          26          13           7          12           4           7           9           5           1          18
# C16           2           6           4           7          13          56          52          21           5          46
# C17         159          81          69         143          89          89         171         164         129         288
# C18         404         284         211         387         153         263         270          92         170          88
# C19         230          92         122         143          94         104          78          53          25          20
# C2           14          17           9          16           8          21          25           8          10         166
# C20         466         220         318         375         372         186         344         226         212         213
# C21         588         374         315         550         332         341         360         300         129         158
# C22           2           2           3           3         110         206         138          72           2           5
# C23         181         124          81         182          99         363         233         231          48          31
# C24          99         100         200         151         223         288         315         339         301           4
# C25           0           0           0           0           0           0           0           0           0         116
# C26         111          68          75         101          61          90          71          64          16           7
# C27        1170         721         643        1069         231        1088         678         684          46          25
# C28           6           8           3           2           2           1           2           1          26           0
# C3            6           6           8           8          35          48          89          51          56          34
# C4          577         344         278         393         155          79          64          69         297          87
# C5            6           4           3           0          27          28          24          23          19          45
# C6          231         153         221         275         228         248         295         222         608         346
# C7          103          87         134         162         109         282         175         197         457         302
# C8          126          80         111         147         372         493         509         570         377         737
# C9          164         117         160         213         476         162         330         260         333         167
# 
# k_E115_WT_2 l_E115_KO_2 m_E115_KO_3
# C1          230         177         445
# C10         974         717        1547
# C11         138          43         142
# C12        1437          92         377
# C13         493         178         568
# C14         632         567        1359
# C15           3           4           3
# C16           8           7          10
# C17         413         105         212
# C18         189         103         218
# C19          37          10          38
# C2           24          11          20
# C20         400         101         221
# C21         291          52         110
# C22           1           0           2
# C23          56          18          49
# C24         180         224         341
# C25           0           0           0
# C26          35           6          28
# C27          48          24          28
# C28          25           1          44
# C3           51          21         158
# C4          146         163         420
# C5           71           9          23
# C6          580         312        1046
# C7          496         320         768
# C8         1126         265         468
# C9          198         158         608

saveArchRProject(ArchRProj = Tbx1_all_aggr, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_all_aggr", load = TRUE)


####################################################
## Proceed to 3. scRNA-seq integration and cluster ID



