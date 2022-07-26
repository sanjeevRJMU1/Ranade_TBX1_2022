###### Tbx1 WT v KO all aggr
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
## make mesoderm + CNC subset for now = subset1

#plots
plotEmbedding(ArchRProj = Tbx1_all_aggr, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
table(Tbx1_all_aggr$Clusters,Tbx1_all_aggr$Sample)

# Subset clusters of interest in the hackiest fucking way 
# clone object
mesoderm_subset1 <- Tbx1_all_aggr

# assign Y/N values to only the clusters you want
library(plyr)
temp_cluster_names <- as.character(revalue(mesoderm_subset1$Clusters, c(
  "C1"="False",
  "C2"="False",
  "C3"="False",
  "C4"="True",
  "C5"="True",
  "C6"="True",
  "C7"="True",
  "C8"="True",
  "C9"="True",
  "C10"="True",
  "C11"="True",
  "C12"="True",
  "C13"="True",
  "C14"="True",
  "C15"="False",
  "C16"="False",
  "C17"="True",
  "C18"="False",
  "C19"="False",
  "C20"="False",
  "C21"="False",
  "C22"="False",
  "C23"="False",
  "C24"="False",
  "C25"="False",
  "C26"="False",
  "C27"="False",
  "C28"="False"
)))
length(unique(temp_cluster_names))
mesoderm_subset1$Clusters2 <- temp_cluster_names
mesoderm_subset1@cellColData@listData[["Clusters2"]] <- temp_cluster_names
table(mesoderm_subset1$Clusters2,mesoderm_subset1$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 j_E115_WT_1
# False        3531        2220        2150        3216        1934        3215        2868        2296        1328        1119
# True         2743        1778        2068        2744        3626        3724        3806        3329        4528        4616
# 
# k_E115_WT_2 l_E115_KO_2 m_E115_KO_3
# False        1578         759        1715
# True         6704        2929        7538

## subset
idxSample <- BiocGenerics::which(mesoderm_subset1$Clusters2 %in% "True")
cellsSample <- mesoderm_subset1$cellNames[idxSample]
subsetArchRProject(
  ArchRProj = mesoderm_subset1,
  cells = cellsSample,
  outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset1_210517", force = T
)
## Dropping ImputeWeights Since You Are Subsetting Cells! ImputeWeights is a cell-x-cell Matrix!













