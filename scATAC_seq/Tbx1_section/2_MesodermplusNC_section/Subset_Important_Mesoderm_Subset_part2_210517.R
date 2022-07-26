#########
# Having run thru clustering + peak calling on Clusters1 + homer annotate, here are the conclusions:
# 1. clusters 11 and 17 is entirely due to a dissection artifact. kick it out.
# 2. i'm over-clustered, clusters 7 and 17 are not real

## kick out cluster 11, 17. for now, leave 2 and 7 and see if you can lower the clustering resolution
### importantly, delete group coverages and peak calls!
### save to new directory and delete old directory all together

### update post clustering on subset2. kick out clusters 2 and 11 and 14

# assign Y/N values to only the clusters you want
library(plyr)
temp_cluster_names <- as.character(revalue(mesoderm_subset$Clusters, c(
  "C1"="True",
  "C2"="False",
  "C3"="True",
  "C4"="True",
  "C5"="True",
  "C6"="True",
  "C7"="True",
  "C8"="True",
  "C9"="True",
  "C10"="True",
  "C11"="False",
  "C12"="True",
  "C13"="True",
  "C14"="False",
  "C15"="True",
  "C16"="True"
)))
length(unique(temp_cluster_names))
mesoderm_subset$Clusters4 <- temp_cluster_names
mesoderm_subset@cellColData@listData[["Clusters4"]] <- temp_cluster_names
table(mesoderm_subset$Clusters4,mesoderm_subset$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 j_E115_WT_1
# False          45          29          42          34          65          62          72          54          75         249
# True         2696        1748        2023        2708        3557        3655        3727        3267        4432        4223
# 
# k_E115_WT_2 l_E115_KO_2 m_E115_KO_3
# False         165          42         112
# True         6513        2877        7275

## subset
idxSample <- BiocGenerics::which(mesoderm_subset1$Clusters3 %in% "True")
cellsSample <- mesoderm_subset$cellNames[idxSample]
subsetArchRProject(
  ArchRProj = mesoderm_subset,
  cells = cellsSample,
  outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517", force = T
)



