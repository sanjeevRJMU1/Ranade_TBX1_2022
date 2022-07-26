
## My manual labeling --> make sure this is a separate metadata column from existing clusters
library(plyr)
temp_cluster_names <- as.character(revalue(wt_atlas$Clusters, c(
  "C1"="Blood",
  "C2"="EndMT",
  "C3"="Endothelium",
  "C4"="Endothelium",
  "C5"="Epicardium",
  "C6"="Epicardium",
  "C7"="Cardiomyocyte",
  "C8"="Cardiomyocyte",
  "C9"="LPM",
  "C10"="Epicardium",
  "C11"="Epicardium",
  "C12"="Ectoderm",
  "C13"="Ectoderm",
  "C14"="Ectoderm",
  "C15"="Ectoderm",
  "C16"="Ectoderm",
  "C17"="Neural_Crest",
  "C18"="Ectoderm",
  "C19"="Endoderm",
  "C20"="Paraxial_Mesoderm",
  "C21"="Neural_Crest",
  "C22"="Neural_Crest",
  "C23"="SHF_Progenitor",
  "C24"="Neural_Crest",
  "C25"="Paraxial_Mesoderm",
  "C26"="Paraxial_Mesoderm",
  "C27"="SHF_Progenitor",
  "C28"="SHF_Progenitor",
  "C29"="SHF_Progenitor",
  "C30"="Neural_Crest",
  "C31"="SHF_Progenitor",
  "C32"="SHF_Progenitor",
  "C33"="Endoderm",
  "C34"="Endoderm",
  "C35"="Endoderm",
  "C36"="Endoderm",
  "C37"="Endoderm",
  "C38"="Endoderm",
  "C39"="Endoderm",
  "C40"="Blood",
  "C41"="Blood"
)))
length(unique(temp_cluster_names))
wt_atlas$Clusters2 <- temp_cluster_names
wt_atlas@cellColData@listData[["Clusters2"]] <- temp_cluster_names
table(wt_atlas$Clusters2,wt_atlas$Sample)

#                     a_E775_1 b_E775_2 c_E825_1 d_E825_2 e_E825_3 f_E825_4 g_E825_5 h_E925_1 i_E925_2 j_E925_3 k_E925_4 l_E925_5 m_E105_1 n_E105_2
# Blood                    2        1        1        1        2        7        0       52      119       38      104      204      181      106
# Cardiomyocyte          161      333      217      195      267      401      202      838      933      721      544      250     1050      856
# Ectoderm               371      406       97      445       60      754      632        2       14       17     1452      797        1        0
# EndMT                    0        0        1        2        0        1        1       17       23       12       14        6      196      212
# Endoderm               412      545      212      400      316      566      390      328     1010      529     1329      781        4      145
# Endothelium             69       80       59       77       74      126       86      196      249      173      229      158      629      396
# Epicardium             101      142       49       87      144      161      114       85      449      271      198      163      437      265
# LPM                     95       92        3       57       88      109       56       69       79      115      110      153        2        0
# Neural_Crest            36       34       10       26       30       74       79       79      455      144      719      308       39       82
# Paraxial_Mesoderm      197      259      111      241       49      404      264       35      217       97      641      597        3       83
# SHF_Progenitor         287      379      152      357      248      584      407      345     1249      583      934      801      196      224
# 
#                     o_E105_3 p_E105_4 q_E115_1 r_E115_2
# Blood                  403      396      154      233
# Cardiomyocyte           84       87      182      281
# Ectoderm              1042      989       63      128
# EndMT                   29       23       43       84
# Endoderm               929      668      416      746
# Endothelium            215      152      269      245
# Epicardium             338      316      161      380
# LPM                     79        3        3        6
# Neural_Crest          1788     1554     2087     2857
# Paraxial_Mesoderm      896      731      561     1150
# SHF_Progenitor         871      706     1796     2172

p16 <- plotEmbedding(
  wt_atlas,
  colorBy = "cellColData",
  name = "Clusters2",
  embedding = "UMAPHarmony")

p16

### save all intermediate plots
plotPDF(p16,
        name = "Plot-UMAP-Manual-Annotations.pdf",
        ArchRProj = wt_atlas,
        addDOC = FALSE, width = , height = 5)

saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)
wt_atlas
# numberOfCells(1): 64956
# medianTSS(1): 12.193
# medianFrags(1): 48811.5


########################
########################
########################
# Proceed to step 4


