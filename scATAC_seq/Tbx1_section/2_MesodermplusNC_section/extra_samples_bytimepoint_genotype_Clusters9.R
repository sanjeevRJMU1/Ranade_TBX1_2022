
# Break up each cluster into timepoint and add as a metadata column as Clusters3
table(mesoderm_subset$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 
# 322         178          76         302        1110        1529        1522        1426         516 
# j_E115_WT_1 k_E115_WT_2 l_E115_KO_2 m_E115_KO_3 
# 1796        2509         299         678 

condition.vec.01 <- substr(mesoderm_subset$Sample,3,8)

unique(condition.vec.01)
#[1] "E105_W" "E105_K" "E115_K" "E115_W" "E925_W" "E925_K"
table(condition.vec.01)
# condition.vec.01
# E105_K E105_W E115_K E115_W E925_K E925_W 
# 7339   7120  14813  11150   4519   4806 

mesoderm_subset$Clusters9 <- as.character(condition.vec.01)
# spot check
table(mesoderm_subset$Clusters9)

p16 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Clusters9", embedding = "UMAPHarmony")


p16

plotPDF(p16,
        name = "Plot-UMAP-Manual-Annotations-by-TimexGenotype.pdf",
        ArchRProj = mesoderm_subset,
        addDOC = FALSE, width = , height = 5)
