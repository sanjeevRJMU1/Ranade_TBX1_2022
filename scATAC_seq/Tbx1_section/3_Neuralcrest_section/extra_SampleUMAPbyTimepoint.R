## label each time point as a sample and plot umap with 3 time points

# Break up each cluster into timepoint and add as a metadata column as Clusters3
table(CNC_subset$Sample)
# a_E925_WT_1 b_E925_KO_1 c_E925_WT_2 d_E925_KO_2 e_E105_KO_1 f_E105_KO_2 g_E105_WT_1 h_E105_WT_2 i_E115_KO_1 
# 322         178          76         302        1110        1529        1522        1426         516 
# j_E115_WT_1 k_E115_WT_2 l_E115_KO_2 m_E115_KO_3 
# 1796        2509         299         678 

condition.vec.01 <- substr(CNC_subset$Sample,3,6)

unique(condition.vec.01)
#[1] "E105" "E115" "E925"
table(condition.vec.01)
# condition.vec.01
# E105 E115 E925 
# 5587 5798  878 

CNC_subset$Clusters8 <- as.character(condition.vec.01)
# spot check
table(CNC_subset$Clusters8)

p16 <- plotEmbedding(ArchRProj = CNC_subset, colorBy = "cellColData", name = "Clusters8", embedding = "UMAPHarmony")


p16

plotPDF(p16,
        name = "Plot-UMAP-Manual-Annotations-by-Timepoint.pdf",
        ArchRProj = CNC_subset,
        addDOC = FALSE, width = , height = 5)




