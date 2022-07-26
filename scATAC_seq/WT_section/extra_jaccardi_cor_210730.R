### jaccard index
library(tidyr)
library(RColorBrewer)
#Now let's look at the overlap between the annotations.

#We will calculate a Jaccard index between all possible overlaps. For plotting, we only keep sets showing at least JI >= 0.15 with another set.

#### Jeeves edit: I renamed columuns and changed color of heatmap (& i wanted to show numbers on heatmap)
# requires you to have already run RNA integration on your atac-seq object


atac_wt_atlas_ji <- Reduce(bind_rows,lapply(unique(wt_atlas$Clusters2), function(SR_label) {
  df <- Reduce(bind_rows,lapply(unique(wt_atlas$predictedGroup), function(CCA_label) {
    ji <- sum(wt_atlas$Clusters2 == SR_label & wt_atlas$predictedGroup == CCA_label)/sum(wt_atlas$Clusters2 == SR_label | wt_atlas$predictedGroup ==CCA_label)
    df <- data.frame(CCA_cell_type=CCA_label, JI=ji, stringsAsFactors = F)
    return(df)
  }))
  df$SR_label_cell_type <- SR_label
  return(df)
}))

atac_wt_atlas_ji_mat <- spread(atac_wt_atlas_ji, key = SR_label_cell_type, value = JI)
row.names(atac_wt_atlas_ji_mat) <- atac_wt_atlas_ji_mat$CCA_cell_type
atac_wt_atlas_ji_mat <- atac_wt_atlas_ji_mat[,2:ncol(atac_wt_atlas_ji_mat)]
summary(apply(atac_wt_atlas_ji_mat, 1, max))
summary(apply(atac_wt_atlas_ji_mat, 2, max))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3140  0.5939  0.8069  0.7456  0.9368  0.9655 

pheatmap::pheatmap(atac_wt_atlas_ji_mat, clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Blues")))(100))
i <- apply(atac_wt_atlas_ji_mat, 1, max) >= 0.15
j <- apply(atac_wt_atlas_ji_mat, 2, max) >= 0.15

pdf(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/jaccard_new.pdf", width = 3, height = 3)
pheatmap::pheatmap(atac_wt_atlas_ji_mat[i,j], clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Reds")))(100), display_numbers = T, number_color = "White",treeheight_row=2, treeheight_col=2, fontsize=6)
dev.off()

#Overall, the two annotations agree well. 
