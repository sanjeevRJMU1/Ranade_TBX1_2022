markerGenes_spotmarkerGenes_spot <- c("Tnnt2","Isl1","Osr1","Fgf8","Tbx1","Wt1","Mab21l2","Meox1","Hand2",
                                      "Tdgf1","Rgs5","Fgf10","Irx4","Vsnl1","Cited1","Nr2f1","Foxd1","Foxc2",
                                      "Alx1","Tlx1","Aplnr","Nrg1","Hand1","Tbx18","Foxf1","Bmp2","Rspo3",
                                      "Dlx1","Dlx2","Dlx3","Dlx4","Dlx5","Dlx6",
                                      "Tdgf1","Irx3","Irx1","Tbx2","Tbx3","Nkx2-5","Gata4",
                                      "Lix1","Tcf21","Lhx1","Lhx2","Gbx1","Otx2","Sox10","Sox2","Myf5","Pax3","Pitx2",
                                      "Twist1","Prrx1","Prrx2","Col1a1","Pdgfra","Ebf1","Ebf2","Lum","Sox9","Msx1","Msx2",
                                      "Pax1","Pax9","Pax5","Wnt4","Hoxa4","Hoxb4","Hoxc4","Hoxd4","Bdnf","Asb4","Zic1","Foxg1","Foxp2")


mesoderm_subset <- addImputeWeights(mesoderm_subset, reducedDims = "Harmony", dimsToUse = 2:20)
markerGenes_spotmarkerGenes_spot_2 <- c("Tagln","Acta2","Des","Vcl")
p14 <- plotEmbedding(
  ArchRProj = mesoderm_subset, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes_spotmarkerGenes_spot_2, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(mesoderm_subset)
)
plotPDF(plotList = p14, 
        name = "Plot-UMAP-Marker-Genes-SpotCheck_2.pdf", 
        ArchRProj = mesoderm_subset, 
        addDOC = FALSE, width = 5, height = 5)



## My manual labeling --> this goes to Clusters6 as rest are already taken
library(plyr)
temp_cluster_names <- as.character(revalue(mesoderm_subset$Clusters, c(
  "C1"="Cardiomyocyte",
  "C2"="Cardiomyocyte",
  "C3"="Cardiopharyngeal_Mesenchyme",
  "C4"="Cranial_Mesenchyme",
  "C5"="Neural_Crest",
  "C6"="Neural_Crest",
  "C7"="Neural_Crest",
  "C8"="Epicardium",
  "C9"="Smooth_Muscle",
  "C10"="Juxtacardiac_Field",
  "C11"="Cardiopharyngeal_Mesenchyme",
  "C12"="Paraxial_Mesoderm",
  "C13"="Cardiopharyngeal_Mesoderm",
  "C14"="Paraxial_Mesoderm",
  "C15"="Posterior_Second_Heart_Field",
  "C16"="Anterior_Second_Heart_Field"
)))
length(unique(temp_cluster_names))
mesoderm_subset$Clusters6 <- temp_cluster_names
mesoderm_subset@cellColData@listData[["Clusters6"]] <- temp_cluster_names
table(mesoderm_subset$Clusters6,mesoderm_subset$Sample)

################################################################

# Break up each cluster into WT and KO, then add this as a new cellcoldata column called Clusters7
condition.vec.01 <- substr(mesoderm_subset$Sample,5,8)
unique(condition.vec.01)
condition.vec.02 <- as.vector(mesoderm_subset$Clusters6)
condition.vec <- paste(condition.vec.02, condition.vec.01, sep = "_")
table(condition.vec)
# condition.vec
# Anterior_Second_Heart_Field_05_K  Anterior_Second_Heart_Field_05_W  Anterior_Second_Heart_Field_15_K 
# 200                               296                               437 
# Anterior_Second_Heart_Field_15_W  Anterior_Second_Heart_Field_25_K  Anterior_Second_Heart_Field_25_W 
# 597                               286                               389 
# Cardiomyocyte_05_K                Cardiomyocyte_05_W                Cardiomyocyte_15_K 
# 255                               163                               778 
# Cardiomyocyte_15_W                Cardiomyocyte_25_K                Cardiomyocyte_25_W 
# 358                               658                               788 
# Cardiopharyngeal_Mesenchyme_05_K  Cardiopharyngeal_Mesenchyme_05_W  Cardiopharyngeal_Mesenchyme_15_K 
# 582                               678                              2284 
# Cardiopharyngeal_Mesenchyme_15_W  Cardiopharyngeal_Mesenchyme_25_K  Cardiopharyngeal_Mesenchyme_25_W 
# 1243                               563                               608 
# Cardiopharyngeal_Mesoderm_05_K    Cardiopharyngeal_Mesoderm_05_W    Cardiopharyngeal_Mesoderm_15_K 
# 872                               741                              1284 
# Cardiopharyngeal_Mesoderm_15_W    Cardiopharyngeal_Mesoderm_25_K    Cardiopharyngeal_Mesoderm_25_W 
# 788                               565                               565 
# Cranial_Mesenchyme_05_K           Cranial_Mesenchyme_05_W           Cranial_Mesenchyme_15_K 
# 500                               511                              1525 
# Cranial_Mesenchyme_15_W           Cranial_Mesenchyme_25_K           Cranial_Mesenchyme_25_W 
# 1057                               289                               356 
# Epicardium_05_K                   Epicardium_05_W                   Epicardium_15_K 
# 779                               565                              1028 
# Epicardium_15_W                   Epicardium_25_K                   Epicardium_25_W 
# 212                               551                               460 
# Juxtacardiac_Field_05_K           Juxtacardiac_Field_05_W           Juxtacardiac_Field_15_K 
# 568                               735                               495 
# Juxtacardiac_Field_15_W           Juxtacardiac_Field_25_K           Juxtacardiac_Field_25_W 
# 1174                               110                                90 
# Neural_Crest_05_K                 Neural_Crest_05_W                 Neural_Crest_15_K 
# 2071                              2213                               998 
# Neural_Crest_15_W                 Neural_Crest_25_K                 Neural_Crest_25_W 
# 3131                               370                               308 
# Paraxial_Mesoderm_05_K            Paraxial_Mesoderm_05_W            Paraxial_Mesoderm_15_K 
# 515                               604                              2358 
# Paraxial_Mesoderm_15_W            Paraxial_Mesoderm_25_K            Paraxial_Mesoderm_25_W 
# 647                               382                               556 
# Posterior_Second_Heart_Field_05_K Posterior_Second_Heart_Field_05_W Posterior_Second_Heart_Field_15_K 
# 935                               550                              3163 
# Posterior_Second_Heart_Field_15_W Posterior_Second_Heart_Field_25_K Posterior_Second_Heart_Field_25_W 
# 1667                               696                               605 
# Smooth_Muscle_05_K                Smooth_Muscle_05_W                Smooth_Muscle_15_K 
# 62                                64                               463 
# Smooth_Muscle_15_W                Smooth_Muscle_25_K                Smooth_Muscle_25_W 
# 276                                49                                81
mesoderm_subset$Clusters7 <- as.character(condition.vec)
# spot check
table(mesoderm_subset$Clusters7)

p16 <- plotEmbedding(ArchRProj = mesoderm_subset, colorBy = "cellColData", name = "Clusters6", embedding = "UMAPHarmony")


p16

plotPDF(p16,
        name = "Plot-UMAP-Manual-Annotations.pdf",
        ArchRProj = mesoderm_subset,
        addDOC = FALSE, width = , height = 5)


saveArchRProject(ArchRProj = mesoderm_subset, outputDirectory = "/Users/sranade/scATAC-seq/Mesoderm_subset3_210517", load = TRUE)





