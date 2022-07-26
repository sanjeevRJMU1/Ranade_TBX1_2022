
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

# load object
wt_atlas <- loadArchRProject(path = "/Users/sranade/scATAC-seq/new_wt_section_211111")
wt_atlas


#### chromvar plots for supp fig 1:
# extract subset of motifs for downstream analysis
motifs <- c("Gata4","Tbx20","Pitx1","Foxa1","ERG","Hoxa2")
markerMotifs <- getFeatures(wt_atlas, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
### This is silly, you can simply give it a chr vector list of motifs (from the exported homer motif matrix)
### and just add z: to the full name of homer motif

# motifs to use:
# blood: GATA.SCL.Zf.bHLH_111
# CM: Mef2d.MADS_164
# ectoderm: Sox2.HMG_261
# EndMT: Atf2.bZIP_11
# Endoderm: Foxa2.Forkhead_89
# endothelium: ERG.ETS_72
# epicardium: TEAD.TEA_294
# LPM: Gata2.Zf_102
# neural crest: RXR.NR..DR1_251
# parax: Ascl1.bHLH_9
# SHF: HOXA2.Homeobox_128




markerMotifs <- c("z:Pax7.Paired.Homeobox..longest_213")

########## 
## Feature plots of motif deviations. VERY USEFUL!!
########## 
# p <- plotEmbedding(
#   ArchRProj = wt_atlas, 
#   colorBy = "MotifMatrix", 
#   name = sort(markerMotifs), 
#   embedding = "UMAPHarmony",
#   imputeWeights = getImputeWeights(wt_atlas)
# )
# p


wt_atlas <- addMotifAnnotations(ArchRProj = wt_atlas, motifSet = "homer", name = "Motif", force = T)

markerMotifs <- c("z:Hoxb1_423")
p <- plotGroups(ArchRProj = wt_atlas, 
                groupBy = "Clusters2", 
                colorBy = "MotifMatrix", 
                maxCells = 25000,
                name = markerMotifs,
                plotAs = "violin",
                imputeWeights = getImputeWeights(wt_atlas)
)
p







df_shf <- p$data
df_shf <- df_shf[df_shf$x == "SHF_Progenitor",]


plotPDF(p, name = "Motif_chromvar_SHF", width = 5, height = 5, ArchRProj = wt_atlas, addDOC = FALSE)


############### get motifs form chromvar
library(ArchR)
library(chromVARmotifs)

data("homer_pwms") #motifs from HOMER

PWMs <- homer_pwms

PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range
#[1] 0.9999996 1.0000004

#Maybe we can just tidy this up a tiny bit

PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range
#[1] 1 1















### ggseqlogo
require(ggplot2)
require(ggseqlogo)

## merged AHF and pharyngeal Mesoderm outputs
setwd("/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs")

# blod
motif_bach1<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif137_blood.motif",sep = '\t',header = T)
motif_bach1<-t(motif_bach1)
rownames(motif_bach1)<-c()
motif_bach1<-motif_bach1[1:4,]
motif_bach1<-as.matrix(motif_bach1)
rownames(motif_bach1)<-c("A","C","G","T")
ggseqlogo(motif_bach1, method = 'prob')
#save as pdf

# CM
motif_Atf2<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif214_CM.motif",sep = '\t',header = T)
motif_Atf2<-t(motif_Atf2)
rownames(motif_Atf2)<-c()
motif_Atf2<-motif_Atf2[1:4,]
motif_Atf2<-as.matrix(motif_Atf2)
rownames(motif_Atf2)<-c("A","C","G","T")
ggseqlogo(motif_Atf2, method = 'prob')
#save as pdf

# ectoderm
motif_Erg<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif338_ectoderm.motif",sep = '\t',header = T)
motif_Erg<-t(motif_Erg)
rownames(motif_Erg)<-c()
motif_Erg<-motif_Erg[1:4,]
motif_Erg<-as.matrix(motif_Erg)
rownames(motif_Erg)<-c("A","C","G","T")
ggseqlogo(motif_Erg, method = 'prob')
#save as pdf

# EndMT
motif_Fox<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif12_endmt.motif",sep = '\t',header = T)
motif_Fox<-t(motif_Fox)
rownames(motif_Fox)<-c()
motif_Fox<-motif_Fox[1:4,]
motif_Fox<-as.matrix(motif_Fox)
rownames(motif_Fox)<-c("A","C","G","T")
ggseqlogo(motif_Fox, method = 'prob')
#save as pdf

# Endoderm
motif_Tcf3<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif111_endoderm.motif",sep = '\t',header = T)
motif_Tcf3<-t(motif_Tcf3)
rownames(motif_Tcf3)<-c()
motif_Tcf3<-motif_Tcf3[1:4,]
motif_Tcf3<-as.matrix(motif_Tcf3)
rownames(motif_Tcf3)<-c("A","C","G","T")
ggseqlogo(motif_Tcf3, method = 'prob')
#save as pdf

# endothelium
motif_hoxa2<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif92_endothelium.motif",sep = '\t',header = T)
motif_hoxa2<-t(motif_hoxa2)
rownames(motif_hoxa2)<-c()
motif_hoxa2<-motif_hoxa2[1:4,]
motif_hoxa2<-as.matrix(motif_hoxa2)
rownames(motif_hoxa2)<-c("A","C","G","T")
ggseqlogo(motif_hoxa2, method = 'prob')
#save as pdf

# epicardium
motif_mef2c<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif380_epicardium.motif",sep = '\t',header = T)
motif_mef2c<-t(motif_mef2c)
rownames(motif_mef2c)<-c()
motif_mef2c<-motif_mef2c[1:4,]
motif_mef2c<-as.matrix(motif_mef2c)
rownames(motif_mef2c)<-c("A","C","G","T")
ggseqlogo(motif_mef2c, method = 'prob')
#save as pdf

# LPM
motif_rxr<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif128_LPM.motif",sep = '\t',header = T)
motif_rxr<-t(motif_rxr)
rownames(motif_rxr)<-c()
motif_rxr<-motif_rxr[1:4,]
motif_rxr<-as.matrix(motif_rxr)
rownames(motif_rxr)<-c("A","C","G","T")
ggseqlogo(motif_rxr, method = 'prob')
#save as pdf

# neural crest
motif_sox2<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif323_neuralcrest.motif",sep = '\t',header = T)
motif_sox2<-t(motif_sox2)
rownames(motif_sox2)<-c()
motif_sox2<-motif_sox2[1:4,]
motif_sox2<-as.matrix(motif_sox2)
rownames(motif_sox2)<-c("A","C","G","T")
ggseqlogo(motif_sox2, method = 'prob')
#save as pdf

# paraxmeso
motif_tead<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif9_parax.motif",sep = '\t',header = T)
motif_tead<-t(motif_tead)
rownames(motif_tead)<-c()
motif_tead<-motif_tead[1:4,]
motif_tead<-as.matrix(motif_tead)
rownames(motif_tead)<-c("A","C","G","T")
ggseqlogo(motif_tead, method = 'prob')
#save as pdf

# SHF
motif_tead<-read.csv(file = "/Users/sranade/scATAC-seq/new_wt_section_211111/Exports/motifs/motif166_SHF.motif",sep = '\t',header = T)
motif_tead<-t(motif_tead)
rownames(motif_tead)<-c()
motif_tead<-motif_tead[1:4,]
motif_tead<-as.matrix(motif_tead)
rownames(motif_tead)<-c("A","C","G","T")
ggseqlogo(motif_tead, method = 'prob')
#save as pdf



#############################################################################################################################
# re-make figure to focus on SHF prog and subsets with unique sets of enrichment
# Ascl1.bHLH_9
# HOXA2.Homeobox_128
# Pitx1.Ebox.Homeobox.bHLH_224
# Tcf21.bHLH_287
# Atoh1.bHLH_15
# ZBTB18.Zf_307
# RXR.NR..DR1_251
# Lhx2.Homeobox_151
# Pax7.Paired.Homeobox..longest_213
# PPARE.NR..DR1_227
# Phox2a.Homeobox_221

markerMotifs <- c("z:Hoxb4.Homeobox_131")

########## 
## Feature plots of motif deviations. VERY USEFUL!!
########## 
p <- plotEmbedding(
  ArchRProj = wt_atlas, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(wt_atlas)
)
p



plotPDF(p, name = "Motif_chromvar_hoxb4", width = 5, height = 5, ArchRProj = wt_atlas, addDOC = FALSE)



#### motifs
setwd("/Users/sranade/Dropbox (Gladstone)/scATAC_Manuscript_v2/Figure_sections/Figure_1/chromvar_plots/Chromvar_new_220316/motifs")

# ascl1
motif_ascl1<-read.csv(file = "ascl1.motif",sep = '\t',header = T)
motif_ascl1<-t(motif_ascl1)
rownames(motif_ascl1)<-c()
motif_ascl1<-motif_ascl1[1:4,]
motif_ascl1<-as.matrix(motif_ascl1)
rownames(motif_ascl1)<-c("A","C","G","T")
ggseqlogo(motif_ascl1, method = 'prob')
#save as pdf


# atoh
motif_atoh<-read.csv(file = "atoh.motif",sep = '\t',header = T)
motif_atoh<-t(motif_atoh)
rownames(motif_atoh)<-c()
motif_atoh<-motif_atoh[1:4,]
motif_atoh<-as.matrix(motif_atoh)
rownames(motif_atoh)<-c("A","C","G","T")
ggseqlogo(motif_atoh, method = 'prob')
#save as pdf

# hoxa2
motif_hoxa2<-read.csv(file = "hoxa2.motif",sep = '\t',header = T)
motif_hoxa2<-t(motif_hoxa2)
rownames(motif_hoxa2)<-c()
motif_hoxa2<-motif_hoxa2[1:4,]
motif_hoxa2<-as.matrix(motif_hoxa2)
rownames(motif_hoxa2)<-c("A","C","G","T")
ggseqlogo(motif_hoxa2, method = 'prob')
#save as pdf

# hoxb4
motif_hoxb4<-read.csv(file = "hoxb4.motif",sep = '\t',header = T)
motif_hoxb4<-t(motif_hoxb4)
rownames(motif_hoxb4)<-c()
motif_hoxb4<-motif_hoxb4[1:4,]
motif_hoxb4<-as.matrix(motif_hoxb4)
rownames(motif_hoxb4)<-c("A","C","G","T")
ggseqlogo(motif_hoxb4, method = 'prob')
#save as pdf

# lhx2
motif_lhx2<-read.csv(file = "lhx2.motif",sep = '\t',header = T)
motif_lhx2<-t(motif_lhx2)
rownames(motif_lhx2)<-c()
motif_lhx2<-motif_lhx2[1:4,]
motif_lhx2<-as.matrix(motif_lhx2)
rownames(motif_lhx2)<-c("A","C","G","T")
ggseqlogo(motif_lhx2, method = 'prob')
#save as pdf

# meis1
motif_meis1<-read.csv(file = "meis1.motif",sep = '\t',header = T)
motif_meis1<-t(motif_meis1)
rownames(motif_meis1)<-c()
motif_meis1<-motif_meis1[1:4,]
motif_meis1<-as.matrix(motif_meis1)
rownames(motif_meis1)<-c("A","C","G","T")
ggseqlogo(motif_meis1, method = 'prob')
#save as pdf


# pax7
motif_pax7<-read.csv(file = "pax7.motif",sep = '\t',header = T)
motif_pax7<-t(motif_pax7)
rownames(motif_pax7)<-c()
motif_pax7<-motif_pax7[1:4,]
motif_pax7<-as.matrix(motif_pax7)
rownames(motif_pax7)<-c("A","C","G","T")
ggseqlogo(motif_pax7, method = 'prob')
#save as pdf

# pitx1
motif_pitx1<-read.csv(file = "pitx1.motif",sep = '\t',header = T)
motif_pitx1<-t(motif_pitx1)
rownames(motif_pitx1)<-c()
motif_pitx1<-motif_pitx1[1:4,]
motif_pitx1<-as.matrix(motif_pitx1)
rownames(motif_pitx1)<-c("A","C","G","T")
ggseqlogo(motif_pitx1, method = 'prob')
#save as pdf

# rbpj
motif_rbpj<-read.csv(file = "rbpj.motif",sep = '\t',header = T)
motif_rbpj<-t(motif_rbpj)
rownames(motif_rbpj)<-c()
motif_rbpj<-motif_rbpj[1:4,]
motif_rbpj<-as.matrix(motif_rbpj)
rownames(motif_rbpj)<-c("A","C","G","T")
ggseqlogo(motif_rbpj, method = 'prob')
#save as pdf

# rxr
motif_rxr<-read.csv(file = "rxr.motif",sep = '\t',header = T)
motif_rxr<-t(motif_rxr)
rownames(motif_rxr)<-c()
motif_rxr<-motif_rxr[1:4,]
motif_rxr<-as.matrix(motif_rxr)
rownames(motif_rxr)<-c("A","C","G","T")
ggseqlogo(motif_rxr, method = 'prob')
#save as pdf

# Zbtb18
motif_Zbtb18<-read.csv(file = "Zbtb18.motif",sep = '\t',header = T)
motif_Zbtb18<-t(motif_Zbtb18)
rownames(motif_Zbtb18)<-c()
motif_Zbtb18<-motif_Zbtb18[1:4,]
motif_Zbtb18<-as.matrix(motif_Zbtb18)
rownames(motif_Zbtb18)<-c("A","C","G","T")
ggseqlogo(motif_Zbtb18, method = 'prob')
#save as pdf


###################################################################################################
#################################
# Subset cells based on chromvar motifs

### subset only SHF progenitor cells
table(wt_atlas$Clusters2)
shf_prog <- wt_atlas
## subset
# assign Y/N values to only the clusters you want
library(plyr)
temp_cluster_names <- as.character(revalue(shf_prog$Clusters2, c(
  "Blood"="False",
  "Cardiomyocyte"="False",
  "Ectoderm"="False",
  "EndMT"="False",
  "Endoderm"="False",
  "Endothelium"="False",
  "Epicardium"="False",
  "LPM"="False",
  "Neural_Crest"="False",
  "Paraxial_Mesoderm"="False",
  "SHF_Progenitor"="True"
)))
length(unique(temp_cluster_names))
shf_prog$Clusters3 <- temp_cluster_names
shf_prog@cellColData@listData[["Clusters3"]] <- temp_cluster_names
table(shf_prog$Clusters3,shf_prog$Sample)
## subset
# idxSample <- BiocGenerics::which(shf_prog$Clusters3 %in% "True")
# cellsSample <- shf_prog$cellNames[idxSample]
# subsetArchRProject(
#   ArchRProj = neuralcrest_subset,
#   cells = cellsSample,
#   outputDirectory = "/Users/sranade/scATAC-seq/shf_prog", force = T
# )


idxSample <- BiocGenerics::which(shf_prog$Clusters2 == "SHF_Progenitor")
cellsSample <- shf_prog$cellNames[idxSample]
subsetArchRProject(
  ArchRProj = shf_prog,
  cells = cellsSample,
  outputDirectory = "/Users/sranade/scATAC-seq/shf_prog", force = T
)
shf_prog <- loadArchRProject(path = "/Users/sranade/scATAC-seq/shf_prog")
shf_prog

### check available matrices
getAvailableMatrices(wt_atlas)

### pull motif matrix from ArchR object
motif_matrix <- getMatrixFromProject(ArchRProj = shf_prog,
                                     useMatrix = "MotifMatrix",
                                     useSeqnames = NULL,
                                     verbose = TRUE,
                                     binarize = FALSE,
                                     threads = getArchRThreads(),
                                     logFile = createLogFile("getMatrixFromProject")
)


### let's dig around and find the deviations matrix
rownames(motif_matrix@assays@data@listData$deviations)


### subset to specific row(s)
## have to find+check all relevant rows
subset_motifMtx <- motif_matrix@assays@data@listData$deviations[rownames(motif_matrix@assays@data@listData$deviations) %like% ".bHLH_9",]
#subset_motifMtx <- motif_matrix@assays@data@listData$deviations[rownames(motif_matrix@assays@data@listData$deviations) == "Ascl1.bHLH_9",]
#subset_motifMtx <- motif_matrix@assays@data@listData$deviations[rownames((motif_matrix@assays@data@listData$deviations) %like% "Fox"),]

rownames(subset_motifMtx)
# [1] "Ascl1.bHLH_9"              "Fox.Ebox.Forkhead.bHLH_90"
subset_motifMtx<- subset_motifMtx[1,]
names(subset_motifMtx)


## ensure this is a named num - deviation values from selected row, named per cell, as data frame
motifDeviations <- data.frame(subset_motifMtx)
## check range if needed
head(motifDeviations)
range(motifDeviations$subset_motifMtx)
motifDeviations$quartile<- ntile(motifDeviations$subset_motifMtx, 10)  

## pull cell names > threshold
# in 10th qtile or top 10%
aboveThresh_cells <- rownames(motifDeviations[motifDeviations$quartile ==10,])


### add metadata
## set label to be below, by default
shf_prog$deviationThreshold <- "Ascl_low"
## flag cells > thresh
shf_prog$deviationThreshold[rownames(shf_prog) %in% aboveThresh_cells] <- "Ascl_high"
## check cell count
table(shf_prog$deviationThreshold)
# Ascl_high  Ascl_low 
# 1229     11062 

# check before DAR
table(shf_prog$Clusters2)

### go on, run your DARs, run your enrichment
markerTest <- getMarkerFeatures(
  ArchRProj = shf_prog, 
  useMatrix = "PeakMatrix",
  groupBy = "deviationThreshold",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Ascl_high",
  bgdGroups = "Ascl_low"
)

pma <- plotMarkers(seMarker = markerTest, name = "Ascl_high", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
pma

shf_prog <- addImputeWeights(shf_prog)

p <- plotGroups(ArchRProj = shf_prog, 
                groupBy = "deviationThreshold", 
                colorBy = "MotifMatrix", 
                name = "z:Atoh1.bHLH_15",
                imputeWeights = getImputeWeights(shf_prog)
)
p

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = shf_prog,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
motifsUp

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = shf_prog,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)


df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo


###
markerList_Ascl1_SHF <- getMarkers(markerTest, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList_Ascl1_SHF_df <- as.data.frame(markerList_Ascl1_SHF@listData[["Ascl_high"]])
markerList_Ascl1_SHF_df_homer <- markerList_Ascl1_SHF_df[order(markerList_Ascl1_SHF_df$Log2FC, decreasing = TRUE),]
markerList_Ascl1_SHF_df_homer$rank <- seq_len(nrow(markerList_Ascl1_SHF_df_homer))
markerList_Ascl1_SHF_df_homer$blank <- ""
markerList_Ascl1_SHF_df_homer_trim <- markerList_Ascl1_SHF_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Ascl1_SHF_df_homer_trim, file = "/Users/sranade/Dropbox (Gladstone)/scATAC_Manuscript_v2/Figure_sections/Figure_1/Ascl1_section/markerList_Ascl1_SHF_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


# 
markerList_Ascl1_SHF_up <- getMarkers(markerTest, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
markerList_Ascl1_SHF_up_df <- as.data.frame(markerList_Ascl1_SHF_up@listData[["Ascl_high"]])
markerList_Ascl1_SHF_up_df_homer <- markerList_Ascl1_SHF_up_df[order(markerList_Ascl1_SHF_up_df$Log2FC, decreasing = F),]
markerList_Ascl1_SHF_up_df_homer$rank <- seq_len(nrow(markerList_Ascl1_SHF_up_df_homer))
markerList_Ascl1_SHF_up_df_homer$blank <-""
markerList_Ascl1_SHF_up_df_homer_trim <- markerList_Ascl1_SHF_up_df_homer[,c(1,2,3,9,10,5)]
write.table(markerList_Ascl1_SHF_up_df_homer_trim, file = "/Users/sranade/Dropbox (Gladstone)/scATAC_Manuscript_v2/Figure_sections/Figure_1/Ascl1_section/markerList_Ascl1_SHF_up_df_homer_trim.bed",
            quote = F, sep = '\t',col.names = F, row.names = F)


###################################################################################################################
###################################################################################################################
# get predicted cells for Ascl1 high vs low and run DGE in seurat

### Pull scRNAseq cell barcodes associated with high or low deviation scores and run DGE...
## load seurat object
library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
scRNA <- readRDS(file = "/Users/sranade/scRNA-seq/WT_aggr_Fig1_new_211119/wt_atlas_rna_annotated.RDS")

### create Seurat metadata to add
## pull unique scRNA cell barcodes
highDeviationCells <- unique(shf_prog@cellColData@listData$predictedCell[shf_prog$deviationThreshold == "Ascl_high"])
# 626 # 405 #227
lowDeviationCells <- unique(shf_prog@cellColData@listData$predictedCell[shf_prog$deviationThreshold == "Ascl_low"])
#1600 #1680 #1741
## exclude cells from low that are in high
# forgot to use "negated" %in% operator when we were on the phone. PEBCAK.
`%!in%` <- Negate(`%in%`)
lowDeviationCells <- lowDeviationCells[lowDeviationCells %!in% highDeviationCells]
# 1169 # 1390 #1568
# sigma = 431 #sigma = 290 #sigma = 173

## pull all cells from Seurat obj
metadata2add <- data.frame(scRNA_barcode = colnames(scRNA), deviation_value = "not_mapped_to_ATAC")
## annotate deviations
metadata2add$deviation_value[metadata2add$scRNA_barcode %in% highDeviationCells] <- "high"
metadata2add$deviation_value[metadata2add$scRNA_barcode %in% lowDeviationCells] <- "low"
## set colnames, drop barcode col
rownames(metadata2add) <- metadata2add$scRNA_barcode
metadata2add$scRNA_barcode <- NULL

### add metadata
scRNA <- AddMetaData(scRNA, metadata = metadata2add)

### dimplot
DimPlot(scRNA, group.by = "deviation_value")

Idents(scRNA) <- "deviation_value"

### DE
## sending it
highVlow_deviationCells <- FindMarkers(scRNA, ident.1 = "high", ident.2 = "low")
## boyscout edition
highVlow_boyscoutEdition <- FindMarkers(scRNA, ident.1 = "high", ident.2 = "low", max.cells.per.ident = min(table(Idents(scRNA))))

## remove not mapped cells
table(scRNA@active.ident)
# low high 
# 1568  227 

scRNA_trim<-subset(x = scRNA, idents = c("not_mapped_to_ATAC"), invert = TRUE)
scRNA_trim@active.ident<-droplevels(scRNA_trim@active.ident)
table(scRNA_trim@active.ident)

my_order <-c("high","low")
scRNA_trim@active.ident <- factor(x = scRNA_trim@active.ident, levels = my_order)




## genes of interest
dge_genes <- c("Msx1","Msx2","Mab21l2","Hand1","Hand2","Actc1","Tbx1","Eya1","Eya4","Tcf21","Six1","Six2","Meox1","Osr1")

pdf(file = "/Users/sranade/Dropbox (Gladstone)/scATAC_Manuscript_v2/Figure_sections/illustrator_files/Fig_parts_220315/fig1/Ascl1hivlow_dotplot.pdf", width = 4, height = 6)
DotPlot(scRNA_trim, features = dge_genes, col.max = 2.5, col.min = -2.5,dot.scale = 6, cols = c("grey","red"))+coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend()
dev.off()

DotPlot(scRNA_trim, features = dge_genes, col.max = 2.5, col.min = -2.5,dot.scale = 6, cols = c("grey","red"))+coord_flip() + RotatedAxis() + xlab('') + ylab('')

#########
VlnPlot(scRNA_trim, c("Ascl1","Atoh1"), pt.size = )







