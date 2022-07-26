###### Tbx1 WT v KO all aggr
###### E9.5 - E11.5, 13 samples total
##### SR38 = E9.5, SR40 = E10.5, SR43 = E11.5
### Chapters 1 - 3: Setup

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 1.0.1

#### 1a. Setup
#
addArchRGenome("mm10")
inputFiles <- c("/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR38A_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR38A_2/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR38A_3/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR38A_4/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR40A_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR40A_2/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR40A_3/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR40A_4/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR43_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR43_2/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR43_3/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR43_4/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/Tbx1_all_aggr/all_fragments_files/SR43_5/fragments.tsv.gz")
is.character(inputFiles)
length(inputFiles)
sampleNames <- c("a_E925_WT_1",
                 "b_E925_KO_1",
                 "c_E925_WT_2",
                 "d_E925_KO_2",
                 "e_E105_KO_1",
                 "f_E105_KO_2",
                 "g_E105_WT_1",
                 "h_E105_WT_2",
                 "i_E115_KO_1",
                 "j_E115_WT_1",
                 "k_E115_WT_2",
                 "l_E115_KO_2",
                 "m_E115_KO_3")
is.character(sampleNames)
length(sampleNames)

#### 1b. Create Arrow Files (with modifications)
setwd("/Users/sranade/scATAC-seq/Tbx1_all_aggr")
starttime <- Sys.time()
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 6, 
  minFrags = 10000,
  maxFrags = 5e+06,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  verbose = T,
  force = T
)
ArrowFiles
endtime <- Sys.time()
## clock
endtime - starttime
#Time difference of 4.416177 hours

#### 2. Doublet Inference with ArchR 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#### 3a. Creating an ArchRProject
Tbx1_all_aggr  <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_all_aggr",
  copyArrows = T
)
Tbx1_all_aggr
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/Tbx1_all_aggr 
# samples(13): g_E105_WT_1 e_E105_KO_1 ... i_E115_KO_1 l_E115_KO_2
# sampleColData names(1): ArrowFiles
# cellColData names(15): Sample TSSEnrichment ... DoubletEnrichment BlacklistRatio
# numberOfCells(1): 83874
# medianTSS(1): 12.175
# medianFrags(1): 48018

saveArchRProject(ArchRProj = Tbx1_all_aggr, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_all_aggr", load = TRUE)

### close R studio, delete fragments files and the initial arrow files, reload project
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 

Tbx1_all_aggr <- loadArchRProject(path = "/Users/sranade/scATAC-seq/Tbx1_all_aggr")
Tbx1_all_aggr
### continue on with analysis

### 3b. Sample Statistics plots for ArchRProject
## ridge plot of TSS enrichment per sample
p1 <- plotGroups(
  ArchRProj = Tbx1_all_aggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

## violin plot for each sample for the TSS enrichment scores
p2 <- plotGroups(
  ArchRProj = Tbx1_all_aggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

## ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
  ArchRProj = Tbx1_all_aggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

## violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
  ArchRProj = Tbx1_all_aggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

## save plots
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = Tbx1_all_aggr, addDOC = FALSE, width = 4, height = 4)

### 3c. Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
p5 <- plotFragmentSizes(ArchRProj = Tbx1_all_aggr)

p6 <- plotTSSEnrichment(ArchRProj = Tbx1_all_aggr)

## save plots
plotPDF(p5,p6, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = Tbx1_all_aggr, addDOC = FALSE, width = 5, height = 5)

### 3d. Filtering Doublets from an ArchRProject
Tbx1_all_aggr <- filterDoublets(Tbx1_all_aggr)
# Filtering 5812 cells from ArchRProject!
#   g_E105_WT_1 : 516 of 7190 (7.2%)
# e_E105_KO_1 : 349 of 5909 (5.9%)
# h_E105_WT_2 : 357 of 5982 (6%)
# f_E105_KO_2 : 562 of 7501 (7.5%)
# m_E115_KO_3 : 1064 of 10317 (10.3%)
# k_E115_WT_2 : 830 of 9112 (9.1%)
# j_E115_WT_1 : 373 of 6108 (6.1%)
# a_E925_WT_1 : 452 of 6726 (6.7%)
# d_E925_KO_2 : 405 of 6365 (6.4%)
# c_E925_WT_2 : 194 of 4412 (4.4%)
# b_E925_KO_1 : 174 of 4172 (4.2%)
# i_E115_KO_1 : 390 of 6246 (6.2%)
# l_E115_KO_2 : 146 of 3834 (3.8%)
Tbx1_all_aggr
# numberOfCells(1): 60427
# medianTSS(1): 11.851
# medianFrags(1): 53758

# Save project
## In this one instance alone, use dropCells = T. Once you do this, you don't need to use it any more and actually using it
## might fuck up the matrix sizes downstream. 
saveArchRProject(ArchRProj = Tbx1_all_aggr, outputDirectory = "/Users/sranade/scATAC-seq/Tbx1_all_aggr", load = TRUE, dropCells = T)
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/Tbx1_all_aggr 
# samples(13): g_E105_WT_1 e_E105_KO_1 ... i_E115_KO_1 l_E115_KO_2
# sampleColData names(1): ArrowFiles
# cellColData names(15): Sample TSSEnrichment ... DoubletEnrichment BlacklistRatio
# numberOfCells(1): 78062
# medianTSS(1): 12.209
# medianFrags(1): 47352.5


####################################################
## Proceed to 2. Clustering script





