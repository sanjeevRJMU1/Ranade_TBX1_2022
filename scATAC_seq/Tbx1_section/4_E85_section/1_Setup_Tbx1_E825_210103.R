###### Tbx1 WT v KO, E8.25. n = 3 KO, n = 4 WT. Embryos were ~6-8 somites, post crescent --> linear heart tube stage
### Chapters 1 - 3: Setup

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 0.9.5

#### 1a. Setup
addArchRGenome("mm10")
inputFiles <- c("/Users/sranade/scATAC-seq/SR39_Tbx1_E825/all_fragments_files/1_ko/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/SR39_Tbx1_E825/all_fragments_files/2_wt/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/SR39_Tbx1_E825/all_fragments_files/3_ko/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/SR39_Tbx1_E825/all_fragments_files/4_ko/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/SR39_Tbx1_E825/all_fragments_files/5_E825_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/SR39_Tbx1_E825/all_fragments_files/6_E825_2/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/SR39_Tbx1_E825/all_fragments_files/7_E825_3/fragments.tsv.gz")
is.character(inputFiles)
length(inputFiles)
sampleNames <- c("a_KO_1","b_WT_1","c_KO_2","d_KO_3","e_WT_2","f_WT_3","g_WT_4")
is.character(sampleNames)
length(sampleNames)

#### 1b. Create Arrow Files (with modifications)
starttime <- Sys.time()
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 6, 
  filterFrags = 10000,
  removeFilteredCells = TRUE,
  maxFrags = 8e+05,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  verbose = T,
  force = T
)
ArrowFiles
endtime <- Sys.time()
## clock
endtime - starttime
#Time difference of 1.118187 hours

#### 2. Doublet Inference with ArchR 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

#### 3a. Creating an ArchRProject
Tbx1_E825  <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/Users/sranade/scATAC-seq/SR39_Tbx1_E825",
  copyArrows = T
)
Tbx1_E825
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/SR39_Tbx1_E825 
# samples(7): c_KO_2 b_WT_1 ... f_WT_3 e_WT_2
# sampleColData names(1): ArrowFiles
# cellColData names(15): Sample TSSEnrichment ... DoubletEnrichment BlacklistRatio
# numberOfCells(1): 15607
# medianTSS(1): 11.359
# medianFrags(1): 51317

saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/SR39_Tbx1_E825", load = TRUE)

### close R studio, delete fragments files and the initial arrow files, reload project
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 0.9.5
Tbx1_E825 <- loadArchRProject(path = "/Users/sranade/scATAC-seq/SR39_Tbx1_E825")
Tbx1_E825
### continue on with analysis

### 3b. Sample Statistics plots for ArchRProject
## ridge plot of TSS enrichment per sample
p1 <- plotGroups(
  ArchRProj = Tbx1_E825, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

## violin plot for each sample for the TSS enrichment scores
p2 <- plotGroups(
  ArchRProj = Tbx1_E825, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

## ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
  ArchRProj = Tbx1_E825, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

## violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
  ArchRProj = Tbx1_E825, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

## save plots
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 4, height = 4)

### 3c. Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
p5 <- plotFragmentSizes(ArchRProj = Tbx1_E825)

p6 <- plotTSSEnrichment(ArchRProj = Tbx1_E825)

## save plots
plotPDF(p5,p6, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = Tbx1_E825, addDOC = FALSE, width = 5, height = 5)

### 3d. Filtering Doublets from an ArchRProject
Tbx1_E825 <- filterDoublets(Tbx1_E825)
# c_KO_2 : 182 of 4274 (4.3%)
# b_WT_1 : 108 of 3295 (3.3%)
# a_KO_1 : 26 of 1616 (1.6%)
# d_KO_3 : 52 of 2283 (2.3%)
# g_WT_4 : 16 of 1294 (1.2%)
# f_WT_3 : 37 of 1925 (1.9%)
# e_WT_2 : 8 of 920 (0.9%)
Tbx1_E825
# numberOfCells(1): 15178
# medianTSS(1): 11.38
# medianFrags(1): 51252

# Save project
## In this one instance alone, use dropCells = T. Once you do this, you don't need to use it any more and actually using it
## might fuck up the matrix sizes downstream. 
saveArchRProject(ArchRProj = Tbx1_E825, outputDirectory = "/Users/sranade/scATAC-seq/SR39_Tbx1_E825", load = TRUE, dropCells = T)

# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/SR26_Tbx1_E825_201214 
# samples(6): e_Het_2 a_WT_1 ... d_Het_1 f_Het_3
# sampleColData names(1): ArrowFiles
# cellColData names(15): Sample TSSEnrichment ... DoubletEnrichment BlacklistRatio
# numberOfCells(1): 12346
# medianTSS(1): 11.913
# medianFrags(1): 46715.5


####################################################
## Proceed to 2. Clustering script





