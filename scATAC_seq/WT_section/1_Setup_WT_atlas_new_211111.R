###### WT atlas, E7.75 - E11.5 (for manuscript)
### Chapters 1 - 3: Setup

library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 4) 
# ArchR : Version 0.9.5

#### 1a. Setup
addArchRGenome("mm10")
inputFiles <- c("/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/1_SR31_1_E775_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/2_SR31_2_E775_2/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/3_SR31_3_E825_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/4_SR33_1_E825_2/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/5_SR33_2_E825_3/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/6_SR39A_2_E825_4/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/7_SR39A_5_E825_5/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/8_SR32_1_E925_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/9_SR32_2_E925_2/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/10_SR32_3_E925_3/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/11_SR38A_1_E925_4/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/12_SR38A_3_E925_5/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/13_SR26_1_E105_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/14_SR26_4_E105_2/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/15_SR40A_3_E105_3/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/16_SR40A_4_E105_4/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/17_SR43_2_E115_1/fragments.tsv.gz",
                "/Users/sranade/scATAC-seq/new_wt_section_211111/all_fragments_files/18_SR43_3_E115_2/fragments.tsv.gz"
                )
is.character(inputFiles)
length(inputFiles)
sampleNames <- c("a_E775_1",
                 "b_E775_2",
                 "c_E825_1",
                 "d_E825_2",
                 "e_E825_3",
                 "f_E825_4",
                 "g_E825_5",
                 "h_E925_1",
                 "i_E925_2",
                 "j_E925_3",
                 "k_E925_4",
                 "l_E925_5",
                 "m_E105_1",
                 "n_E105_2",
                 "o_E105_3",
                 "p_E105_4",
                 "q_E115_1",
                 "r_E115_2")
is.character(sampleNames)
length(sampleNames)

#### 1b. Create Arrow Files (with modifications)
starttime <- Sys.time()
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 6, 
  minFrags = 10000,
  maxFrags = 3e+06,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  verbose = T,
  force = T
)
ArrowFiles
endtime <- Sys.time()
## clock
endtime - starttime
# Check Time difference 

#### 2. Doublet Inference with ArchR 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1, force = T
)

#### 3a. Creating an ArchRProject
wt_atlas  <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111",
  copyArrows = T
)
wt_atlas
# class: ArchRProject 
# outputDirectory: /Users/sranade/scATAC-seq/wt_atlas_new_211027 
# samples(16): m_E105_3 n_E105_4 ... b_E825_2 a_E825_1
# sampleColData names(1): ArrowFiles
# cellColData names(15): Sample TSSEnrichment ... DoubletEnrichment BlacklistRatio
# numberOfCells(1): 64409
# medianTSS(1): 12.066
# medianFrags(1): 49373

saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)

# ### close R studio, delete fragments files and the initial arrow files, reload project
# library(ArchR)
# library(dplyr)
# library(SeuratWrappers)
# addArchRThreads(threads = 4) 
# # ArchR : Version 0.9.5
# wt_atlas <- loadArchRProject(path = "/Users/sranade/scATAC-seq/wt_atlas_new_211027")
# wt_atlas
### continue on with analysis

### 3b. Sample Statistics plots for ArchRProject
## ridge plot of TSS enrichment per sample
p1 <- plotGroups(
  ArchRProj = wt_atlas, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

## violin plot for each sample for the TSS enrichment scores
p2 <- plotGroups(
  ArchRProj = wt_atlas, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

## ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
  ArchRProj = wt_atlas, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

## violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
  ArchRProj = wt_atlas, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

## save plots
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = wt_atlas, addDOC = FALSE, width = 4, height = 4)

### 3c. Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
p5 <- plotFragmentSizes(ArchRProj = wt_atlas)

p6 <- plotTSSEnrichment(ArchRProj = wt_atlas)

## save plots
plotPDF(p5,p6, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = wt_atlas, addDOC = FALSE, width = 5, height = 5)

### 3d. Filtering Doublets from an ArchRProject
wt_atlas <- filterDoublets(wt_atlas)
wt_atlas
# numberOfCells(1): 64956
# medianTSS(1): 12.193
# medianFrags(1): 48811.5

# Save project
saveArchRProject(ArchRProj = wt_atlas, outputDirectory = "/Users/sranade/scATAC-seq/new_wt_section_211111", load = TRUE)

####################################################
## Proceed to 2. Clustering script





