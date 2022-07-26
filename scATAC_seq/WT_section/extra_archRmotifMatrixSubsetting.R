###### E105_DM project
#### Subset to NC
library(ArchR)
library(dplyr)
library(data.table)
library(SeuratWrappers)
addArchRThreads(threads = 12) 
set.seed(1)
# ArchR : Version 1.0.1

######## 
# Load ArchR project!
E105_DM_Harmony <- loadArchRProject(path = "~/magic_local/manuscript_ArchR_efforts_mar2022/data/NC/")

### check available matrices
getAvailableMatrices(E105_DM_Harmony)

### ensure MotifMatrix is available
getAvailableMatrices(E105_DM_Harmony)

### add matrix if needed
if("Motif" %ni% names(E105_DM_Harmony@peakAnnotation)){
  E105_DM_Harmony <- addMotifAnnotations(ArchRProj = E105_DM_Harmony, motifSet = "cisbp", name = "Motif")
}

E105_DM_Harmony <- addBgdPeaks(E105_DM_Harmony)

E105_DM_Harmony <- addDeviationsMatrix(
  ArchRProj = E105_DM_Harmony, 
  peakAnnotation = "Motif",
  force = TRUE
)


### pull motif matrix from ArchR object
motif_matrix <- getMatrixFromProject(ArchRProj = E105_DM_Harmony,
                                     useMatrix = "MotifMatrix",
                                     useSeqnames = NULL,
                                     verbose = TRUE,
                                     binarize = FALSE,
                                     threads = getArchRThreads(),
                                     logFile = createLogFile("getMatrixFromProject")
                                     )

# saveRDS(motif_matrix, "../results/NC_motif_matrix_summarizedExperiment.RDS")


### let's dig around and find the deviations matrix
head(motif_matrix@assays@data@listData$deviations)


### subset to specific row(s)
## have to find+check all relevant rows
subset_motifMtx <- motif_matrix@assays@data@listData$deviations[rownames(motif_matrix@assays@data@listData$deviations) %like% "Pitx",]
# rownames(subset_motifMtx)
# "Pitx1_434" "Pitx3_452" "Pitx2_464"


### split cells based on 1 or more rows of interest
## ensure this is a named num - deviation values from selected row, named per cell
motifDeviations <- subset_motifMtx[rownames(subset_motifMtx) == "Pitx1_434",]
## check range if needed
# range(motifDeviations)
## pull cell names > threshold
deviationThreshold <- 0.06 # arbitrary
# greater than
aboveThresh_cells <- names(motifDeviations[motifDeviations > deviationThreshold])
# less than or equal to... if you want that for some reason
# belowThresh_cells <- names(motifDeviations[motifDeviations <= deviationThreshold])


# ### NOTE: need to account for more rows? rinse & repeat ###
motifDeviations2 <- subset_motifMtx[rownames(subset_motifMtx) == "Pitx3_452",]
## pull cell names > threshold
deviationThreshold2 <- 0.06
# greater than
aboveThresh_cells2 <- names(motifDeviations2[motifDeviations2 > deviationThreshold2])
## check intersection of cells
CellsAboveAllThresholds <- intersect(aboveThresh_cells, aboveThresh_cells2)


### NOTE: if you need to run a third row, fourth, and so on... just keep looping through the above code
## can keep common cells to latest row check + all previous w/ this intersect
# CellsAboveAllThresholds <- intersect(CellsAboveAllThresholds, aboveThresh_cells2)


### add metadata
## set label to be below, by default
E105_DM_Harmony$deviationThreshold <- "below_6e-2_Pitx1_434and452"
## flag cells > thresh
E105_DM_Harmony$deviationThreshold[rownames(E105_DM_Harmony) %in% CellsAboveAllThresholds] <- "above_6e-2_Pitx1_434and452"
## check cell count
table(E105_DM_Harmony$deviationThreshold)

### go on, run your DARs, run your enrichment


# saveArchRProject(E105_DM_Harmony, outputDirectory = "~/magic_local/manuscript_ArchR_efforts_mar2022/data/NCsubset_withMotifMatrixandMetadataAdded_mar172022")
