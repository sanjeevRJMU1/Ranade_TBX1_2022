### so you've added your motif deviation score metadata and have run DAR analysis on high v low deviations
### do this next
### requires deviationThreshold metadta in ArchR project!




### Pull scRNAseq cell barcodes associated with high or low deviation scores and run DGE...
## load seurat object
YOURseuratObj <- readRDS("./data/_source_RDS_objects/NC_labeled_07-17-2021.RDS")

### create Seurat metadata to add
## pull unique scRNA cell barcodes
highDeviationCells <- unique(E105_DM_Harmony@cellColData@listData$predictedCell_NC[E105_DM_Harmony$deviationThreshold == "above_6e-2_Pitx1_434and452"])
lowDeviationCells <- unique(E105_DM_Harmony@cellColData@listData$predictedCell_NC[E105_DM_Harmony$deviationThreshold == "below_6e-2_Pitx1_434and452"])
## exclude cells from low that are in high
# forgot to use "negated" %in% operator when we were on the phone. PEBCAK.
`%!in%` <- Negate(`%in%`)
lowDeviationCells <- lowDeviationCells[lowDeviationCells %!in% highDeviationCells]

## pull all cells from Seurat obj
metadata2add <- data.frame(scRNA_barcode = colnames(YOURseuratObj), deviation_value = "not_mapped_to_ATAC")
## annotate deviations
metadata2add$deviation_value[metadata2add$scRNA_barcode %in% highDeviationCells] <- "high"
metadata2add$deviation_value[metadata2add$scRNA_barcode %in% lowDeviationCells] <- "low"
## set colnames, drop barcode col
rownames(metadata2add) <- metadata2add$scRNA_barcode
metadata2add$scRNA_barcode <- NULL

### add metadata
YOURseuratObj <- AddMetaData(YOURseuratObj, metadata = metadata2add)

### dimplot
DimPlot(YOURseuratObj, group.by = "deviation_value")

Idents(YOURseuratObj) <- "deviation_value"

### DE
## sending it
highVlow_deviationCells <- FindMarkers(YOURseuratObj, ident.1 = "high", ident.2 = "low")
## boyscout edition
highVlow_boyscoutEdition <- FindMarkers(YOURseuratObj, ident.1 = "high", ident.2 = "low", max.cells.per.ident = min(table(Idents(YOURseuratObj))))

