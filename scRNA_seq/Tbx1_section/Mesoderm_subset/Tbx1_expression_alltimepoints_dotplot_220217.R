### subset individual time points, wt only
wt_subset_e95 <- SubsetData(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E925"))
wt_subset_e105 <- SubsetData(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E105"))
wt_subset_e115 <- SubsetData(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E115"))

wt_cells <- OldWhichCells(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E925","WT_E105","WT_E115"))
wt_subset_tbx1_expr <- SubsetData(Tbx1_mesoderm_cnc, cells =wt_cells)
table(wt_subset_tbx1_expr@active.ident)

## get cell barcodes for each ident per time point
e95_Posterior_Second_Heart_Field_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Posterior_Second_Heart_Field")
e105_Posterior_Second_Heart_Field_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Posterior_Second_Heart_Field")
e115_Posterior_Second_Heart_Field_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Posterior_Second_Heart_Field")

e95_Paraxial_Mesoderm_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Paraxial_Mesoderm")
e105_Paraxial_Mesoderm_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Paraxial_Mesoderm")
e115_Paraxial_Mesoderm_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Paraxial_Mesoderm")

e95_Cardiopharyngeal_Mesoderm_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Cardiopharyngeal_Mesoderm")
e105_Cardiopharyngeal_Mesoderm_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Cardiopharyngeal_Mesoderm")
e115_Cardiopharyngeal_Mesoderm_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Cardiopharyngeal_Mesoderm")

e95_Cardiopharyngeal_Mesenchyme_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Cardiopharyngeal_Mesenchyme")
e105_Cardiopharyngeal_Mesenchyme_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Cardiopharyngeal_Mesenchyme")
e115_Cardiopharyngeal_Mesenchyme_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Cardiopharyngeal_Mesenchyme")

e95_Neural_Crest_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Neural_Crest")
e105_Neural_Crest_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Neural_Crest")
e115_Neural_Crest_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Neural_Crest")

e95_Cranial_Mesenchyme_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Cranial_Mesenchyme")
e105_Cranial_Mesenchyme_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Cranial_Mesenchyme")
e115_Cranial_Mesenchyme_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Cranial_Mesenchyme")

e95_Smooth_Muscle_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Smooth_Muscle")
e105_Smooth_Muscle_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Smooth_Muscle")
e115_Smooth_Muscle_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Smooth_Muscle")

e95_Epicardium_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Epicardium")
e105_Epicardium_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Epicardium")
e115_Epicardium_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Epicardium")

e95_Anterior_Second_Heart_Field_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Anterior_Second_Heart_Field")
e105_Anterior_Second_Heart_Field_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Anterior_Second_Heart_Field")
e115_Anterior_Second_Heart_Field_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Anterior_Second_Heart_Field")

e95_Cardiomyocyte_cells <- OldWhichCells(wt_subset_e95, ident.keep = "Cardiomyocyte")
e105_Cardiomyocyte_cells <- OldWhichCells(wt_subset_e105, ident.keep = "Cardiomyocyte")
e115_Cardiomyocyte_cells <- OldWhichCells(wt_subset_e115, ident.keep = "Cardiomyocyte")


###
Idents(object = wt_subset_tbx1_expr, cells= e95_Posterior_Second_Heart_Field_cells) <- "e95_Posterior_Second_Heart_Field"
Idents(object = wt_subset_tbx1_expr, cells= e105_Posterior_Second_Heart_Field_cells) <- "e105_Posterior_Second_Heart_Field"
Idents(object = wt_subset_tbx1_expr, cells= e115_Posterior_Second_Heart_Field_cells) <- "e115_Posterior_Second_Heart_Field"

Idents(object = wt_subset_tbx1_expr, cells= e95_Paraxial_Mesoderm_cells) <- "e95_Paraxial_Mesoderm"
Idents(object = wt_subset_tbx1_expr, cells= e105_Paraxial_Mesoderm_cells) <- "e105_Paraxial_Mesoderm"
Idents(object = wt_subset_tbx1_expr, cells= e115_Paraxial_Mesoderm_cells) <- "e115_Paraxial_Mesoderm"

Idents(object = wt_subset_tbx1_expr, cells= e95_Cardiopharyngeal_Mesoderm_cells) <- "e95_Cardiopharyngeal_Mesoderm"
Idents(object = wt_subset_tbx1_expr, cells= e105_Cardiopharyngeal_Mesoderm_cells) <- "e105_Cardiopharyngeal_Mesoderm"
Idents(object = wt_subset_tbx1_expr, cells= e115_Cardiopharyngeal_Mesoderm_cells) <- "e115_Cardiopharyngeal_Mesoderm"

Idents(object = wt_subset_tbx1_expr, cells= e95_Cardiopharyngeal_Mesenchyme_cells) <- "e95_Cardiopharyngeal_Mesenchyme"
Idents(object = wt_subset_tbx1_expr, cells= e105_Cardiopharyngeal_Mesenchyme_cells) <- "e105_Cardiopharyngeal_Mesenchyme"
Idents(object = wt_subset_tbx1_expr, cells= e115_Cardiopharyngeal_Mesenchyme_cells) <- "e115_Cardiopharyngeal_Mesenchyme"

Idents(object = wt_subset_tbx1_expr, cells= e95_Neural_Crest_cells) <- "e95_Neural_Crest"
Idents(object = wt_subset_tbx1_expr, cells= e105_Neural_Crest_cells) <- "e105_Neural_Crest"
Idents(object = wt_subset_tbx1_expr, cells= e115_Neural_Crest_cells) <- "e115_Neural_Crest"

Idents(object = wt_subset_tbx1_expr, cells= e95_Cranial_Mesenchyme_cells) <- "e95_Cranial_Mesenchyme"
Idents(object = wt_subset_tbx1_expr, cells= e105_Cranial_Mesenchyme_cells) <- "e105_Cranial_Mesenchyme"
Idents(object = wt_subset_tbx1_expr, cells= e115_Cranial_Mesenchyme_cells) <- "e115_Cranial_Mesenchyme"

Idents(object = wt_subset_tbx1_expr, cells= e95_Smooth_Muscle_cells) <- "e95_Smooth_Muscle"
Idents(object = wt_subset_tbx1_expr, cells= e105_Smooth_Muscle_cells) <- "e105_Smooth_Muscle"
Idents(object = wt_subset_tbx1_expr, cells= e115_Smooth_Muscle_cells) <- "e115_Smooth_Muscle"

Idents(object = wt_subset_tbx1_expr, cells= e95_Epicardium_cells) <- "e95_Epicardium"
Idents(object = wt_subset_tbx1_expr, cells= e105_Epicardium_cells) <- "e105_Epicardium"
Idents(object = wt_subset_tbx1_expr, cells= e115_Epicardium_cells) <- "e115_Epicardium"

Idents(object = wt_subset_tbx1_expr, cells= e95_Anterior_Second_Heart_Field_cells) <- "e95_Anterior_Second_Heart_Field"
Idents(object = wt_subset_tbx1_expr, cells= e105_Anterior_Second_Heart_Field_cells) <- "e105_Anterior_Second_Heart_Field"
Idents(object = wt_subset_tbx1_expr, cells= e115_Anterior_Second_Heart_Field_cells) <- "e115_Anterior_Second_Heart_Field"

Idents(object = wt_subset_tbx1_expr, cells= e95_Cardiomyocyte_cells) <- "e95_Cardiomyocyte"
Idents(object = wt_subset_tbx1_expr, cells= e105_Cardiomyocyte_cells) <- "e105_Cardiomyocyte"
Idents(object = wt_subset_tbx1_expr, cells= e115_Cardiomyocyte_cells) <- "e115_Cardiomyocyte"

table(wt_subset_tbx1_expr@active.ident)

DotPlot(wt_subset_tbx1_expr, features = c("Tbx1"), split.by = )



wt_cells <- OldWhichCells(Tbx1_mesoderm_cnc, subset.name = "gem.group", accept.value = c("WT_E925","WT_E105","WT_E115"))
wt_subset_tbx1_expr <- SubsetData(Tbx1_mesoderm_cnc, cells =wt_cells)
table(wt_subset_tbx1_expr@active.ident)
table(wt_subset_tbx1_expr$gem.group)

e95_Posterior_Second_Heart_Field_cells <- OldWhichCells(wt_subset_tbx1_expr, ide)