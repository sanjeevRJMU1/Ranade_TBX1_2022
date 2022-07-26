# script to generate the input files for monocle from Seurat cardiac cluster
library(data.table)

## 
Tbx1_CNC[["seurat_idents"]] <- Idents(Tbx1_CNC)

# subset just the e10.5 time point
E105_cnc <- SubsetData(Tbx1_CNC, subset.name = "gem.group",accept.value = c("WT_E105","KO_E105"))
table(E105_cnc@active.ident)
table(E105_cnc$gem.group)

#export raw counts
data_to_write_out <- as.matrix(E105_cnc@assays$SCT@counts)
#saveRDS(cardiac_matrix, file = "/Users/sranade/scRNA-seq/2019_scRNA-seq/iPSC/all_aggr_monocle/cardiac_matrix.RDS")

#generate cell netadata
cell_metadata <- as.data.frame(E105_cnc@meta.data)
#saveRDS(cell_metadata, file = "/Users/sranade/scRNA-seq/2019_scRNA-seq/iPSC/all_aggr_monocle/cell_metadata.RDS")

#generate gene_metadata 
gene_short_name <- row.names(data_to_write_out)
gene_metadata <- data.frame(gene_short_name)
row.names(gene_metadata) <- gene_short_name
#saveRDS(gene_metadata, file = "/Users/sranade/scRNA-seq/2019_scRNA-seq/iPSC/all_aggr_monocle/gene_metadata.RDS")

# start monocle
library(monocle3)
library(ggplot2)
library(dplyr)

#load input files
#cardiac_matrix <- readRDS("~/scRNA-seq/2019_scRNA-seq/iPSC/all_aggr_monocle/cardiac_matrix.RDS")
#cell_metadata <- readRDS("~/scRNA-seq/2019_scRNA-seq/iPSC/all_aggr_monocle/cell_metadata.RDS")
#gene_metadata <- readRDS("~/scRNA-seq/2019_scRNA-seq/iPSC/all_aggr_monocle/gene_metadata.RDS")

#cardiac_matrix <- as.matrix(cardiac_matrix)

# Make the CDS object
cds <- new_cell_data_set(data_to_write_out,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)
# pre-process data
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds, verbose = T, reduction_method = "UMAP", preprocess_method = "PCA")
plot_cells(cds)
plot_cells(cds, color_cells_by="seurat_idents")
plot_cells(cds, color_cells_by="gem.group")
plot_cells(cds, genes=c("Hoxa3"))

#check for and correct batch effect --> skipped bc I think this is over-fitting
# plot_cells(cds, color_cells_by="gem.group", label_cell_groups=FALSE)
cds = align_cds(cds, num_dim = 100, alignment_group = "gem.group")
cds = reduce_dimension(cds)
plot_cells(cds)
plot_cells(cds, color_cells_by="seurat_idents")
plot_cells(cds, color_cells_by="gem.group")

# Cluster cells
cds = cluster_cells(cds, reduction_method = "UMAP", resolution=1e-3)
plot_cells(cds)

plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

# # Find marker genes
# marker_test_res <- top_markers(cds, group_cells_by="partition", 
#                                reference_cells=1000, cores=8)
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
# #Dot plot of top marker genes
# top_specific_markers <- marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(3, pseudo_R2)
# 
# top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
# 
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="partition",
#                     ordering_type="cluster_row_col",
#                     max.size=3)
# at this point, I skipped over to single cell trajectories
# learn the trajectory graph
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=T,
           label_branch_points=FALSE)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

#DGE, skip the fit models section and go to graph test
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=4)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

# modules of co-expressed genes
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
### Run this later!
#cell_group_df <- tibble::tibble(cell=row.names(colData(neurons_cds)), 
#                                cell_group=partitions(cds)[colnames(neurons_cds)])
#agg_mat <- aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
#row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
#colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
#
#pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
#                   scale="column", clustering_method="ward.D2",
#                   fontsize=6)
#
#Finding genes that change as a function of pseudotime
plot_cells(cds,
           color_cells_by = "seurat_idents",
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=T)

avc_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids_prin_graph <- row.names(subset(avc_cds_pr_test_res, q_value < 0.05))

plot_cells(cds, genes=c("Dlx5","Hand2","Hoxa3","Dlx1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

avc_genes_5 = c("Hand2")
avc_cds_5 = cds[rowData(cds)$gene_short_name %in% avc_genes_5]

plot_genes_in_pseudotime(avc_cds_5,
                         color_cells_by="gem.group",
                         cell_size = 0.2,
                         min_expr = 0.5
)

