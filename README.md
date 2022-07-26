## Introduction

This is all code used by Ranade et. al. for the submission "Single Cell Epigenetics Reveal Cell-Cell Communication Networks in Normal and Abnormal Cardiac Morphogenesis" to [awaiting journal decision].

Manuscript is currently available on bioRxiv: [DOI 2022.07.25.501458](https://www.biorxiv.org/content/10.1101/2022.07.25.501458v1)

## Analysis
All data was processed and analyzed using Cellranger, Seurat, ArchR and supporting packages as detailed in provided scripts. See 10x Genomics documenation for [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) and [Cellranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac) usage.

Analysis order:

1. Process scRNA/scATAC 10x Genomics Cellranger v5.0.1 & Cellranger-atac v1.2.0 pipelines:
   - `cellranger count` & `cellranger-atac count`
   - `cellranger aggr`
2. Process and analyze scRNA seq data with Seurat v3 and v4 using scripts `scRNA-seq/*/` folder:
3. Process and analyze scATAC data with ArchR v1.0.1 using scripts in `scATAC-seq/*` folder

Please flag any issues and contact Dr. Sanjeev Ranade at `sanjeev.ranade@gladstone.ucsf.edu`.

## Data Availability

All sequencing data is available via GEO/SRA upon publication: [GSE198567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198567)