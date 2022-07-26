## Introduction

This is all code used by Ranade et. al. for the submission "Single Cell Epigenetics Reveal Cell-Cell Communication Networks in Normal and Abnormal Cardiac Morphogenesis" to [awaiting journal decision].

Manuscript is currently available on bioRxiv: [DOI 2022.07.25.501458](https://www.biorxiv.org/content/10.1101/2022.07.25.501458v1)

## Analysis
All data was processed and analyzed using Cellranger, Seurat, ArchR and supporting packages as detailed in provided scripts. See 10x Genomics documenation for [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) and [Cellranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac) usage.

Analysis order:

1. Process scRNA/scATAC 10x Genomics Cellranger v5.x.x & Cellranger-atac v2.x.x pipelines:
   - `cellranger count` & `cellranger-atac count`
   - `cellranger aggr` & `cellranger-atac aggr`
2. Process and analyze scRNA seq data with Seurat v4.x.x using scripts `scRNA-seq/*/` folder:
3. Process and analyze scATAC data with ArchR v1.x.x using scripts in `scATAC-seq/*` folder

Please flag any issues and contact Dr. Sanjeev Ranade at `sanjeev.ranade@gladstone.ucsf.edu`.

## Data Availability

All sequencing data will be available via GEO/SRA upon publication: [link-provided-when-data-released](https://www.ncbi.nlm.nih.gov/sra)