# Analysis of single-cell retina data

This repository contains some code for the analysis of single-cell data of human and mouse retina

## Files description

- `annotate_cell_types.R`: R script for adding manual cell type annotation to an already saved [Seurat](https://satijalab.org/seurat/) object

- `cluster-celltypes`: TSV file containing manually annotated clusters for human data

- `cluster-celltypes-mouse`: TSV file containing manually annotated clusters for mouse data

- `human_preprocessing.R`: R script for preprocessing human retina single-cell data according to a publication

- `mouse_preprocessing.R`: R script for preprocessing mouse retina single-cell data according to a publication

- `norm_counts_n_genes.R`: R script for generating normalized tables of gene counts per sample and number of cells expressing each gene per sample