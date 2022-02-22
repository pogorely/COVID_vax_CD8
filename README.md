## Data and scripts for study Minervina, Pogorelyy et al. 2021

This R code uses aggregated _cellranger_ output (folder _libs_aggregate_rev_) and adds TCR, UMAP GEX coordinates, GEX cluster, and epitope specificity info for each CD8+ T cell from the study. 

## How to run
Download and unzip _libs_aggregate_rev_ folder from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6232103.svg)](https://doi.org/10.5281/zenodo.6232103)

Install neccesary R packages (_Seurat_ version 4.0.4, _data.table_, _igraph_, _stringr_, _viridis_). Run:

    source("postprocessing.R")

This should take a few hours to run all of the analysis and will result in _cd8_only_dextr_rev_clean.tsv_ file (table with aggregated TCR, GEX and MHC-multimer specificity) in the working directory.

Please note, that different _Seurat_ versions might result in slightly different GEX clustering. To reproduce paper results exactly, install _Seurat_ v. 4.0.4