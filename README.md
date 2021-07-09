## Data and scripts for study Minervina, Pogorelyy et al. 2021

This R code uses aggregated _cellranger_ output (folder _libs_aggregate_) and adds TCR, UMAP GEX coordinates, GEX cluster, and epitope specificity info for each CD8+ T cell from the study. 

## How to run
Install neccesary R packages (_Seurat_ version 3.2.3, _data.table_, _igraph_, _stringr_, _viridis_). Run:

    source("postprocessing.R")

This should take about an hour to run all of the analysis and will result in _cd8_only_dextr.tsv_ file (table with aggregated TCR, GEX and MHC-multimer specificity) in the working directory.

Please note, that different _Seurat_ versions might result in slightly different GEX clustering. To reproduce paper results exactly, install _Seurat_ v. 3.2.3

    remotes::install_version("spatstat", version = "1.64-1")

Depending on your setup you might need to downgrade _spatstat_ dependency
first:

    remotes::install_version("spatstat", version = "1.64-1")
    
