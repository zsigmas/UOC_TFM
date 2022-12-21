#!/usr/local/bin/Rscript
options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
BiocManager::install("Rsamtools", update = FALSE)
BiocManager::install("reactome.db", update = FALSE)
BiocManager::install("DelayedArray", update = FALSE)
BiocManager::install("BSgenome", update = FALSE)
