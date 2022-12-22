#!/usr/local/bin/Rscript
# Install independently to avoid installing errors
options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
BiocManager::install("Rsamtools", update = FALSE)
BiocManager::install("reactome.db", update = FALSE)
BiocManager::install("DelayedArray", update = FALSE)
BiocManager::install("BSgenome", update = FALSE)
pak::repo_add(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
pak::pkg_install("uebvhir/maUEB")
pak::pkg_install("clariomsmousetranscriptcluster.db")
