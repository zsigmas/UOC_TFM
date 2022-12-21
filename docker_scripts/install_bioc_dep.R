#!/usr/local/bin/Rscript
# Install independently to avoid installing errors
# BiocManager::install("Rsamtools", update = FALSE)
# BiocManager::install("reactome.db", update = FALSE)
options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
pak::repo_add(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")
pak::pkg_install("uebvhir/maUEB")
