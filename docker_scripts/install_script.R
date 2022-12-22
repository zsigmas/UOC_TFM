#!/usr/local/bin/Rscript
install.packages("pak")
pak::pak("targets")
pak::pak("tarchetypes")

# Otherwise gsea filt throws an error
pak::pak("ggnewscale")
