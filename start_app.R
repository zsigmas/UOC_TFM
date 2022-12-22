if (!requireNamespace("brightPlots")) {
    devtools::install_github("zsigmas/brightPlots")
}

print(brightPlots::app(targets::tar_read(toptable_BH)))
