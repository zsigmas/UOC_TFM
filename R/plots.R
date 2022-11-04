file_path_must_work <- function(...) {
    fn <- file.path(...)
    stopifnot(file.exists(fn))
    fn
}

qc_boxplot_alias <- function(...) {
    args <- list(...)
    maUEB::qc_boxplot(...)
    file_path_must_work(args$outputDir, paste0("Boxplot", args$label, ".pdf"))
}

qc_pca1_alias <- function(...) {
    args <- list(...)
    structure(
        maUEB::qc_pca1(...),
        files = file_path_must_work(args$outputDir, paste0("PCA", args$label, ".pdf"))
    )
}

qc_hc_alias <- function(...) {
    args <- list(...)
    maUEB::qc_hc(...)
    file_path_must_work(args$outputDir, paste0("QC", args$label, ".pdf"))
}

arrayQualityMetrics_alias <- function(...) {
    args <- list(...)
    obj <- arrayQualityMetrics::arrayQualityMetrics(...)
    stopifnot(dir.exists(args$outdir))
    structure(
        obj,
        files = list.files(args$outdir, full.names = TRUE)
    )
}

normalization_alias <- function(...) {
    args <- list(...)
    structure(
        maUEB::normalization(...),
        files = file_path_must_work(args$outputDir, paste0(args$outputFN, ".csv"))
    )
}

save_annotations_alias <- function(...) {
    args <- list(...)
    obj <- maUEB::save_annotations(...)
    structure(
        obj,
        files = c(
            if (args$saveHTML) file_path_must_work(args$outputDir, paste0(args$outputFN, ".html")) else NULL,
            file_path_must_work(args$outputDir, paste0(args$outputFN, ".csv"))
        )
    )
}

remove_samples <- function(dataset, outlier_samples) {
    dataset[,-which(colnames(dataset) %in% outlier_samples)]
}
