
as_list_of_files <- function(fn, outputDir_args) {
    function(...) {
        fn(...)
        args <- list(...)
        unlist(purrr::map(args[outputDir_args], ~ list.files(.x, full.names = TRUE)))
    }
}

append_args_as_attr <- function(fn, appendable_args = NULL) {
    function(...) {
        args <- list(...)
        rlang::exec(structure, .Data = fn(...), !!!args[appendable_args])
    }
}

remove_samples <- function(dataset, outlier_samples) {
    dataset[, -which(colnames(dataset) %in% outlier_samples)]
}

abs_gsea_ReactomePAunfilt_alias <- function(...) {
    args <- list(...)
    annotPackage <<- args$annotPackage # Trick to place in global environment not recommended to push this forward. Fix package
    args <- args[names(args) != "annotPackage"]
    do.call(maUEB::abs_gsea_ReactomePAunfilt, args)
}
