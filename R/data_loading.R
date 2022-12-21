get_cel_files <- function(targets, cel_dir) {
    checkmate::assert_character(targets$FileName, any.missing = FALSE)
    file.path(cel_dir, targets$FileName)
}

read_targets_alias <- function(target_file, target_fact) {
    target_dir <- dirname(target_file)
    target_file <- basename(target_file)
    maUEB::read_targets(inputDir = target_dir,
        targetsFN = target_file,
        targets.fact = target_fact)
}

read_celfiles_alias  <- function(cel_files, targets) {
    cel_dir <- unique(dirname(cel_files))
    maUEB::read_celfiles(inputDir = cel_dir, targets = targets)
}