# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Parameters

get_parameters <- function() {
  data_dir <- "data"
  cel_dir <- "celfiles"
  cel_files <- list.files("celfiles", full.names = TRUE)
  targets_file_name <- "targets.RSRCHR.STUDY.csv"
  targets_file <- file.path(data_dir, targets_file_name)
  results_dir <- "results"
  targets_fact <- c("Group", "Batch")
  colorlist <- c(
    "firebrick1", "blue", "darkgoldenrod1",
    "darkorchid1", "lightblue", "orange4",
    "seagreen", "aquamarine3", "red4",
    "chartreuse2", "plum", "darkcyan",
    "darkolivegreen3", "forestgreen", "gold",
    "khaki1", "lightgreen", "gray67",
    "deeppink"
    ) # list of colors for some plots, alternative to "Colors" column  stopifnot(file.exists(targets_file))
  summary_fn <- "ResultsSummary_Load-QC-Norm-Filt.txt"

  boxplot <- list(
    raw = list(
      group = "Group",
      group_color = "Colors",
      sample_names = "ShortName",
      output_dir = results_dir,
      label = "RawData"
    ),
    norm = list(
      group = "Group",
      group_color = "Colors",
      sample_names = "ShortName",
      output_dir = results_dir,
      label = "NormData"
    )
  )

  pca <- list(
    raw = list(
      fact = targets_fact,
      scale = FALSE,
      col_group = "Colors",
      color_list = colorlist,
      output_dir = results_dir,
      label = "RawData"
    ),
    norm = list(
      fact = targets_fact,
      scale = FALSE,
      col_group = "Colors",
      label = "NormData",
      color_list = colorlist,
      output_dir = results_dir,
      batch_remove = TRUE,
      batch_factors = "Batch",
      size = 1.5,
      glineas = 0.25
    )
  )

  batch_pca <- list(
    norm = c(
      pca$norm,
      list(
        batch_remove = TRUE,
        batch_factors = "Batch",
        size = 1.5,
        glineas = 0.25
      )
    )
  )

  pvca <- list(
    norm = list(
      pct_threshold = .6,
      targets_file_name = targets_file_name,
      fact = targets_fact,
      label = "NormData",
      summary_fn = summary_fn
    )
  )

  hc <- list(
    raw = list(
      method = "average",
      cex_row = 0.6,
      cex_col = 0.6,
      rect = TRUE,
      num_clusters = 2,
      output_dir = results_dir,
      label = "RawData"
    ),
    norm = list(
      method = "average",
      cex_row = 0.6,
      cex_col = 0.6,
      rect = TRUE,
      num_clusters = 2,
      output_dir = results_dir,
      label = "NormData"
    )
  )

  aqm <- list(
    raw = list(
      out_dir = file.path(results_dir, "QCDir.Raw"),
      force = TRUE,
      int_group = targets_fact
    ),
    norm = list(
      out_dir = file.path(results_dir, "QCDir.Norm"),
      force = TRUE,
      int_group = targets_fact
    )
  )

  norm <- list(
    annot_pkg = "clariomsmousetranscriptcluster.db",
    method = "rma",
    output_fn = "Normalized.All",
    annot = list(
      output_fn = "Annotations.AllGenes",
      saveHTML = TRUE,
      title = "Annotations for all genes"
    )
  )

  outlier_samples <- "CMP.2"

  list(
    targets_file = targets_file,
    cel_dir = cel_dir,
    targets_fact = targets_fact,
    results_dir = results_dir,
    color_list = colorlist,
    boxplot = boxplot,
    pca = pca,
    batch_pca = batch_pca,
    pvca = pvca,
    hc = hc,
    aqm = aqm,
    norm = norm,
    outlier_samples = outlier_samples
  )
}

par <- get_parameters()

# require(par$norm$annot_pkg, character.only = TRUE)
# require(annaffy)


# Set target options:
tar_option_set(
  packages = c("tibble", "oligoClasses", "annaffy", par$norm$annot_pkg), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Parameters

get_parameters <- function() {
  data_dir <- "data"
  cel_dir <- "celfiles"
  cel_files <- list.files("celfiles", full.names = TRUE)
  targets_file_name <- "targets.RSRCHR.STUDY.csv"
  targets_file <- file.path(data_dir, targets_file_name)
  results_dir <- "results"
  targets_fact <- c("Group", "Batch")
  colorlist <- c(
    "firebrick1", "blue", "darkgoldenrod1",
    "darkorchid1", "lightblue", "orange4",
    "seagreen", "aquamarine3", "red4",
    "chartreuse2", "plum", "darkcyan",
    "darkolivegreen3", "forestgreen", "gold",
    "khaki1", "lightgreen", "gray67",
    "deeppink"
    ) # list of colors for some plots, alternative to "Colors" column  stopifnot(file.exists(targets_file))
  summary_fn <- "ResultsSummary_Load-QC-Norm-Filt.txt"

  boxplot <- list(
    raw = list(
      group = "Group",
      group_color = "Colors",
      sample_names = "ShortName",
      output_dir = results_dir,
      label = "RawData"
    ),
    norm = list(
      group = "Group",
      group_color = "Colors",
      sample_names = "ShortName",
      output_dir = results_dir,
      label = "NormData"
    ),
    norm_f = list(
      group = "Group",
      group_color = "Colors",
      sample_names = "ShortName",
      output_dir = results_dir,
      label = "NormFData"
    )
  )

  pca <- list(
    raw = list(
      fact = targets_fact,
      scale = FALSE,
      col_group = "Colors",
      color_list = colorlist,
      output_dir = results_dir,
      label = "RawData"
    ),
    norm = list(
      fact = targets_fact,
      scale = FALSE,
      col_group = "Colors",
      label = "NormData",
      color_list = colorlist,
      output_dir = results_dir,
      batch_remove = TRUE,
      batch_factors = "Batch",
      size = 1.5,
      glineas = 0.25
    ),
    norm_f = list(
      fact = targets_fact,
      scale = FALSE,
      col_group = "Colors",
      label = "NormFData",
      color_list = colorlist,
      output_dir = results_dir,
      batch_remove = TRUE,
      batch_factors = "Batch",
      size = 1.5,
      glineas = 0.25
    )
  )

  batch_pca <- list(
    norm = c(
      pca$norm,
      list(
        batch_remove = TRUE,
        batch_factors = "Batch",
        size = 1.5,
        glineas = 0.25
      )
    ),
    norm_f = c(
      pca$norm_f,
      list(
        batch_remove = TRUE,
        batch_factors = "Batch",
        size = 1.5,
        glineas = 0.25
      )
    )
  )

  pvca <- list(
    norm = list(
      pct_threshold = .6,
      targets_file_name = targets_file_name,
      fact = targets_fact,
      label = "NormData",
      summary_fn = summary_fn
    ),
    norm_f = list(
      pct_threshold = .6,
      targets_file_name = targets_file_name,
      fact = targets_fact,
      label = "NormFData",
      summary_fn = summary_fn
    )
  )

  hc <- list(
    raw = list(
      method = "average",
      cex_row = 0.6,
      cex_col = 0.6,
      rect = TRUE,
      num_clusters = 2,
      output_dir = results_dir,
      label = "RawData"
    ),
    norm = list(
      method = "average",
      cex_row = 0.6,
      cex_col = 0.6,
      rect = TRUE,
      num_clusters = 2,
      output_dir = results_dir,
      label = "NormData"
    ),
    norm_f = list(
      method = "average",
      cex_row = 0.6,
      cex_col = 0.6,
      rect = TRUE,
      num_clusters = 2,
      output_dir = results_dir,
      label = "NormFData"
    )
  )

  aqm <- list(
    raw = list(
      out_dir = file.path(results_dir, "QCDir.Raw"),
      force = TRUE,
      int_group = targets_fact
    ),
    norm = list(
      out_dir = file.path(results_dir, "QCDir.Norm"),
      force = TRUE,
      int_group = targets_fact
    ),
    norm_f = list(
      out_dir = file.path(results_dir, "QCDir.NormF"),
      force = TRUE,
      int_group = targets_fact
    )
  )

  norm <- list(
    annot_pkg = "clariomsmousetranscriptcluster.db",
    method = "rma",
    output_fn = "Normalized.All",
    annot = list(
      output_fn = "Annotations.AllGenes",
      saveHTML = TRUE,
      title = "Annotations for all genes"
    )
  )

  norm_f <- list(
    annot_pkg = "clariomsmousetranscriptcluster.db",
    method = "rma",
    output_fn = "NormalizedF.All",
    annot = list(
      output_fn = "AnnotationsF.AllGenes",
      saveHTML = TRUE,
      title = "Annotations for all genes"
    )
  )

  outlier_samples <- "CMP.2"

  list(
    targets_file = targets_file,
    cel_dir = cel_dir,
    targets_fact = targets_fact,
    results_dir = results_dir,
    color_list = colorlist,
    boxplot = boxplot,
    pca = pca,
    batch_pca = batch_pca,
    pvca = pvca,
    hc = hc,
    aqm = aqm,
    norm = norm,
    norm_f = norm_f,
    outlier_samples = outlier_samples
  )
}

par <- get_parameters()

require(par$norm$annot_pkg, character.only = TRUE)
require(annaffy)

# Resolvemos antes y así no recalculamos todo el rato.

# Extraer como targets los paths a los archivoes que están almacenados como atributos

evaler <- function(x) {
  rlang::eval_tidy(rlang::enexpr(x))
}

# Replace the target list below with your own:
list(
  # FILES ----
  tar_target(
    name = targets_file,
    command = !!par$targets_file,
    format = "file"
  ),
  tar_target(
    name = targets,
    command = read_targets_alias(
      targets_file,
      !!par$targets_fact
      )
  ),
  tar_target(
    name = cel_files,
    command = get_cel_files(
      targets,
      !!par$cel_dir
      ),
    format = "file"
  ),
  # RAW ----
  tar_target(
    name = eset_raw,
    command = read_celfiles_alias(
      cel_files,
      targets
      )
  ),
  tar_target(
    name = qc_raw_boxplot_file,
    command = qc_boxplot_alias(
      data = eset_raw,
      group = !!par$boxplot$raw$group,
      group.color = !!par$boxplot$raw$group_color,
      samplenames = !!par$boxplot$raw$sample_names,
      outputDir = !!par$pca$raw$output_dir,
      label = !!par$boxplot$raw$label
      ),
    format = "file"
  ),
  tar_target(
    name = qc_raw_load_pca1,
    command = qc_pca1_alias(
      data = Biobase::exprs(eset_raw),
      scale = !!par$pca$raw$scale,
      pca2D_factors = !!par$pca$raw$fact,
      targets = targets,
      col.group = !!par$pca$raw$col_group,
      colorlist = !!par$pca$raw$color_list,
      names = targets$ShortName,
      outputDir = !!par$pca$raw$output_dir,
      label = !!par$pca$raw$group_color
      )
  ),
  tar_target(
    name = qc_raw_hc_file,
    command = qc_hc_alias(
      data = Biobase::exprs(eset_raw),
      hclust.method = !!par$hc$raw$method,
      names = targets$ShortName,
      cexRow = !!par$hc$raw$cex_row,
      cexCol = !!par$hc$raw$cex_col,
      rect = !!par$hc$raw$rect,
      numclusters = !!par$hc$raw$num_clusters,
      outputDir = !!par$hc$raw$output_dir,
      label = !!par$hc$raw$label
      ),
    format = "file"
  ),
  tar_target(
    name = qc_raw_arrayQM,
    command = arrayQualityMetrics_alias(
      eset_raw,
      outdir = !!par$aqm$raw$out_dir,
      force = !!par$aqm$raw$force,
      intgroup = !!par$aqm$raw$int_group
      )
  ),
  # NORM ----
  tar_target(
    name = eset_norm,
    command = normalization_alias(
      data = eset_raw,
      norm.method = !!par$norm$method,
      annotPkg = !!par$norm$annot_pkg,
      outputFN = !!par$norm$output_fn,
      outputDir = !!par$results_dir
      )
  ),
  tar_target(
    name = norm_annot,
    command = save_annotations_alias(
      data = rownames(eset_norm),
      annotPkg = !!par$norm$annot_pkg,
      outputFN = !!par$norm$annot$output_fn,
      saveHTML = !!par$norm$annot$saveHTML,
      title = !!par$norm$annot$title,
      outputDir = !!par$results_dir
      )
  ),
  tar_target(
    name = qc_norm_boxplot_file,
    command = qc_boxplot_alias(
      data = eset_norm,
      group = !!par$boxplot$norm$group,
      group.color = !!par$boxplot$norm$group_color,
      samplenames = !!par$boxplot$norm$sample_names,
      outputDir = !!par$pca$norm$output_dir,
      label = !!par$boxplot$norm$label
      ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_load_pca1,
    command = qc_pca1_alias(
      data = Biobase::exprs(eset_norm),
      scale = !!par$pca$norm$scale,
      pca2D_factors = !!par$pca$norm$fact,
      targets = targets,
      col.group = !!par$pca$norm$col_group,
      colorlist = !!par$pca$norm$color_list,
      names = targets$ShortName,
      outputDir = !!par$pca$norm$output_dir,
      label = !!par$pca$norm$group_color
      )
  ),
  tar_target(
    name = qc_norm_corr_batch_load_pca1,
    command = qc_pca1_alias(
      data = Biobase::exprs(eset_norm),
      scale = !!par$batch_pca$norm$scale,
      pca2D_factors = !!par$batch_pca$norm$fact,
      targets = targets,
      col.group = !!par$batch_pca$norm$col_group,
      colorlist = !!par$batch_pca$norm$color_list,
      names = targets$ShortName,
      outputDir = !!par$batch_pca$norm$output_dir,
      label = !!par$pca$norm$group_color,
      batchRemove = !!par$batch_pca$norm$batch_remove,
      batchFactors = !!par$batch_pca$norm$batch_factors,
      size = !!par$batch_pca$norm$size,
      glineas = !!par$batch_pca$norm$glineas
      )
  ),
  tar_target(
    name = qc_norm_hc_file,
    command = qc_hc_alias(
      data = Biobase::exprs(eset_norm),
      hclust.method = !!par$hc$norm$method,
      names = targets$ShortName,
      cexRow = !!par$hc$norm$cex_row,
      cexCol = !!par$hc$norm$cex_col,
      rect = !!par$hc$norm$rect,
      numclusters = !!par$hc$norm$num_clusters,
      outputDir = !!par$hc$norm$output_dir,
      label = !!par$hc$norm$label
      ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_arrayQM,
    command = arrayQualityMetrics_alias(
      eset_norm,
      outdir = !!par$aqm$norm$out_dir,
      force = !!par$aqm$norm$force,
      intgroup = !!par$aqm$norm$int_group
      )
  ),
  # BLOCKED UNTIL https://github.com/uebvhir/maUEB/issues/3 IS SOLVED
  # tar_target(
  #   name = qc_norm_pcpva,
  #   command = evaler(maUEB::qc_pvca(data = eset_norm, factors = !!par$pvca$norm$fact, targets = !!par$pvca$norm$targets_file_name, pct_threshold = !!par$pvca$norm$threshold, label = !!par$pvca$norm$label, outputDir = !!par$pvca$norm$output_dir, summaryFN = !!par$pvca$norm$summary_fn))
  # ),
  # SAMPLE REMOVAL ----
  tar_target(
    name = eset_raw_f,
    command = remove_samples(
      eset_raw,
      par$outlier_samples
      )
  ),
  tar_target(
    name = targets_f,
    command = pData(
      eset_raw_f
      )
  ),
  # SAMPLE REMOVAL NORM ----
  tar_target(
    name = eset_norm_f,
    command = normalization_alias(
      data = eset_raw_f,
      norm.method = !!par$norm_f$method,
      annotPkg = !!par$norm_f$annot_pkg,
      outputFN = !!par$norm_f$output_fn,
      outputDir = !!par$results_dir
      )
  ),
  # tar_target( #MAYBE REMOVED
  #   name = norm_annot,
  #   command = save_annotations_alias(
  #     data = rownames(eset_norm),
  #     annotPkg = !!par$norm$annot_pkg,
  #     outputFN = !!par$norm$annot$output_fn,
  #     saveHTML = !!par$norm$annot$saveHTML,
  #     title = !!par$norm$annot$title,
  #     outputDir = !!par$results_dir
  #     )
  # ),
  tar_target(
    name = qc_norm_boxplot_file_f,
    command = qc_boxplot_alias(
      data = eset_norm_f,
      group = !!par$boxplot$norm_f$group,
      group.color = !!par$boxplot$norm_f$group_color,
      samplenames = !!par$boxplot$norm_f$sample_names,
      outputDir = !!par$pca$norm_f$output_dir,
      label = !!par$boxplot$norm_f$label
      ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_load_pca1_f,
    command = qc_pca1_alias(
      data = Biobase::exprs(eset_norm_f),
      scale = !!par$pca$norm_f$scale,
      pca2D_factors = !!par$pca$norm_f$fact,
      targets = targets_f,
      col.group = !!par$pca$norm_f$col_group,
      colorlist = !!par$pca$norm_f$color_list,
      names = targets$ShortName,
      outputDir = !!par$pca$norm_f$output_dir,
      label = !!par$pca$norm_f$group_color
      )
  ),
  tar_target(
    name = qc_norm_corr_batch_load_pca1_f,
    command = qc_pca1_alias(
      data = Biobase::exprs(eset_norm_f),
      scale = !!par$batch_pca$norm_f$scale,
      pca2D_factors = !!par$batch_pca$norm_f$fact,
      targets = targets_f,
      col.group = !!par$batch_pca$nor_fm$col_group,
      colorlist = !!par$batch_pca$norm_f$color_list,
      names = targets_f$ShortName,
      outputDir = !!par$batch_pca$norm_f$output_dir,
      label = !!par$pca$norm_f$group_color,
      batchRemove = !!par$batch_pca$norm_f$batch_remove,
      batchFactors = !!par$batch_pca$norm_f$batch_factors,
      size = !!par$batch_pca$norm_f$size,
      glineas = !!par$batch_pca$norm_f$glineas
      )
  ),
  tar_target(
    name = qc_norm_hc_file_f,
    command = qc_hc_alias(
      data = Biobase::exprs(eset_norm_f),
      hclust.method = !!par$hc$norm_f$method,
      names = targets_f$ShortName,
      cexRow = !!par$hc$norm_f$cex_row,
      cexCol = !!par$hc$norm_f$cex_col,
      rect = !!par$hc$norm_f$rect,
      numclusters = !!par$hc$norm_f$num_clusters,
      outputDir = !!par$hc$norm_f$output_dir,
      label = !!par$hc$norm_f$label
      ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_arrayQM_f,
    command = arrayQualityMetrics_alias(
      eset_norm_f,
      outdir = !!par$aqm$norm_f$out_dir,
      force = !!par$aqm$norm_f$force,
      intgroup = !!par$aqm$norm_f$int_group
      )
  )

)

## CONTINUE AT SAMPLE EXCLUSION

