# Load packages required to define the pipeline:
library(targets)
# devtools::install_local("../maUEB", force = TRUE)
# devtools::load_all("../maUEB")
requireNamespace("purrr")
requireNamespace("rlang")

# Helper function to create result directories on the fly
create_dir <- function(...) {
  path <- file.path(...)
  if (dir.exists(path)) {
    return(path)
  }
  success <- dir.create(path, showWarnings = TRUE, recursive = TRUE)
  if (success) {
    return(path)
  } else {
    stop("Error creating: ", path)
  }
}

# Possible improvements:
# Script level parameters for some directories

# Features:
# Parameters can be specified per stage
# A change of one parameter does not relaunch the whole script

# Cons:
# High coupling between data and presenatation in maUEB library that makes not the ideal candidate
# to use with targets

# Comments:
# The idea is that the system is pretty flexible the only requisit is that in the end each target
# receives a list with the arguments for the function call. How these parameters are stored can be
# decided by the user, external files, list structure such as the one below, etc. The division in blocks
# just follows the provided Rmd but it is far from compulsory and any other division can be made, from
# a single parameter structure per call to a list to covers the whole pipeline.


# Quality parameters ----

get_load_qc_parameters <- function() {
  data_dir <- "data"
  cel_dir <- "celfiles"
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

  raw <- list()
  norm <- list()
  norm_f <- list()
  norm_filt <- list()

  # boxplot ----

  raw[["boxplot"]] <- list(
    group = "Group",
    group.color = "Colors",
    samplenames = "ShortName",
    cex.lab = .6,
    outputDir = create_dir(results_dir, "raw", "boxplot"),
    label = "RawData"
  )

  norm[["boxplot"]] <- purrr::list_modify(
    raw[["boxplot"]],
    label = "NormData",
    outputDir = create_dir(results_dir, "norm", "boxplot"),
  )

  norm_f[["boxplot"]] <- purrr::list_modify(
    raw[["boxplot"]],
    label = "NormFData",
    outputDir = create_dir(results_dir, "norm_f", "boxplot"),
  )

  # pca ----

  raw[["pca"]] <- list(
    pca2D_factors = targets_fact,
    scale = FALSE,
    col.group = "Colors",
    colorlist = colorlist,
    outputDir = create_dir(results_dir, "raw", "pca"),
    label = "Raw"
  )

  norm[["pca"]] <- purrr::list_modify(
    raw[["pca"]],
    label = "Normalized",
    outputDir = create_dir(results_dir, "norm", "pca")
  )

  norm_f[["pca"]] <- purrr::list_modify(
    raw[["pca"]],
    label = "NormalizedF",
    outputDir = create_dir(results_dir, "norm_f", "pca")
  )

  # batch pca ----

  norm[["batch_pca"]] <- purrr::list_modify(
    raw[["pca"]], # We copy from pca as we are only adding the batch
    label = "NormDataBatch",
    batchRemove = TRUE,
    batchFactors = "Batch",
    size = 1.5,
    glineas = 0.25,
    outputDir = create_dir(results_dir, "norm", "batch_pca")
  )

  norm_f[["batch_pca"]] <- purrr::list_modify(
    norm[["batch_pca"]],
    label = "NormFData",
    outputDir = create_dir(results_dir, "norm_f", "batch_pca")
  )

  # pvca ----

  norm[["pvca"]] <- list(
    pct_threshold = .6,
    factors = targets_fact,
    label = "NormData",
    inputDir = data_dir,
    outputDir = create_dir(results_dir, "norm", "pvca"),
    summaryFN = summary_fn
  )

  norm_f[["pvca"]] <- purrr::list_modify(
    norm[["pvca"]],
    label = "NormFData",
    outputDir = create_dir(results_dir, "norm_f", "pvca"),
  )

  # hc ----

  raw[["hc"]] <- list(
    hclust.method = "average",
    cexRow = 0.6,
    cexCol = 0.6,
    rect = TRUE,
    numclusters = 2,
    outputDir = create_dir(results_dir, "raw", "hc"),
    label = "RawData"
  )

  norm[["hc"]] <- purrr::list_modify(
    raw[["hc"]],
    label = "NormData",
    outputDir = create_dir(results_dir, "norm", "hc")
  )

  norm_f[["hc"]] <- purrr::list_modify(
    raw[["hc"]],
    label = "NormFData",
    outputDir = create_dir(results_dir, "norm_f", "hc")
  )

  # aqm ----

  raw[["aqm"]] <- list(
    outdir = create_dir(results_dir, "raw", "QCDir.Raw"),
    force = TRUE,
    intgroup = targets_fact
  )

  norm[["aqm"]] <- purrr::list_modify(
    raw[["aqm"]],
    outdir = create_dir(results_dir, "norm", "QCDir.Norm"),
  )

  norm_f[["aqm"]] <- purrr::list_modify(
    raw[["aqm"]],
    outdir = create_dir(results_dir, "norm_f", "QCDir.NormF"),
  )

  # normalize ----

  norm[["normalize"]] <- list(
    annotPkg = "clariomsmousetranscriptcluster.db",
    norm.method = "rma",
    outputFN = "Normalized.All",
    outputDir = create_dir(results_dir, "norm", "normalize")
  )
  norm_f[["normalize"]] <- purrr::list_modify(
    norm[["normalize"]],
    outputFN = "NormalizedF.All",
    outputDir = create_dir(results_dir, "norm_f", "normalize")
  )


  norm[["annotate"]] <- list(
    annotPkg = "clariomsmousetranscriptcluster.db",
    outputFN = "Annotations.AllGenes",
    saveHTML = TRUE,
    title = "Annotations for all genes",
    outputDir = create_dir(results_dir, "norm", "annotate")
  )
  norm_f[["annotate"]] <- purrr::list_modify(
    norm[["annotate"]],
    outputFN = "AnnotationsF.AllGenes",
    title = "Annotations for all genes (No outliers)",
    outputDir = create_dir(results_dir, "norm_f", "annotate")
  )
  norm_filt[["annotate"]] <- purrr::list_modify(
    norm[["annotate"]],
    outputFN = "Annotations.Filtered",
    title = "Annotations for Filtered genes",
    outputDir = create_dir(results_dir, "norm_filt", "annotate")
  )

  filtering <- list(
    outputFN = "Normalized.Filtered",
    outputDir = create_dir(results_dir, "norm_filt", "filtering"),
    summaryFN = "ResultsSummary_Load-QC-Norm-Filt.txt",
    feature.exclude = "^AFFX",
    require.entrez = TRUE,
    remove.dupEntrez = TRUE,
    var.filter = FALSE,
    var.func = "IQR",
    var.cutoff = 0.65,
    filterByQuantile = TRUE,
    require.GOBP = FALSE,
    require.GOCC = FALSE,
    require.GOMF = FALSE
  )

  outlier_samples <- "CMP.2"

  sdplot <- list(
    var_thr = 0.65,
    label = "NormFData",
    outputDir = create_dir(results_dir, "norm_filtered", "sdplot")
  )

  list(
    targets_file = targets_file,
    cel_dir = cel_dir,
    targets_fact = targets_fact,
    results_dir = results_dir,
    color_list = colorlist,
    raw = raw,
    norm = norm,
    norm_f = norm_f,
    norm_filt = norm_filt,
    outlier_samples = outlier_samples,
    sdplot = sdplot,
    filtering = filtering
  )
}

qcl_par <- get_load_qc_parameters()

# dea parameters ----

get_dea_par <- function() {
  dataDirN <- "data" # name of data directory
  resultsDir <- "results" # name of results directory
  annotPackage <- "clariomsmousetranscriptcluster.db"
  resultsSummFN <- "ResultsSummary_DEA.txt"

  # Contrasts matrix
  contrastsv <- c(
    "PD1vsCTL = PD1 - CTL",
    "CMPvsCTL = CMP - CTL",
    "PD1vsCMP = PD1 - CMP"
  ) # Comparisons to be performed

  ncomp <- length(contrastsv)

  lmdesign <- list(
    sampleNames = "ShortName",
    group = "Group", # Main variable containing group information to include in Design matrix (eg. "Group")
    cov.fact = "Batch", # Covariates to be included in the linear model as cofactors (eg. "Batch", "Sex", ...). If not to consider any variable, set NULL.
    cov.num = NULL, # Covariates to be included in the linear model as numeric variables (eg. "Age", ...). If not to consider any variable, set NULL.
    fmla = NULL,
    summaryFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "lm", "lmdesign")
  )

  lmfit <- list(
    block.cor = NULL,
    summaryFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "lm", "lmfit")
  )

  lmcompare <- list()
  lmcompare[["main"]] <- list(
    contrasts = contrastsv,
    moderated = TRUE,
    summaryFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "lm", "lmcompare")
  )

  toptab <- list()
  toptab[["BH"]] <- list(
    padjust.method = "fdr",
    html_report = FALSE, # It does not work properly in maUEB
    html_ntop = 500,
    html_group = "Group",
    html_padjust_method = c("none", "BH", "BH"), # p-value adjust method for HTML top table (eg. "BH", "none") for each coefficient. Select "none" if you want to show the raw pvalue in html table
    outputDir = create_dir(resultsDir, "lm", "toptab_BH")
  )

  toptab[["FDR"]] <- purrr::list_modify(toptab[["BH"]], padjust.method = "fdr", html_padjust_method = purrr::zap(), outputDir = create_dir(resultsDir, "lm", "toptab_FDR"))

  sum_ngc1 <- list(
    B_thr = c(0),
    Pval_thr = c(0.01, 0.05, 0.1),
    adjPval_thr = c(0.01, 0.05, 0.15, 0.25),
    logFC_thr = 0,
    extended = TRUE,
    B_thr.e = c(0),
    Pval_thr.e = c(0.01, 0.05, 0.1),
    adjPval_thr.e = c(0.01, 0.05, 0.15, 0.25),
    logFC_thr.e = c(0, 0.5, 1, 1.5, 2),
    outputDir = create_dir(resultsDir, "lm", "sum_ngc1")
  )

  volcano <- list(
    volc_logFC = rep(1, ncomp),
    volc_pval = rep("adj.P.Val", ncomp),
    volc_pval.thr = rep(0.05, ncomp),
    volc_x0 = rep(-3, ncomp),
    volc_x1 = rep(+3, ncomp),
    volc_y0 = rep(0, ncomp),
    volc_y1 = rep(10, ncomp),
    n = 6,
    cols = 2,
    outputDir = create_dir(resultsDir, "lm", "volcano"),
    label = ""
  )

  list(
    lmdesign = lmdesign,
    lmfit = lmfit,
    lmcompare = lmcompare,
    toptab = toptab,
    sum_ngc1 = sum_ngc1,
    volcano = volcano
  )
}

dea_par <- get_dea_par()
# Multiple Comparison parameters ----
get_mc_par <- function() {
  dataDirN <- "data" # name of data directory
  resultsDir <- "results" # name of results directory

  # Resultsfile
  resultsSummFN <- "ResultsSummary_MC.txt" # A summary of results obtained in each step will be saved in a text file

  venn_upset <- list(
    label = "EffectTreatment", # name of set of comparisons (without spaces/dots/_)
    venn_comparNames = "EffectTreatment",
    colFeat = "Gene.Symbol",
    colPval = "adj.P.Val", # "P.Value" or "adj.P.Val" (for each set of comparisons)
    pval = 0.05,
    colFC = "logFC",
    FC = 1,
    include = "abs",
    pltR = TRUE,
    pltPdf = TRUE,
    pltPng = FALSE,
    venn = TRUE,
    eul = FALSE,
    saveTables = TRUE,
    upsetPlot = FALSE,
    saveTables_extended = TRUE,
    trans = 0.5,
    cex1 = 1,
    rotation = 0,
    position = c(0, 0, 180),
    cex2 = 1,
    resultsSummFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "mc", "venn_upset")
  )

  heatmap <- list(
    featureCol = "Gene.Symbol",
    hm_comparNames = c("EffectTreatment", "EffectPD1"), # name of set of comparisons (without spaces/dots/_). Usually the same as used for Venn (eg. venn_comparNames),
    hm_compar = list(1:3, 1), # list with set of comparisons to select genes to do the Heatmap, positions (or names) according to 'listofcoef' or 'listofTopNames'. Usually the same as used for Venn (eg. venn_compar),
    hm_groupexclude = list(c(), c("CMP")), # samples to exclude in each set of comparisons for Heatmap (length must be equal to number of comparisons). If no sample is to be excluded put c().
    hm_pval = c("adj.P.Val", "P.Value"), # "P.Value" or "adj.P.Val" (for each set of comparisons)
    hm_pval.thr = c(0.05, 0.05), # pvalue threshold  (for each set of comparisons)
    hm_logFC = "logFC",
    hm_logFC.thr = c(1, 1), # logFC threshold (for each set of comparisons)
    hm_palette = colorRampPalette(c("blue", "white", "red"))(n = 199), # si valor alto no colores en la key
    hm_clustCol = TRUE, # whether to perform hierarchical clustering of samples in  heatmap (column dendogram) (eg. TRUE/FALSE)
    hm_clustCol.dist = "euclidian", # required if hm_clustCol=TRUE: method for distance calculation in hierarchical clustering of samples (default in heatmap.2 is "euclidian")
    hm_clustCol.method = "complete", # required if hm_clustCol=TRUE: method for hierarchical clustering (default in heatmap.2 is "complete", other good may be "average")
    hm_clustRow.cor = TRUE, # whether to cluster gene expression based on distance correlation (otherwise distance matrix is calcluated based on euclidian distance (default in heatmap.2 function))
    batcheffect = FALSE, # whether to correct expression values for batch effect (eg. TRUE/FALSE)
    batchcolName = NULL,
    batchcolName2 = NULL,
    batchNumcolName = NULL,
    hm_plots_interactive = TRUE,
    resultsSummFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "mc", "heatmap")
  )
  list(
    venn_upset = venn_upset,
    heatmap = heatmap
  )
}

mc_par <- get_mc_par()
# abs parameters ----
get_abs_par <- function() {
  # Directories
  dataDir <- "data" # name of data directory
  resultsDir <- "results" # name of results directory

  # Resultsfile
  resultsSummFN <- "ResultsSummary_ABS-GSEA.txt" # A summary of results obtained in each step will be saved in a text file

  gsea <- list(
    go = list(),
    reactome = list()
  )

  gsea[["go"]][["unfilt"]] <- list(
    annotPackage = "clariomsmousetranscriptcluster.db",
    organism_annot = "org.Mm.eg.db", # organism annotation package (eg. "org.Hs.eg.db", "org.Mm.eg.db") (required for GO)
    ranking_metric = "logFC",
    geneColname = "EntrezID",
    keyType = "ENTREZID", # amb majuscules, corresponent a camp d taula d'anotacions orgdb
    GO_minSetSize = 3,
    GO_maxSetSize = 500,
    resultsSummFN = resultsSummFN,
    GO_categories = "BP",
    saveRda = TRUE,
    saveGMT = TRUE,
    label = "",
    outputDir = create_dir(resultsDir, "abs", "gsea", "go", "unfilt"),
    rdaDir = dataDir,
    gmtDir = create_dir(resultsDir, "abs", "gsea", "go", "unfilt")
  )
  gsea[["go"]][["filt"]] <- list(
    GO_pvalcutoff = list(c(BP = .05)), # list(rep(0.05, 3), rep(0.05, 3), rep(0.05, 3)), # in this order: BP, CC, MF
    GO_pAdjustMethod = list(c(BP = "BH")), # list(rep("BH", 3), rep("BH", 3), rep("BH", 3)), # (eg. "BH", "bonferroni" or "none")
    saveTables = TRUE,
    GOPlots = TRUE,
    label = "",
    resultsSummFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "abs", "gsea", "go", "filt")
  )
  gsea[["go"]][["GO_comparisons"]] <- 1 # c(1:3) # comparisons to be analysed, indicated by the corresponding position in 'listofcoef' or 'listofTopNames' (eg. 1:10 or c(1,2,8))

  gsea[["reactome"]][["unfilt"]] <- list(
    annotPackage = "clariomsmousetranscriptcluster.db",
    organism = "mouse",
    organism_annot = "org.Mm.eg.db",
    Reac_minSetSize = 3,
    Reac_maxSetSize = 500,
    resultsSummFN = resultsSummFN,
    saveRda = TRUE,
    saveGMT = TRUE,
    label = "",
    outputDir = create_dir(resultsDir, "abs", "gsea", "reactome", "unfilt"),
    rdaDir = dataDir,
    gmtDir = create_dir(resultsDir, "abs", "gsea", "reactome", "unfilt")
  )

  gsea[["reactome"]][["reactome_comparisons"]] <- 1

  gsea[["reactome"]][["filt"]] <- list(
    Reac_pvalcutoff = .05,
    Reac_pAdjustMethod = "BH",
    saveTables = TRUE,
    ReacPlots = TRUE,
    label = "",
    resultsSummFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "abs", "gsea", "reactome", "filt")
  )

  # ORA ----
  ora <- list(
    go = list(),
    reactome = list()
  )
  ora[["go"]][["GO_comparisons"]] <- 1
  ora[["go"]][["unfilt"]] <- list(
    universe = NULL,
    geneColname = "EntrezID",
    readable = TRUE,
    selected.pval.go = c("P.Value"),
    selected.pvalthr.go = c(0.05),
    selected.logFC.go = c(1),
    col_logFC.go = c("logFC"),
    sign_logFC.go = c("abs"),
    organism_annot = "org.Mm.eg.db",
    keyType = "ENTREZID",
    GO_categories = list("BP"),
    GO_minSetSize = 3,
    GO_maxSetSize = 500,
    resultsSummFN = resultsSummFN,
    saveRda = TRUE,
    saveGMT = TRUE,
    label = "",
    outputDir = create_dir(resultsDir, "abs", "ora", "go", "unfilt"),
    rdaDir = dataDir,
    gmtDir = create_dir(resultsDir, "abs", "ora", "go", "unfilt")
  )

  ora[["go"]][["filt"]] <- list(
    sign = c(1, -1),
    GO_pvalcutoff = list(c(BP = 0.15)),
    GO_pAdjustMethod = list(c(BP = "BH")),
    saveTables = TRUE,
    GOPlots = TRUE,
    label = "label",
    resultsSummFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "abs", "ora", "go", "filt")
  )

  ora[["reactome"]][["unfilt"]] <- list(
    universe = NULL,
    col_entrez = "EntrezID",
    readable = TRUE,
    selected.pval.reac = c("P.Value"),
    selected.pvalthr.reac = c(0.05),
    selected.logFC.reac = c(1),
    col_logFC.reac = c("logFC"),
    sign_logFC.reac = c("abs"),
    organism = "mouse",
    organism_annot = "org.Mm.eg.db",
    keyType = "ENTREZID",
    Reac_minSetSize = 3,
    Reac_maxSetSize = 500,
    resultsSummFN = resultsSummFN,
    saveRda = TRUE,
    saveGMT = TRUE,
    label = "label",
    outputDir = create_dir(resultsDir, "abs", "ora", "reactome", "unfilt"),
    rdaDir = dataDir,
    gmtDir = create_dir(resultsDir, "abs", "ora", "reactome", "unfilt")
  )

  ora[["reactome"]][["reactome_comparisons"]] <- 1
  ora[["reactome"]][["filt"]] <- list(
    sign = c(1, -1),
    Reac_pvalcutoff = 0.15,
    Reac_pAdjustMethod = "BH",
    saveTables = TRUE,
    ReacPlots = TRUE,
    label = "",
    resultsSummFN = resultsSummFN,
    outputDir = create_dir(resultsDir, "abs", "ora", "reactome", "filt")
  )

  list(
    gsea = gsea,
    ora = ora
  )
}

abs_par <- get_abs_par()


# Set target options: ----
tar_option_set(
  packages = c(
    "tibble", "oligoClasses", "annaffy",
    "Biostrings", "XVector", "GenomeInfoDb",
    "GGally", qcl_par$norm$annotate$annotPkg, "oligo",
    "limma", "maUEB", "grid"
  ), # packages that your targets need to run
  format = "rds", # default storage format
  # Set other options as needed.
  memory = "transient", garbage_collection = TRUE,
  error = "null"
)

options(clustermq.scheduler = "multicore")
# Run the R scripts in the R/ folder with your custom functions:
tar_source()


# Helper for debugging do calls with !!
evaler <- function(x) {
  rlang::eval_tidy(rlang::enexpr(x))
}

# target list ----

# parameters in the call from parameter lists are preceded by !!. This is done so the target depend exclusively on its particular parameters
# and not on ht whole list. e.g.: targets_file target depends only on the entry !!qcl_par$targets_file and not on the whole qcl_par. This
# avoids that changing the parameter list in another entry (e.g.: cel_files) invalidates the targets_file target.

# Running specific targets can be done modifying the targets::tar_make call, indicating the specific targets that we want to make while avoiding
# creating the whole list. (e.g.: If we want to create just the volcano plot we can pass that specific target and targets will run the minimum amount
# of targets needed to make the volcano plot)

list(
  # FILES ----
  tar_target(
    name = targets_file,
    command = !!qcl_par$targets_file,
    format = "file"
  ),
  tar_target(
    name = targets,
    command = read_targets_alias(
      targets_file,
      !!qcl_par$targets_fact
    )
  ),
  tar_target(
    name = cel_files,
    command = get_cel_files(
      targets,
      !!qcl_par$cel_dir
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
    command = do.call(
      what = maUEB::qc_boxplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$raw$boxplot,
        data = eset_raw
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_raw_load_pca1_file,
    command = do.call(
      what = maUEB::qc_pca1 %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$raw$pca,
        data = Biobase::exprs(eset_raw),
        targets = targets,
        names = targets$ShortName
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_raw_hc,
    command = do.call(
      what = maUEB::qc_hc %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$raw$hc,
        data = Biobase::exprs(eset_raw),
        names = targets$ShortName
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_raw_arrayQM,
    command = do.call(
      what = arrayQualityMetrics::arrayQualityMetrics %>% append_args_as_attr("outdir"),
      args = purrr::list_modify(
        !!qcl_par$raw$aqm,
        expressionset = eset_raw,
        reporttitle = "arrayQualityMetrics report for eset_raw"
      )
    )
  ),
  # NORM ----
  tar_target(
    name = eset_norm,
    command = do.call(
      what = maUEB::normalization %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm$normalize,
        data = eset_raw
      )
    )
  ),
  tar_target(
    name = norm_annot,
    command = do.call(
      what = maUEB::save_annotations %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm$annotate,
        data = rownames(eset_norm)
      )
    )
  ),
  tar_target(
    name = qc_norm_boxplot_file,
    command = do.call(
      what = maUEB::qc_boxplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm$boxplot,
        data = eset_norm
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_load_pca1_file,
    command = do.call(
      what = maUEB::qc_pca1 %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm$pca,
        data = Biobase::exprs(eset_norm),
        targets = targets,
        names = targets$ShortName
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_corr_batch_load_pca1,
    command = do.call(
      what = maUEB::qc_pca1 %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm$batch_pca,
        data = Biobase::exprs(eset_norm),
        targets = targets,
        names = targets$ShortName
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_hc,
    command = do.call(
      what = maUEB::qc_hc %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm$hc,
        data = Biobase::exprs(eset_norm),
        names = targets$ShortName
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_arrayQM,
    command = do.call(
      what = arrayQualityMetrics::arrayQualityMetrics %>% append_args_as_attr("outdir"),
      args = purrr::list_modify(
        !!qcl_par$norm$aqm,
        expressionset = eset_norm,
        reporttitle = "arrayQualityMetrics report for eset_norm"
      )
    )
  ),
  # BLOCKED UNTIL https://github.com/uebvhir/maUEB/issues/3 IS SOLVED
  # tar_target(
  #   name = qc_norm_pcpva,
  #   command = do.call(
  #     what = maUEB::qc_pvca %>% as_list_of_files("outputDir"),
  #     args = purrr::list_modify(
  #       !!qcl_par$norm$pvca,
  #       data = eset_norm
  #     )
  #   ),
  #   format = "file"
  # ),
  # SAMPLE REMOVAL ----
  tar_target(
    name = eset_raw_f,
    command = remove_samples(
      eset_raw,
      !!qcl_par$outlier_samples
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
    command = do.call(
      what = maUEB::normalization %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm_f$normalize,
        data = eset_raw_f
      )
    )
  ),
  # tar_target( # DO WE ANNOTATE WITH NO OUTLIERS?
  #   name = norm_annot_f,
  #   command = do.call(
  #     what = save_annotations_alias,
  #     args = purrr::list_modify(
  #       !!qcl_par$annotate$norm_f,
  #       data = rownames(eset_norm_f)
  #     )
  #   )
  # ),
  tar_target(
    name = qc_norm_f_boxplot_file,
    command = do.call(
      what = maUEB::qc_boxplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm_f$boxplot,
        data = eset_norm_f
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_f_load_pca1_file,
    command = do.call(
      what = maUEB::qc_pca1 %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm_f$pca,
        data = Biobase::exprs(eset_norm_f),
        targets = targets_f,
        names = targets_f$ShortName
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_f_corr_batch_load_pca1_file,
    command = do.call(
      what = maUEB::qc_pca1 %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm_f$batch_pca,
        data = Biobase::exprs(eset_norm_f),
        targets = targets_f,
        names = targets_f$ShortName
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_f_hc_file,
    command = do.call(
      what = maUEB::qc_hc %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm_f$hc,
        data = Biobase::exprs(eset_norm_f),
        names = targets_f$ShortName
      )
    ),
    format = "file"
  ),
  tar_target(
    name = qc_norm_f_arrayQM,
    command = do.call(
      what = arrayQualityMetrics::arrayQualityMetrics %>% append_args_as_attr("outdir"),
      args = purrr::list_modify(
        !!qcl_par$norm_f$aqm,
        expressionset = eset_norm_f,
        reporttitle = "arrayQualityMetrics report for eset_norm_f"
      )
    )
  ),
  # BLOCKED UNTIL https://github.com/uebvhir/maUEB/issues/3 IS SOLVED
  # tar_target(
  #   name = qc_norm_f_pcpva,
  #   command = evaler(maUEB::qc_pvca(data = eset_norm, factors = !!par$pvca$norm$fact, targets = !!par$pvca$norm$targets_file_name, pct_threshold = !!par$pvca$norm$threshold, label = !!par$pvca$norm$label, outputDir = !!par$pvca$norm$output_dir, summaryFN = !!par$pvca$norm$summary_fn)) # nolintr
  # ),

  # FILTERING ----
  tar_target(
    name = sdplot_file,
    command = do.call(
      what = maUEB::sdplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$sdplot,
        data = exprs(eset_norm_f)
      )
    ),
    format = "file"
  ),
  tar_target(
    name = eset_filtered,
    command = do.call(
      what = maUEB::filtering %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$filtering,
        data = eset_norm_f
      )
    )
  ),
  tar_target(
    name = norm_filtered_annot,
    command = do.call(
      what = maUEB::save_annotations %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!qcl_par$norm_filt$annotate,
        data = rownames(eset_norm_f)
      )
    )
  ),
  # BLOCK 02 ----
  #
  tar_target(
    name = design,
    command = do.call(
      what = maUEB::dea_lmdesign %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!dea_par$lmdesign,
        data = eset_filtered,
        targets = pData(eset_filtered)
      )
    )
  ),
  tar_target(
    name = fit,
    command = do.call(
      what = maUEB::dea_lmfit %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!dea_par$lmfit,
        data = eset_filtered,
        targets = pData(eset_filtered),
        design = design
      )
    )
  ),
  tar_target(
    name = compare_main,
    command = do.call(
      what = maUEB::dea_compare %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!dea_par$lmcompare$main,
        fit = fit,
        design = design
      )
    )
  ),
  tar_target(
    name = toptable_BH,
    command = do.call(
      what = maUEB::dea_toptab %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!dea_par$toptab$BH,
        fit.main = compare_main,
        eset = eset_filtered,
        listofcoef = colnames(compare_main)
      )
    )
  ),
  tar_target(
    name = toptable_FDR,
    command = do.call(
      what = maUEB::dea_toptab1 %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!dea_par$toptab$FDR,
        fit.main = compare_main,
        eset = eset_filtered,
        listofcoef = colnames(compare_main)
      )
    )
  ),
  tar_target(
    name = numGenes_all,
    command = do.call(
      what = maUEB::dea_summary_ngc1 %>% append_args_as_attr("outputDir"),
      args = purrr::list_modify(
        !!dea_par$sum_ngc1,
        listofcoef = colnames(compare_main),
        listofcsv = toptable_BH
      )
    )
  ),
  tar_target(
    name = volcano_file,
    command = do.call(
      what = maUEB::dea_volcanoplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!dea_par$volcano,
        listofcoef = colnames(compare_main),
        listofcsv = toptable_BH
      )
    ),
    format = "file"
  ),
  # MC ----
  tar_target(
    name = mc_venn_upset_file,
    command = do.call(
      what = maUEB::mc_venn_upset %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!mc_par$venn_upset,
        listofcsv = toptable_BH,
        namescomp = names(toptable_BH),
        colors = rainbow(length(names(toptable_BH)))
      )
    ),
    format = "file"
  ),
  tar_target(
    name = mc_heatmap_file,
    command = do.call(
      what = maUEB::mc_hm %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!mc_par$heatmap,
        listofcsv = toptable_BH,
        hm_dataset = toptable_BH[[1]],
        targets = targets_f
      )
    ),
    format = "file"
  ),
  # ABS GSEA ----
  tar_target(
    name = abs_gsea_go_unfilt,
    command = do.call(
      what = maUEB::abs_gsea_GOunfilt %>% append_args_as_attr(c("outputDir", "gmtDir")),
      args = purrr::list_modify(
        !!abs_par$gsea$go$unfilt,
        listofTopNamed = toptable_BH,
        namescomp = names(toptable_BH)[!!abs_par$gsea$go$GO_comparisons]
      )
    )
  ),
  tar_target(
    name = abs_gsea_go_filt_file,
    command = do.call(
      what = maUEB::abs_gsea_GOfiltplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!abs_par$gsea$go$filt,
        GSEAresGO = abs_gsea_go_unfilt
      )
    ),
    format = "file"
  ),
  tar_target(
    name = abs_gsea_reactome_pa_unfilt,
    command = do.call(
      what = abs_gsea_ReactomePAunfilt_alias %>% append_args_as_attr(c("outputDir", "gmtDir")),
      args = purrr::list_modify(
        !!abs_par$gsea$reactome$unfilt,
        listofTopNamed = toptable_BH,
        namescomp = names(toptable_BH)[abs_par$gsea$reactome$reactome_comparisons]
      )
    )
  ),
  tar_target(
    name = abs_gsea_reactome_pa_filt_file,
    command = do.call(
      what = maUEB::abs_gsea_ReactomePAfiltplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!abs_par$gsea$reactome$filt,
        GSEAresReac = abs_gsea_reactome_pa_unfilt
      )
    ),
    format = "file"
  ),
  # ABS ORA ----
  tar_target(
    name = abs_ora_go_unfilt,
    command = do.call(
      what = maUEB::abs_ora_GOunfilt %>% append_args_as_attr(c("outputDir", "gmtDir")),
      args = purrr::list_modify(
        !!abs_par$ora$go$unfilt,
        listofTopNamed = toptable_BH,
        namescomp = names(toptable_BH)[!!abs_par$ora$go$GO_comparisons]
      )
    )
  ),
  tar_target(
    name = abs_ora_go_filt_file,
    command = do.call(
      what = maUEB::abs_ora_GOfiltplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!abs_par$ora$go$filt,
        ORAresGO = abs_ora_go_unfilt
      )
    ),
    format = "file"
  ),
  tar_target(
    name = abs_ora_reactome_pa_unfilt,
    command = do.call(
      what = maUEB::abs_ora_ReactomePAunfilt %>% append_args_as_attr(c("outputDir", "gmtDir")),
      args = purrr::list_modify(
        !!abs_par$ora$reactome$unfilt,
        listofTopNamed = toptable_BH,
        namescomp = names(toptable_BH)[abs_par$ora$reactome$reactome_comparisons]
      )
    )
  ),
  tar_target(
    name = abs_ora_reactome_pa_filt_file,
    command = do.call(
      what = maUEB::abs_ora_ReactomePAfiltplot %>% as_list_of_files("outputDir"),
      args = purrr::list_modify(
        !!abs_par$ora$reactome$filt,
        ORAresReac = abs_ora_reactome_pa_unfilt
      )
    ),
    format = "file"
  ),
  tarchetypes::tar_render(report, "report.Rmd", output_dir = "results")
)
