#' Calculate p-value for provided data
#'
#' @param df_in data.frame
#' @param data_col Column containing data to use for test
#' @param cluster_col Column containing clusters to compare
#' @param group_col Column containing labels to use for grouping data, a
#' p-value will be calculated for each group. Each group should have at least
#' two clusters from `cluster_col`
#' @param method Method to use for calculating p-value, possible values include:
#'
#' * 't', two sample t-test
#' * 'wilcox', Wilcoxon rank sum test
#' * 'anova', Kruskal-Wallis rank sum test
#'
#' @noRd
.calc_pvalue <- function(df_in, data_col, cluster_col, group_col,
                         method = NULL, adj_method = "bonferroni",
                         file = NULL) {

  # Set method based on number of clusters for comparison
  n_clsts <- n_distinct(df_in[[cluster_col]])

  if (n_clsts > 2 && !identical(method, "kruskal")) {
    cli::cli_warn(
      "p-values will be calculated using the Kruskal-Wallis test
       since more than two groups are being compared."
    )
  }

  method <- dplyr::case_when(
    n_clsts == 1 ~ "none",
    n_clsts > 2  ~ "kruskal",
    TRUE         ~ method %||% "t"
  )

  # Calculate p-values
  if (identical(method, "edgeR")) {
    res <- .calc_edgeR(df_in, data_col, cluster_col, group_col)

    adj_method <- NULL

  } else {
    res <- .calc_p(df_in, data_col, cluster_col, group_col, method = method)
  }

  # Multiple testing correction
  p_clmn <- "p_value"

  if (!is.null(adj_method)) {
    p_res <- dplyr::distinct(res, !!!syms(c(group_col, p_clmn)))

    p_res <- dplyr::mutate(
      p_res, p_adj = stats::p.adjust(!!sym(p_clmn), method = adj_method)
    )

    p_res  <- purrr::set_names(p_res$p_adj, p_res[[group_col]])
    res    <- dplyr::mutate(res, p_adj = p_res[as.character(!!sym(group_col))])
    p_clmn <- "p_adj"
  }

  # Write table
  if (!is.null(file)) {
    p_tbl <- dplyr::arrange(res, !!!syms(c(p_clmn, group_col, cluster_col)))
    p_tbl <- dplyr::mutate(p_tbl, method = method)

    readr::write_csv(p_tbl, file, progress = FALSE)
  }

  res <- dplyr::rename(res, p_value = !!sym(p_clmn))

  res
}


#' Calculate p-value for specified groups
#' @noRd
.calc_p <- function(df_in, data_col, cluster_col, group_col,
                    method = "t") {

  # Helper to perform selected test
  .calc_test <- function(x, g, method) {
    p <- split(x, g, drop = TRUE)

    if (length(p) != 2) {
      cli::cli_abort(
        "There must be two groups present to perform test.", .internal = TRUE
      )
    }

    p <- method(p[[1]], p[[2]])
    p <- p$p.value

    p
  }

  # Return NAs if method is "none"
  if (identical(method, "none")) {
    res <- dplyr::mutate(df_in, p_value = as.numeric(NA))

    return(res)
  }

  # Set testing method
  p_fn <- switch(
    method,
    t       = function(x, ...) .calc_test(x, ..., method = stats::t.test),
    wilcox  = function(x, ...) .calc_test(x, ..., method = stats::wilcox.test),
    kruskal = function(x, ...) (stats::kruskal.test(x, ...))$p.value,
  )

  # Calculate p-values
  df_in <- dplyr::group_by(df_in, !!sym(group_col))

  res <- dplyr::mutate(
    df_in, p_value = p_fn(!!sym(data_col), !!sym(cluster_col))
  )

  res <- dplyr::ungroup(res)

  if (any(is.na(res$p_value))) {
    cli::cli_warn("`NA`s were returned when calculating p-values.")
  }

  res
}


#' Test for differential abundance using edgeR
#' @noRd
.calc_edgeR <- function(df_in, data_col, cluster_col, group_col) {

  # Format input data
  dat <- dplyr::select(df_in, all_of(c(data_col, cluster_col, group_col)))
  dat <- dplyr::filter(dat, !is.na(as.character(!!sym(group_col))))

  dat <- dplyr::group_by(dat, !!!syms(c(group_col, cluster_col)))

  dat <- dplyr::mutate(
    dat,
    .pool   = row_number(),
    .sample = paste0(!!sym(cluster_col), "_", .pool)
  )

  dat <- dplyr::ungroup(dat)

  # data.frame of sample info
  sam_info <- dplyr::distinct(dat, !!sym(cluster_col), .sample, .pool)

  # Counts matrix
  mat <- dplyr::select(dat, all_of(c(group_col, ".sample", data_col)))

  mat <- tidyr::pivot_wider(
    mat, names_from = .sample,
    values_from = !!sym(data_col), values_fill = 0
  )

  mat <- tibble::column_to_rownames(mat, group_col)
  mat <- as.matrix(mat)

  # Create DE list
  de_l <- edgeR::DGEList(mat, samples = sam_info)

  # Remove groups with too few cells
  keep <- edgeR::filterByExpr(de_l, group = de_l$samples[[cluster_col]])
  de_l <- de_l[keep, ]

  # Model formula
  frm <- paste0("~ factor(.pool) + factor(", cluster_col, ")")

  design <- stats::model.matrix(stats::as.formula(frm), de_l$samples)

  # Estimage dispersion
  de_l <- edgeR::estimateDisp(de_l, design, trend = "none")

  # Fit glm
  de_fit <- edgeR::glmQLFit(de_l, design, robust = TRUE, abundance.trend = FALSE)

  # Test for differential abundance
  pvals <- edgeR::glmQLFTest(de_fit, coef = ncol(design))
  pvals <- pvals$table
  pvals <- purrr::set_names(pvals$PValue, rownames(pvals))

  res <- dplyr::mutate(
    df_in, p_value = pvals[as.character(!!sym(group_col))]
  )

  if (any(is.na(res$p_value))) {
    cli::cli_warn("`NA`s were returned when calculating p-values.")
  }

  res
}


#' Format p-values for plotting
#' @noRd
.format_pvalue <- function(p, digits = 1, cutoffs = NULL) {

  # Set p label based on vector of cutoffs
  if (!is.finite(p)) return(as.character(NA))

  if (!is.null(cutoffs)) {
    if (any(duplicated(cutoffs))) {
      cli::cli_abort("Cutoff values for p_label must be unique.")
    }

    # Set default labels when not provided by user
    if (is.null(names(cutoffs))) {
      cutoffs <- sort(cutoffs, decreasing = TRUE)

      names(cutoffs) <- purrr::imap_chr(
        cutoffs, ~ paste0(rep("*", .y), collapse = "")
      )
    }

    cutoffs <- sort(cutoffs)
    p_label <- as.character(NA)

    for (val in names(cutoffs)) {
      if (p < cutoffs[val]) {
        p_label <- val

        break()
      }
    }

    # Treat "value" as a keyword that will allow user to display actual
    # p-value for a certain cutoff
    # All custom labels need to be wrapped in quotes for parsing
    if (!identical(p_label, "value")) {
      if (!is.na(p_label)) p_label <- paste0("\'", p_label, "\'")

      return(p_label)
    }
  }

  # Format p-value label
  p <- scales::label_scientific(digits = digits)(p)

  ex <- .str_extract_all(p, "[+\\-][0-9]+$")

  p <- sub(paste0("\\", ex, "$"), "", p)

  ex <- as.numeric(ex)
  ex <- as.character(ex)

  p <- sub("e", "*x*10^", p)
  p <- paste0(p, ex)

  p
}


#' Parse label params
#' divide params based on prefix, e.g. 'n.' or 'p.'
#' this allows user to adjust aesthetics for a specific label
#' @noRd
.parse_label_params <- function(l, prfxs = c("n", "p"), sep = ".") {
  nms   <- names(l)
  prfxs <- purrr::set_names(paste0("^", prfxs, "\\", sep), prfxs)
  idxs  <- purrr::map(prfxs, ~ grepl(.x, nms))

  shared <- purrr::reduce(idxs, ~ !.x & !.y)
  shared <- l[shared]

  res <- purrr::imap(idxs, ~ {
    prfx <- .y
    prms <- l[.x]

    names(prms) <- purrr::map_chr(names(prms), ~ sub(prfxs[[prfx]], "", .x))

    prms <- append(prms, shared)
    prms <- .standardize_aes(prms)

    prms
  })

  res
}







# library(edgeR)
# library(djvdj)
# library(tidyverse)

# # Bioconductor example ----
#
# abundances <- table(merged$celltype.mapped, merged$sample)
# abundances <- unclass(abundances)
# head(abundances)
#
# # Attaching some column metadata.
# extra.info <- colData(merged)[match(colnames(abundances), merged$sample),]
#
# y.ab <- DGEList(abundances, samples=extra.info)
# y.ab
#
# # Remove cell groups with few cells
# keep <- filterByExpr(y.ab, group=y.ab$samples$tomato)
# y.ab <- y.ab[keep,]
# summary(keep)
#
# # Model formula
# design <- model.matrix(~factor(pool) + factor(tomato), y.ab$samples)
#
# y.ab <- estimateDisp(y.ab, design, trend="none")
# summary(y.ab$common.dispersion)
#
# plotBCV(y.ab, cex=1)
#
# fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
# summary(fit.ab$var.prior)
#
# summary(fit.ab$df.prior)
#
# plotQLDisp(fit.ab, cex=1)
#
# res <- glmQLFTest(fit.ab, coef=ncol(design))
# summary(decideTests(res))


# Relevant example ----

# dat <- test_vdj %>%
#   calc_frequency("seurat_clusters", "sample", return_df = TRUE) %>%
#   distinct(sample, orig.ident, seurat_clusters, seurat_clusters_freq)
#
# extra.info <- dat %>%
#   distinct(orig.ident, sample) %>%
#   mutate(pool = str_extract(sample, "[0-9]+$"))
#
# abnd <- dat %>%
#   dplyr::select(-orig.ident) %>%
#   pivot_wider(names_from = "sample", values_from = "seurat_clusters_freq", values_fill = 0) %>%
#   column_to_rownames("seurat_clusters") %>%
#   as.matrix()
#
# y.ab <- DGEList(abnd, samples = extra.info)
#
# # Remove cell groups with few cells
# keep <- filterByExpr(y.ab, group = y.ab$samples$orig.ident)
# y.ab <- y.ab[keep, ]
#
# summary(keep)
#
# # Model formula
# design <- model.matrix(~ factor(pool) + factor(orig.ident), y.ab$samples)
#
# y.ab <- estimateDisp(y.ab, design, trend = "none")
# summary(y.ab$common.dispersion)
#
# plotBCV(y.ab, cex = 1)
#
# fit.ab <- glmQLFit(y.ab, design, robust = TRUE, abundance.trend = FALSE)
# summary(fit.ab$var.prior)
#
# summary(fit.ab$df.prior)
#
# plotQLDisp(fit.ab, cex=1)
#
# res <- glmQLFTest(fit.ab, coef = ncol(design))
#
# summary(decideTests(res))
