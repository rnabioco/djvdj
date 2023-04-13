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
                         method = NULL, file = NULL) {

  # Set method based on number of clusters for comparison
  if (is.null(method)) {
    n_clsts <- n_distinct(df_in[[cluster_col]])

    method <- dplyr::case_when(
      n_clsts == 1 ~ "none",
      n_clsts > 2  ~ "kruskal",
      TRUE         ~ "t"
    )
  }

  p_fn <- switch(
    method,
    t       = function(x, ...) .calc_p(x, ..., method = stats::t.test),
    wilcox  = function(x, ...) .calc_p(x, ..., method = stats::wilcox.test),
    kruskal = function(x, ...) (stats::kruskal.test(x, ...))$p.value,
    none    = function(x, ...) as.numeric(NA)
  )

  # Calculate p-values
  p_clmn <- paste0(method, "_pval")

  df_in <- dplyr::group_by(df_in, !!sym(group_col))

  res <- dplyr::mutate(
    df_in,
    !!sym(p_clmn) := p_fn(!!sym(data_col), !!sym(cluster_col))
  )

  res <- dplyr::ungroup(res)

  # Multiple testing correction
  p_res <- dplyr::distinct(res, !!sym(group_col), !!sym(p_clmn))

  p_res <- dplyr::mutate(
    p_res, p_adj = stats::p.adjust(!!sym(p_clmn), method = "bonferroni")
  )

  p_res <- purrr::set_names(p_res$p_adj, p_res[[group_col]])

  res <- dplyr::mutate(res, p_adj = p_res[as.character(!!sym(group_col))])

  # Write table
  if (!is.null(file)) {
    p_tbl <- dplyr::arrange(res, p_adj, !!!syms(c(group_col, cluster_col)))
    # p_tbl <- dplyr::relocate(p_tbl, !!!syms(c(group_col, cluster_col)))

    readr::write_csv(p_tbl, file, progress = FALSE)
  }

  res
}

.calc_p <- function(x, g, method) {
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

#' Format p-values plotting
#' @noRd
.format_pvalue <- function(p, digits = 1, cutoffs = NULL) {

  # Set p label based on vector of cutoffs
  if (!is.finite(p)) return(as.character(NA))

  if (!is.null(cutoffs)) {
    if (any(duplicated(cutoffs))) {
      cli::cli_abort("If p_label is a vector each cutoff value must be unique.")
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

    return(p_label)
  }

  # Format p-value label
  # if (p >= 0.01) return(as.character(round(p, digits = 2)))

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
#
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
#
#
# # Relevant example ----
#
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
# res <- glmQLFTest(fit.ab, coef=ncol(design))
# summary(decideTests(res))
