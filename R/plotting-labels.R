#' Calculate p-value for provided data
#'
#' @param df_in data.frame
#' @param data_col Column containing data to use for test
#' @param cluster_col Column containing clusters to compare
#' @param group_col Column containing labels to use for grouping data, a
#' p-value will be calculated for each group. Each group should have at least
#' two clusters from `cluster_col`
#' @param p_method Method to use for calculating p-value, possible values
#' include:
#'
#' * 't', two sample t-test
#' * 'wilcox', Wilcoxon rank sum test
#' * 'kruskal', Kruskal-Wallis rank sum test
#' * 'edgeR', test for differential abundance using edgeR
#'
#' @param adj_method Method to use for multiple testing correction
#' @param file File path to use for saving p-value csv
#' @param ... Not used, absorbs extra arguments
#' @noRd
.calc_pvalue <- function(df_in, data_col, cluster_col, group_col = NULL,
                         p_method = NULL, adj_method = "bonferroni",
                         file = NULL) {

  # Set method based on number of clusters for comparison
  n_clsts <- n_distinct(df_in[[cluster_col]])

  if (n_clsts > 2 && !is.null(p_method) && !identical(p_method, "kruskal")) {
    cli::cli_warn(
      "p-values will be calculated using the Kruskal-Wallis test
       since more than two groups are being compared."
    )
  }

  p_method <- dplyr::case_when(
    n_clsts == 1 ~ "none",
    n_clsts > 2  ~ "kruskal",
    TRUE         ~ p_method %||% "t"
  )

  # Calculate p-values
  if (identical(p_method, "edgeR")) {
    p_fn <- .calc_edgeR

    adj_method <- NULL

  } else {
    p_fn <- .calc_p
  }

  res <- p_fn(
    df_in, data_col = data_col,
    cluster_col = cluster_col, group_col = group_col, p_method = p_method
  )

  # Multiple testing correction
  p_clmn <- "p_value"

  if (!is.null(adj_method)) {
    p_res <- dplyr::distinct(res, !!!syms(c(group_col, p_clmn)))

    p_res <- dplyr::mutate(
      p_res, p_adj = stats::p.adjust(!!sym(p_clmn), method = adj_method)
    )

    res <- dplyr::left_join(res, p_res, by = c("p_value", group_col))

    p_clmn <- "p_adj"
  }

  # Write table
  if (!is.null(file)) {
    p_clmns <- unique(c(group_col, cluster_col, data_col, "p_value", p_clmn))

    p_tbl <- dplyr::arrange(res, !!!syms(c(p_clmn, group_col, cluster_col)))
    p_tbl <- dplyr::mutate(p_tbl, method = p_method)
    p_tbl <- dplyr::select(p_tbl, all_of(p_clmns), .data$method)

    readr::write_csv(p_tbl, file, progress = FALSE)
  }

  res <- dplyr::mutate(res, p_value = !!sym(p_clmn))

  res
}

.calc_p <- function(df_in, data_col, cluster_col, group_col = NULL,
                    p_method = "t") {

  # Helper to perform selected test
  .calc_test <- function(x, g, p_method) {
    p <- split(x, g, drop = TRUE)

    if (length(p) != 2) {
      cli::cli_abort(
        "There must be two groups present to perform test.", .internal = TRUE
      )
    }

    p <- p_method(p[[1]], p[[2]])
    p <- p$p.value

    p
  }

  # Return NAs if p_method is "none"
  if (identical(p_method, "none")) {
    res <- dplyr::mutate(df_in, p_value = as.numeric(NA))

    return(res)
  }

  # Set testing p_method
  p_fn <- switch(
    p_method,
    t       = function(x, ...) .calc_test(x, ..., p_method = stats::t.test),
    wilcox  = function(x, ...) .calc_test(x, ..., p_method = stats::wilcox.test),
    kruskal = function(x, ...) (stats::kruskal.test(x, ...))$p.value,
  )

  # Calculate p-values
  if (!is.null(group_col)) df_in <- dplyr::group_by(df_in, !!sym(group_col))

  res <- dplyr::mutate(
    df_in, p_value = p_fn(!!sym(data_col), !!sym(cluster_col))
  )

  res <- dplyr::ungroup(res)

  if (any(is.na(res$p_value))) {
    cli::cli_warn("`NA`s were returned when calculating p-values.")
  }

  res
}

.calc_edgeR <- function(df_in, data_col, cluster_col, group_col, ...) {

  .check_packages("edgeR", db = "Bioconductor")

  # Format input data
  dat <- dplyr::select(df_in, all_of(c(data_col, cluster_col, group_col)))
  dat <- dplyr::filter(dat, !is.na(as.character(!!sym(group_col))))

  dat <- dplyr::group_by(dat, !!!syms(c(group_col, cluster_col)))

  dat <- dplyr::mutate(
    dat,
    .pool   = row_number(),
    .sample = paste0(!!sym(cluster_col), "_", .data$.pool)
  )

  dat <- dplyr::ungroup(dat)

  # data.frame of sample info
  sam_info <- dplyr::distinct(
    dat, !!sym(cluster_col), .data$.sample, .data$.pool
  )

  # Counts matrix
  mat <- dplyr::select(dat, all_of(c(group_col, ".sample", data_col)))

  mat <- tidyr::pivot_wider(
    mat, names_from = .data$.sample,
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
  de_fit <- edgeR::glmQLFit(
    de_l, design, robust = TRUE, abundance.trend = FALSE
  )

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
  if (p >= 0.1) return(as.character(round(p, 1)))

  p <- scales::label_scientific(digits = digits)(p)

  ex <- .str_extract_all(p, "[+\\-][0-9]+$")

  p <- sub(paste0("\\", ex, "$"), "", p)

  ex <- as.numeric(ex)
  ex <- as.character(ex)

  p <- sub("e", "*x*10^", p)
  p <- paste0(p, ex)

  p
}


#' Add n label to plot
#'
#' @param gg_in ggplot2 object
#' @param df_in data.frame to use for counting number of values plotted. For
#' `.add_n_label()` can also be named list providing separate data.frames for
#' corner, axis, and legend labels.
#' @param grp Variable to use for grouping data when counting the number of
#' values
#' @param n_label Vector indicating where n labels should be added
#' @param crnr_col Column in `df_in` containing groups that will be shown for
#' the corner label
#' @param axis_col Column in `df_in` containing groups that will be shown on the
#' x-axis
#' @param lgnd_col Column in `df_in` containing groups that will be shown in
#' the legend
#' @param lgnd_clrs Colors to pass to `ggplot2::scale_fill_manual()` and
#' `ggplot2::scale_color_manual()`, only applicable when adding legend label
#' @param na_clr Color to use for `NA` value, only applicable when adding legend
#' label
#' @param y_exp y-axis expansion, only applicable when adding corner label
#' @param n_fn Function to use for calculating number of values plotted. By
#' default this is `dplyr::n` which will just count the number of rows for
#' each group. If another function is provided, it should take a vector as
#' input. The function will be applied to the `n_col` column in `df_in`. e.g.
#' `sum` will sum the values in `n_col` for each group, this is useful if the
#' number of cells has already been counted for each group.
#' @param n_col Column to store n values, this also specifies the column that
#' should be modified when a function is provided to the `n_fn`.
#' @param lab_args named list with aesthetic parameters to used for modifying
#' n label
#' @param axis Should label be added to the x- or y-axis, only applicable when
#' adding axis label
#' @param ... Absorbs unused arguments passed to label functions
#' @param sep Separator to use when creating n label
#' @return ggplot object with n labels added
#' @noRd
.add_n_label <- function(gg_in, df_in, n_label, crnr_col = NULL,
                         axis_col = NULL, lgnd_col = NULL, lgnd_clrs = NULL,
                         na_clr = "grey80", y_exp = .n_label_expansion,
                         n_fn = dplyr::n, lab_args = list()) {

  n_label <- unique(c("none", n_label))

  lab_args <- .parse_label_params(lab_args)$n

  if (!is.data.frame(df_in) && is.list(df_in)) {
    crnr_dat <- df_in$corner
    axis_dat <- df_in$axis
    lgnd_dat <- df_in$legend

  } else {
    crnr_dat <- axis_dat <- lgnd_dat <- df_in
  }

  # Named list containing possible label functions, group columns, and data
  lab_fns <- list(
    none   = list(.add_no_label, NULL, NULL),
    corner = list(.add_corner_label, crnr_col, crnr_dat)
  )

  if (!is.null(axis_col) && !is.null(axis_dat)) {
    lab_fns$axis <- list(.add_axis_label, axis_col, axis_dat)
  }

  if (!is.null(lgnd_col) && !is.null(lgnd_dat)) {
    lab_fns$legend <- list(.add_legend_label, lgnd_col, lgnd_dat)
  }

  if (any(!n_label %in% names(lab_fns))) {
    cli::cli_abort("`n_label` can be any combination of {names(lab_fns)}")
  }

  lab_fns <- lab_fns[unique(n_label)]

  res <- gg_in

  for (fn in lab_fns) {
    f <- fn[[1]]
    g <- fn[[2]]
    d <- fn[[3]]

    dat <- .calc_n(df_in = d, grp = g, n_fn = n_fn)

    res <- f(
      res, dat,
      grp       = g,
      lgnd_clrs = lgnd_clrs,
      na_clr    = na_clr,
      y_exp     = y_exp,
      lab_args  = lab_args
    )
  }

  res
}

.add_corner_label <- function(gg_in, df_in, lab_args,
                              y_exp = .n_label_expansion, ...) {

  # Automatically adjust label position based on length of string
  dat <- .format_n_label(df_in)

  dat <- dplyr::mutate(
    dat,
    label = lab_args$label %||% .data$label,
    hjust = .get_label_just(.data$label)
  )

  lab_args$size <- lab_args$size %||% global$base_size * 0.8
  lab_args$size <- lab_args$size / ggplot2::.pt

  lab_args$mapping     <- ggplot2::aes(label = .data$label, hjust = .data$hjust)
  lab_args$data        <- dat
  lab_args$inherit.aes <- FALSE
  lab_args$x           <- lab_args$x %||% Inf
  lab_args$y           <- lab_args$y %||% Inf
  lab_args$vjust       <- lab_args$vjust %||% 1.5

  res <- gg_in +
    lift(ggplot2::geom_text)(lab_args)

  if (!is.null(y_exp)) {
    res <- res +
      ggplot2::scale_y_continuous(expand = y_exp)
  }

  res
}

.add_axis_label <- function(gg_in, df_in, grp, axis = "x", lab_args, ...) {

  if (is.null(grp)) return(gg_in)

  dat <- .format_n_label(df_in, grp)

  dat_labs <- purrr::set_names(dat$label, dat[[grp]])

  if (identical(axis, "x")) {
    res <- gg_in +
      ggplot2::scale_x_discrete(labels = dat_labs)

  } else if (identical(axis, "y")) {
    res <- gg_in +
      ggplot2::scale_y_discrete(labels = dat_labs)

  } else {
    cli::cli_abort("`axis` must be x or y")
  }

  res
}

.add_legend_label <- function(gg_in, df_in, grp, lab_args, lgnd_clrs = NULL,
                              na_clr = "grey80", ...) {

  if (is.null(grp)) return(gg_in)

  dat <- .format_n_label(df_in, grp)

  dat_labs <- purrr::set_names(dat$label, dat[[grp]])

  # Scale arguments
  gg_args <- list(labels = dat_labs)

  if (!is.null(lgnd_clrs)) {
    gg_args$values   <- lgnd_clrs
    gg_args$na.value <- na_clr

    res <- gg_in +
      lift(ggplot2::scale_color_manual)(gg_args) +
      lift(ggplot2::scale_fill_manual)(gg_args)

  } else {
    res <- gg_in +
      lift(ggplot2::scale_color_discrete)(gg_args) +
      lift(ggplot2::scale_fill_discrete)(gg_args)
  }

  res
}

.add_no_label <- function(gg_in, ...) gg_in

.format_n_label <- function(df_in, grp = NULL, sep = "\n", n_col = ".n") {
  res <- dplyr::mutate(
    df_in,
    label = scales::label_comma()(!!sym(n_col)),
    label = paste0("n = ", .data$label)
  )

  if (!is.null(grp)) {
    res <- dplyr::mutate(
      res,
      label = paste0(!!sym(grp), sep, .data$label)
    )
  }

  res
}

.calc_n <- function(df_in, grp = NULL, n_fn = dplyr::n, n_col = ".n") {
  res <- df_in

  if (is.null(df_in)) return(df_in)

  if (!is.null(grp)) res <- dplyr::group_by(res, !!sym(grp))

  if (identical(n_fn, dplyr::n)) {
    res <- dplyr::summarize(res, !!sym(n_col) := n_fn(), .groups = "drop")

  } else {
    res <- dplyr::summarize(
      res, !!sym(n_col) := n_fn(!!sym(n_col)), .groups = "drop"
    )
  }

  res
}

.n_label_expansion <- ggplot2::expansion(c(0.05, 0.1))


#' Parse label params
#' divide params based on prefix, e.g. 'n.' or 'p.'
#' this allows user to adjust aesthetics for a specific label
#' @noRd
.parse_label_params <- function(l, prfxs = c("n", "p", "clone"), sep = ".") {
  nms   <- names(l)
  prfxs <- purrr::set_names(paste0("^", prfxs, "\\", sep), prfxs)
  idxs  <- purrr::map(prfxs, ~ grepl(.x, nms))

  # Identify params that do not have any of the prefixes and
  # should be included in all
  shared <- purrr::map(idxs, `!`)       # do not have the prefix
  shared <- purrr::reduce(shared, `&`)  # do not have any of the prfxs
  shared <- l[shared]

  res <- purrr::imap(idxs, ~ {
    prfx <- .y
    prms <- l[.x]

    names(prms) <- purrr::map_chr(names(prms), ~ sub(prfxs[[prfx]], "", .x))

    prms <- append(prms, shared)
    prms <- .standardize_aes(prms)

    prms
  })

  res$other <- shared

  res
}


#' Get axis label based on axis units
#'
#' @param units Units for axis
#' @param sffx Suffix to include in label
#' @return Axis label
#' @noRd
.get_axis_label <- function(units, sfx = "cells") {
  res <- switch(
    units,
    frequency = paste0("number of ", sfx),
    percent   = paste0("% of ", sfx)
  )

  res
}

#' Set label justification based on length of label
#'
#' This allows labels of different lengths to be positioned the same distance
#' from the edge of the plot area
#'
#' @param x character vector of labels
#' @param base base justification, increase this value to increase spacing
#' @return numeric vector of justification values
#' @noRd
.get_label_just <- function(x, base = 0.5) {

  # Assumed height x width ratio
  char_h_w <- 1.5

  len <- nchar(x)
  res <- 1 + (1 / len * (base * char_h_w))

  res
}

#' Trim long labels
#'
#' @param x Character vector containing labels to trim
#' @param max_len Maximum number of characters to allow
#' @param ellipsis Ellipsis to add to indicate label has been trimmed
#' @noRd
trim_lab <- function(x, max_len = 25, ellipsis = "...") {
  len <- nchar(x)

  trim_me <- len > max_len

  x[trim_me] <- strtrim(x[trim_me], max_len)
  x[trim_me] <- paste0(x[trim_me], ellipsis)

  x
}

