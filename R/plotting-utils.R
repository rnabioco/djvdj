#' ggplot2 imports
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_histogram
#' @importFrom ggplot2 geom_density geom_tile geom_boxplot geom_violin geom_col
#' @importFrom ggplot2 scale_x_discrete scale_y_continuous scale_x_continuous
#' @importFrom ggplot2 position_dodge2 scale_color_manual scale_fill_manual
#' @importFrom ggplot2 scale_color_gradientn scale_fill_gradientn stat_summary
#' @importFrom ggplot2 facet_wrap guides guide_legend labs theme element_blank
#' @importFrom ggplot2 element_text element_line expansion after_stat
#' @noRd
NULL


#' Theme for djvdj plotting functions
#'
#' @param base_size base font size in pts
#' @param base_family base font family
#' @param base_line_size base size for line elements
#' @param base_rect_size base size for rect element
#' @param line_color color for line elements
#' @return ggplot theme
#' @importFrom ggplot2 %+replace%
#' @aliases %+replace%
#' @export %+replace%
#' @examples
#'
#' plot_scatter(vdj_so, data_col = "seurat_clusters") +
#'   djvdj_theme()
#'
#' @export
djvdj_theme <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22,
                        line_color = "grey85") {

  ggplot2::theme_classic(
    base_size      = base_size,
    base_family    = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    ggplot2::theme(
      strip.background  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(size = base_size),

      panel.border      = ggplot2::element_rect(fill = NA, color = line_color),
      panel.background  = ggplot2::element_blank(),

      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),

      axis.line         = ggplot2::element_blank(),
      axis.ticks        = ggplot2::element_line(color = line_color),

      complete = TRUE
    )
}


#' Create ggplot
#'
#' @param df_in data.frame
#' @param fn ggplot2 function to use for generating plot,
#' e.g. `ggplot2::geom_point`
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param grp Varible to use for `ggplot2::facet_wrap()`
#' @param .color Variable to use for color, or logical indicating whether `clrs`
#' should be used to set color, if `NULL`, `clrs` will be used
#' @param .fill Variable to use for fill, or logical indicating whether `clrs`
#' should be used to set fill, if `NULL`, `clrs` will be used
#' @param clrs Vector of colors for plotting
#' @param na_clr Color to use for `NA` values
#' @param trans_x Method to use for transforming x-axis
#' @param trans_y Method to use for transforming y-axis
#' @param tranx_clr Method to use for transforming color
#' @param nrow Number of rows to use for arranging plot facets
#' @param scales Specification controlling facet scales to pass to
#' `ggplot2::facet_wrap()`
#' @param n_label Vector indicating where n labels should be added
#' @param label_params Named list providing additional parameters to modify
#' n label aesthetics, e.g. list(size = 4, color = "red")
#' @param label_data data.frame to use for generating n label
#' @param ... Additional parameters to pass to `fn`
#' @return ggplot object
#' @noRd
.create_plot <- function(df_in, fn, x = NULL, y = NULL, grp = NULL,
                         .color = NULL, .fill = NULL, clrs = NULL,
                         na_clr = "grey80", trans_x = "identity",
                         trans_y = "identity", trans_clr = "identity",
                         nrow = NULL, scales = "fixed", n_label = NULL,
                         label_params = list(), n_fn = dplyr::n,
                         label_data = df_in, ...) {

  # Remove rows with missing values
  # Need to do this so n labels are correct
  # Only remove NAs for numeric columns since character or factor NAs can be
  # plotted
  chk_na <- purrr::set_names(c(x, y))
  chk_na <- map_lgl(chk_na, ~ is.numeric(df_in[[.x]]))
  chk_na <- names(chk_na[chk_na])
  n_orig <- nrow(df_in)

  df_in <- dplyr::filter(df_in, dplyr::if_all(all_of(chk_na), ~ !is.na(.x)))

  n_rm <- n_orig - nrow(df_in)

  if (n_rm > 0) {
    cli::cli_warn("Removed {n_rm} row{?s} containing missing value{?s}")
  }

  # Check inputs
  scale_clr  <- is.character(.color)
  scale_fill <- is.character(.fill)

  if (!scale_clr && !scale_fill) clrs <- clrs[1]

  if (scale_clr && scale_fill && !identical(.fill, .color)) {
    cli::cli_abort(
      "If both `.color` and `.fill` are provided,
       they must refer to the same variable", .internal = TRUE
    )
  }

  # Set aesthetics and geom arguments
  num_clr  <- scale_clr  && is.numeric(df_in[[.color]])
  num_fill <- scale_fill && is.numeric(df_in[[.fill]])

  gg_aes  <- ggplot2::aes()
  gg_args <- .standardize_aes(list(...))

  if (is.null(x)) gg_aes$x <- y %||% "x"
  else            gg_aes$x <- sym(x)

  if (!is.null(y)) gg_aes$y      <- sym(y)
  if (scale_fill)  gg_aes$fill   <- sym(.fill)
  if (scale_clr)   gg_aes$colour <- sym(.color)

  if (!scale_fill && !scale_clr) {
    .color <- is.null(.color) || (is.logical(.color) && .color)
    .fill  <- is.null(.fill)  || (is.logical(.fill) && .fill)

    if (.color) gg_args$colour <- gg_args$colour %||% clrs
    if (.fill)  gg_args$fill   <- gg_args$fill   %||% clrs
  }

  # Allow alpha to be adjusted by passing new mapping
  # override alpha specification in gg_args otherwise
  # new alpha mapping will be ignored
  if ("alpha" %in% names(gg_args$mapping)) gg_args$alpha <- NULL

  # Create plot
  res <- ggplot2::ggplot(df_in, gg_aes) +
    lift(fn)(gg_args) +
    djvdj_theme()

  # Transform x-axis
  if (!identical(trans_x, "identity")) {
    res <- res +
      ggplot2::scale_x_continuous(trans = trans_x)
  }

  # Transform y-axis
  y_args <- list()

  if ("corner" %in% n_label)           y_args$expand <- .n_label_expansion
  if (!identical(trans_y, "identity")) y_args$trans  <- trans_y

  if (!purrr::is_empty(y_args)) {
    res <- res +
      lift(ggplot2::scale_y_continuous)(y_args)
  }

  # Set colors
  if (num_clr || num_fill) {
    clrs <- clrs %||% c("#132B43", "#56B1F7")

    if (num_clr) {
      res <- res +
        ggplot2::scale_color_gradientn(
          colors   = clrs,
          na.value = na_clr,
          trans    = trans_clr
        )
    }

    if (num_fill) {
      res <- res +
        ggplot2::scale_fill_gradientn(
          colors   = clrs,
          na.value = na_clr,
          trans    = trans_clr
        )
    }

  } else if (!is.null(clrs) && !"legend" %in% n_label) {
    res <- res +
      ggplot2::scale_fill_manual(values = clrs, na.value = na_clr) +
      ggplot2::scale_color_manual(values = clrs, na.value = na_clr)
  }

  # Create facets
  # Need to filter label data to only include labels for grps that are plotted
  if (!is.null(grp)) {
    res <- res +
      ggplot2::facet_wrap(
        stats::as.formula(paste("~", grp)),
        nrow = nrow, scales = scales
      )

    if (!is.data.frame(label_data) && is.list(label_data)) {
      label_data <- purrr::map(
        label_data, dplyr::filter, !!sym(grp) %in% df_in[[grp]]
      )

    } else {
      label_data <- dplyr::filter(label_data, !!sym(grp) %in% df_in[[grp]])
    }
  }

  # Add n label
  if (is.null(x)) n_label <- n_label[n_label != "axis"]

  .chk_num <- function(df_in, clmn) {
    if (!is.null(clmn) && !is.numeric(df_in[[clmn]])) return(clmn)
    else                                              return(NULL)
  }

  # can't use `%||%` here since .fill and .color can be logical
  if (scale_fill)     lgnd_col <- .fill
  else if (scale_clr) lgnd_col <- .color
  else                lgnd_col <- NULL

  res <- .add_n_label(
    res, label_data,
    n_label   = n_label,
    crnr_col  = grp,
    axis_col  = .chk_num(df_in, x),
    lgnd_col  = .chk_num(df_in, lgnd_col),
    lgnd_clrs = clrs,
    na_clr    = na_clr,
    n_fn      = n_fn,
    y_exp     = NULL,
    lab_args  = label_params
  )

  # Adjust theme
  if (is.null(x)) {
    res <- res +
      ggplot2::theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  res
}


#' Create a grouped plot summarizing replicates
#'
#' Currently this can only be used to plot frequency data
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param p_y Variable to use for calculating p-values
#' @param clst Variable containing cluster IDs, e.g. healthy-1, healthy-2,
#' disease-1, disease-2
#' @param grp Variable to use for grouping clusters IDs, e.g. healthy and
#' disease
#' @param method Method to use for generating plot, can be 'bar' or 'boxplot'
#' @param n_label n label specification
#' @param p_label Should p-values be shown on plot.
#' @param p_method Method to calculate p-values
#' @param p_grp Variable to use for grouping samples when calculating p-values.
#' A separate p-value will be calculated for each label in p_grp.
#' @param p_x Manually set x coordinate for p-label, provide x specification or
#' 'right', 'left', or 'center'.
#' @param p_file File path to save p-value csv
#' @param label_params Named list with specifications to modify label
#' aesthetics
#' @param show_points Should data points be shown on boxplot
#' @param add_zeros If `TRUE` zeros will be added for missing variables, this
#' should only be used when plotting frequency data. This is useful if missing
#' data should be plotted as a zero.
#' @param show_zeros If `TRUE` cell labels that are missing from a cluster will
#' still be shown on the plot.
#' @param ... Additional arguments to pass to plotting function
#' @return ggplot object
#' @noRd
.create_grouped_plot <- function(df_in, x, y, clst, grp, method = "bar",
                                 n_label = NULL, p_label = c(value = 0.05),
                                 p_y = y, p_method = NULL, p_grp = x, p_file = NULL,
                                 p_x = NULL, label_params = list(),
                                 show_points = TRUE, add_zeros = TRUE,
                                 show_zeros = TRUE, ...) {

  # Check arguments
  if (!is.numeric(p_label)) {
    .check_args(
      p_label = list(Class = "character", len_one = TRUE)
    )

    .check_possible_values(p_label = c("all", "none"))
  }

  .check_possible_values(method = c("bar", "boxplot"))

  # Add zeros for missing groups
  # this is only necessary when plotting frequency
  # if plotting another metric, e.g. diversity, this should be FALSE
  if (add_zeros) {
    df_in <- .add_missing_zeros(
      df_in,
      dat_cols   = c(y, p_y),
      expand_col = x,
      clst_col   = clst,
      grp_col    = grp
    )
  }

  # Calculate/format p-value labels
  # by default only show significant p-values
  if (identical(p_label, "all")) p_label <- c(value = Inf)

  add_p <- !identical(p_label, "none")

  if (add_p) {
    p <- .calc_pvalue(
      df_in,
      data_col = p_y, cluster_col = grp, group_col = p_grp,
      p_method = p_method, file = p_file
    )

    # Determine min and max y for adding p-value labels
    p <- dplyr::group_by(p, !!!syms(c(p_grp, "p_value")))

    p <- dplyr::summarize(
      p,
      y_min = min(!!sym(y)),
      y_max = max(!!sym(y)),
      gap   = (max(.data$y_max) - min(.data$y_min)) * 0.05,
      y     = .data$y_max + .data$gap,
      n_x   = dplyr::n_distinct(!!sym(x)),  # number of x-axis groups for each
      .groups = "drop"                      # p-value plotted, used for p_x
    )

    # Format p-values
    p <- dplyr::rowwise(p)

    p <- dplyr::mutate(
      p, p_lab = .format_pvalue(.data$p_value, cutoffs = p_label)
    )

    p <- dplyr::filter(p, !is.na(.data$p_lab))
    p <- dplyr::ungroup(p)

    if (nrow(p) == 0) add_p <- FALSE

    # Adjust p_x
    if (is.null(p_x) && is.null(p_grp)) p_x <- "center"

    if (!is.null(p_x)) {
      p <- dplyr::mutate(
        p,
        !!sym(x) := switch(
          as.character(p_x),  # EXPR should evaluate to character
          right  = Inf,
          left   = -Inf,
          center = (n_x / 2) + 0.5,
          p_x
        )
      )
    }
  }

  # Remove labels in grp that have all zeros
  # only should be used when plotting frequency data
  if (!show_zeros) {
    df_in <- dplyr::group_by(df_in, !!!syms(c(x, grp)))
    df_in <- dplyr::filter(df_in, !all(!!sym(y) == 0))
    df_in <- dplyr::ungroup(df_in)
  }

  # Plot arguments
  gg_args <- list(
    df_in        = df_in,
    x            = x,
    y            = y,
    n_label      = n_label,
    label_params = label_params,
    ...
  )

  # Create boxplots
  if (identical(method, "boxplot")) {
    gg_args$show_points <- show_points

    res <- lift(.create_boxes)(gg_args) +
      ggplot2::theme(legend.position = "right")

  # Create bar graphs
  } else {
    df_in <- dplyr::group_by(df_in, !!!syms(unique(c(x, grp, p_grp))))

    df_in <- dplyr::summarize(
      df_in,
      .sd       = stats::sd(!!sym(y)),
      !!sym(y) := mean(!!sym(y)),
      .groups   = "drop"
    )

    # Need to adjust y_max based on mean and sd to position p-value label
    if (add_p && is.null(p_x)) {
      y_dat <- dplyr::group_by(df_in, !!sym(p_grp))

      y_dat <- dplyr::mutate(
        y_dat,
        y_min = 0,
        y_max = !!sym(y) + .data$.sd,
        y_max = max(.data$y_max)
      )

      y_dat <- dplyr::ungroup(y_dat)
      y_dat <- dplyr::distinct(y_dat, !!sym(p_grp), .data$y_min, .data$y_max)

      p <- dplyr::select(p, -dplyr::all_of(c("y_max", "y_min")))
      p <- dplyr::left_join(p, y_dat, by = p_grp)
    }

    gg_args$df_in <- df_in
    gg_args$err   <- ".sd"

    gg_args$position <- gg_args$position %||%
      ggplot2::position_dodge2(preserve = "single", width = 0.8)

    res <- lift(.create_bars)(gg_args)
  }

  # Add p-values
  if (add_p) {
    all_sym <- !any(grepl("[a-zA-Z0-9]", names(p_label)))

    label_params <- .parse_label_params(label_params)$p

    label_params$mapping <- ggplot2::aes(
      x = !!sym(x), y = .data$y, label = .data$p_lab, fill = NULL
    )

    # If only one p-value calculated for plot, position center of panel
    # include 'p = ' when setting hjust for p label
    if (!is.null(p_x)) {
      label_params$mapping$y <- Inf
      label_params$vjust     <- label_params$vjust %||% 1.5

      p_just <- ifelse(all_sym, p$p_lab, paste0("p = ", p$p_lab))

      p <- dplyr::mutate(
        p,
        hjust = case_when(
          orig.ident == Inf  ~ .get_label_just(p_just),
          orig.ident == -Inf ~ .get_label_just(p_just, side = "left"),
          TRUE               ~ 0.5
        )
      )

      label_params$mapping$hjust <- sym("hjust")

      if (!all_sym) {
        p <- dplyr::mutate(p, p_lab = paste0("italic(p) == ", .data$p_lab))
      }
    }

    label_params$data   <- p
    label_params$parse  <- TRUE
    label_params$colour <- label_params$colour %||% "black"
    label_params$vjust  <- label_params$vjust  %||% 0

    # Set default label size
    # set larger font size for symbols, e.g. '*'
    if (is.null(label_params$size)) {
      label_params$size <- ifelse(
        all_sym,
        global$base_size * 1.3,
        global$base_size * 0.8
      )
    }

    label_params$size <- label_params$size / ggplot2::.pt

    res <- res +
      lift(ggplot2::geom_text)(label_params)

    if (!is.null(p_x) && !"corner" %in% n_label) {
      res <- res +
        ggplot2::scale_y_continuous(expand = .n_label_expansion)
    }
  }

  res
}


#' Create ggplot scatter plot
#'
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param .color Variable to use for color
#' @param n_label Vector indicating where n labels should be added
#' @param ... Additional arguments to pass to `.create_plot()`
#' @return ggplot object
#' @noRd
.create_scatter <- function(df_in, x, y, .color = NULL, n_label = NULL, ...) {

  # Check input values
  if (identical(x, y)) cli::cli_abort("`x` and `y` must be different")

  if (!is.numeric(df_in[[x]]) || !is.numeric(df_in[[y]])) {
    cli::cli_abort("`x` and `y` must refer to numeric columns")
  }

  # Set n label
  if (is.null(n_label)) {
    n_label <- "corner"

    if (!is.null(.color) && !is.numeric(df_in[[.color]])) {
      n_label <- c(n_label, "legend")
    }
  }

  # Create scatter plot
  res <- .create_plot(df_in, x, y, .color = .color, n_label = n_label, ...)

  res
}


#' Create ggplot bargraph
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param err Variable to use for error bars
#' @param trans Method to use for transforming y-axis
#' @param y_ttl Title for y-axis
#' @param x_ang Angle of x-axis text
#' @param x_hjst Horizontal justification for x-axis text
#' @param n_label Vector indicating where n labels should be added
#' @param ... Additional arguments to pass to `.create_plot()`
#' @return ggplot object
#' @noRd
.create_bars <- function(df_in, x = NULL, y, err = NULL, trans = "identity",
                         y_ttl = y, x_ang = NULL, x_hjst = NULL, n_label = NULL,
                         ...) {

  # Check input values
  if (!is.null(x) && is.numeric(df_in[[x]])) {
    cli::cli_abort("`x` cannot refer to a numeric column", call = NULL)
  }

  # Set n label
  if (is.null(n_label)) {
    if (is.null(x)) n_label <- "corner"
    else            n_label <- "axis"
  }

  # Create bar graph
  gg_args <- list(
    df_in   = df_in,
    fn      = ggplot2::geom_col,
    x       = x,
    y       = y,
    trans_y = trans,
    n_label = n_label,
    ...
  )

  gg_args$alpha <- gg_args$alpha %||% 0.5

  if (is.null(gg_args$position)) {
    gg_args$position <- ggplot2::position_dodge2(
      preserve = "single", padding = 0.05
    )
  }

  res <- lift(.create_plot)(gg_args)

  # Add error bars
  if (!is.null(err)) {
    gg_args <- .standardize_aes(gg_args)

    err_args <- list(
      position = ggplot2::position_dodge2(width = 0.9),
      show.legend = FALSE
    )

    err_args$colour    <- gg_args$colour
    err_args$linewidth <- gg_args$linewidth %||% 1

    err_args$mapping <- ggplot2::aes(
      !!sym(x), ymin = !!sym(y) - !!sym(err), ymax = !!sym(y) + !!sym(err)
    )

    res <- res +
      lift(ggplot2::geom_linerange)(err_args)
  }

  # Adjust theme
  if (!is.null(x) && dplyr::n_distinct(df_in[[x]]) > 6) {
    x_ang  <- x_ang  %||% 45
    x_hjst <- x_hjst %||% 1

  } else {
    x_ang  <- x_ang  %||% 0
    x_hjst <- x_hjst %||% 0.5
  }

  res <- res +
    ggplot2::labs(y = y_ttl) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  if (!is.null(x)) {
    res <- res +
      theme(axis.text.x = ggplot2::element_text(angle = x_ang, hjust = x_hjst))
  }

  res
}


#' Create ggplot boxplot
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param method Method to use for plotting, can be 'boxplot' or 'violin'
#' @param trans Method to use for transforming y-axis
#' @param y_ttl Title for y-axis
#' @param x_ang Angle of x-axis text
#' @param x_hjst Horizontal justification for x-axis text
#' @param show_points If `TRUE` data points will be shown on boxplots, or a
#' point indicating the median will be shown on violin plots
#' @param n_label Vector indicating where n labels should be added
#' @param ... Additional arguments to pass to `.create_plot()`
#' @return ggplot object
#' @noRd
.create_boxes <- function(df_in, x = NULL, y, method = "boxplot",
                          trans = "identity", y_ttl = y, x_ang = NULL,
                          x_hjst = NULL, n_label = NULL, show_points = TRUE,
                          point.size = NULL, point.color = NULL, ...) {

  # Set n label
  if (is.null(n_label)) {
    if (!is.null(x)) n_label <- "axis"
    else             n_label <- "corner"
  }

  # Set plotting function
  plt_fns <- list(
    boxplot = ggplot2::geom_boxplot,
    violin  = ggplot2::geom_violin
  )

  .check_possible_values(method = names(plt_fns))

  plt_fn <- plt_fns[[method]]

  gg_args <- list(
    df_in    = df_in,
    fn       = plt_fn,
    x        = x,
    y        = y,
    trans_y  = trans,
    n_label  = n_label,
    ...
  )

  gg_args <- .standardize_aes(gg_args)

  gg_args$alpha <- gg_args$alpha %||% 0.5

  if (identical(method, "violin")) pos <- ggplot2::position_dodge()
  else pos <- ggplot2::position_dodge2(preserve = "single")

  gg_args$position <- gg_args$position %||% pos

  # Create boxplots
  # allow user to set point size and color
  if (show_points && identical(method, "boxplot")) gg_args$outlier.color <- NA

  res <- lift(.create_plot)(gg_args)

  # Add additional points
  if (show_points) {
    pt_args        <- list()
    pt_args$alpha  <- 1
    pt_args$size   <- point.size %||% 1
    pt_args$colour <- point.color %||% gg_args$colour

    if (identical(method, "violin")) {
      pt_args$geom <- "point"
      pt_args$fun  <- stats::median

      res <- res +
        lift(ggplot2::stat_summary)(pt_args)

    } else {
      pt_args$position <- ggplot2::position_jitterdodge(jitter.width = 0.05)

      res <- res +
        lift(ggplot2::geom_jitter)(pt_args)
    }
  }

  # Adjust theme
  if (!is.null(x) && dplyr::n_distinct(df_in[[x]]) > 6) {
    x_ang  <- x_ang  %||% 45
    x_hjst <- x_hjst %||% 1

  } else {
    x_ang  <- x_ang  %||% 0
    x_hjst <- x_hjst %||% 0.5
  }

  res <- res +
    ggplot2::labs(y = y_ttl) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x    = ggplot2::element_blank(),
      axis.text.x     = ggplot2::element_text(angle = x_ang, hjust = x_hjst)
    )

  res
}


#' Create ggplot histogram
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param .color Variable to use for color
#' @param .fill Variable to use for fill
#' @param method Method to use for plotting, can be 'histogram' or 'density'
#' @param units Units to use for y-axis when plotting histogram. Use
#' 'frequency' to show the raw number of values or 'percent' to show
#' the percentage of total values.
#' @param trans Method to use for transforming x-axis
#' @param y_ttl Title for y-axis
#' @param n_label Vector indicating where n labels should be added
#' @param ... Additional arguments to pass to `.create_plot()`
#' @return ggplot object
#' @noRd
.create_hist <- function(df_in, x, .color = NULL, .fill = NULL,
                         method = "histogram", units = "frequency",
                         trans = "identity", y_ttl = units, n_label = NULL,
                         ...) {

  # Set plotting function
  plt_fns <- list(
    histogram = ggplot2::geom_histogram,
    density   = ggplot2::geom_density
  )

  # Check inputs
  .check_possible_values(
    method = names(plt_fns),
    units  = c("frequency", "percent")
  )

  # Set n label
  if (is.null(n_label)) {
    n_label <- "corner"

    n_label <- dplyr::case_when(
      !is.null(.color) && !is.numeric(df_in[[.color]]) ~ c(n_label, "legend"),
      !is.null(.fill)  && !is.numeric(df_in[[.fill]])  ~ c(n_label, "legend"),
      TRUE ~ n_label
    )
  }

  # Set aesthetics
  # Only plot percent for histogram
  gg_aes <- ggplot2::aes()

  if (identical(units, "percent") && identical(method, "histogram")) {
    gg_aes <- ggplot2::aes(y = ggplot2::after_stat((count / nrow(df_in)) * 100))
  }

  gg_aes$x <- sym(x)

  # Create histogram
  gg_args <- list(
    df_in   = df_in,
    fn      = plt_fns[[method]],
    mapping = gg_aes,
    x       = x,
    .color  = .color,
    .fill   = .fill,
    trans_x = trans,
    n_label = n_label,
    ...
  )

  gg_args$alpha <- gg_args$alpha %||% 0.5

  res <- lift(.create_plot)(gg_args)

  if (identical(method, "histogram")) {
    res <- res +
      ggplot2::labs(y = y_ttl)
  }

  res
}


#' Create ggplot heatmap
#'
#' @param df_in data.frame
#' @param x Variable to plot on the x-axis
#' @param y Variable to plot on the y-axis
#' @param .fill Variable to use for the fill color
#' @param clrs Vector of colors for plotting
#' @param na_color Color to use for missing values
#' @param trans Method to use for transforming data
#' @param plt_ttl Plot title
#' @param lgd_ttl Legend title
#' @param x_ang Angle of x-axis text
#' @param x_hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_gg_heatmap <- function(df_in, x = NULL, y, .fill, clrs = NULL,
                               plt_ttl = ggplot2::waiver(), lgd_ttl = .fill,
                               trans = "identity", x_ang = NULL, x_hjst = NULL,
                               ...) {

  # Set colors
  clrs <- clrs %||% "#619CFF"

  if (length(clrs) == 1) clrs <- c("white", clrs)

  # Create heatmap
  res <- .create_plot(
    df_in,
    fn        = ggplot2::geom_tile,
    x         = x,
    y         = y,
    .fill     = .fill,
    clrs      = clrs,
    trans_clr = trans,
    ...
  )

  # Adjust theme
  if (!is.null(x) && dplyr::n_distinct(df_in[[x]]) > 6) {
    x_ang  <- x_ang  %||% 45
    x_hjst <- x_hjst %||% 1

  } else {
    x_ang  <- x_ang  %||% 0
    x_hjst <- x_hjst %||% 0.5
  }

  res <- res +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.line   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = x_ang, hjust = x_hjst)
    ) +
    labs(title = plt_ttl, fill = lgd_ttl)

  res
}


#' Create ComplexHeatmap heatmap
#'
#' @param mat_in Matrix
#' @param clrs Vector of colors for plotting
#' @param na_color Color to use for missing values
#' @param lvls Vector specifying level order
#' @param lgd_ttl Legend title
#' @param rm_upper If TRUE, upper triangle for heatmap will not be shown and
#' rows/columns will not be clustered.
#' @param rm_diag If TRUE, diagonal for heatmap will not be shown and
#' rows/columns will not be clustered.
#' @param ... Aditional arguments to pass to ComplexHeatmap::Heatmap()
#' @importFrom grid grid.rect
#' @importFrom ComplexHeatmap Heatmap
#' @noRd
.create_heatmap <- function(mat_in, clrs = NULL, na_color = NA, lvls = NULL,
                            lgd_ttl = NULL, rm_upper = FALSE, rm_diag = FALSE,
                            cluster = TRUE, ...) {

  # Set plot levels
  plt_args <- list(...)

  plt_args$cell_fun <- NULL

  r_nms <- rownames(mat_in)
  c_nms <- colnames(mat_in)

  if (!is.null(lvls) || !cluster) {
    plt_args$cluster_rows    <- FALSE
    plt_args$cluster_columns <- FALSE
    plt_args$row_names_side  <- plt_args$row_names_side %||% "left"
  }

  if (is.null(lvls)) {
    lvls <- union(r_nms, c_nms)
    lvls <- sort(lvls)
  }

  r_lvls <- lvls[lvls %in% r_nms]
  c_lvls <- lvls[lvls %in% c_nms]
  mat_in <- mat_in[r_lvls, c_lvls]

  # Remove upper triangle and/or diagonal
  if (rm_upper || rm_diag) {
    mat_rm <- .remove_upper_triangle(mat_in, rm_upper, rm_diag)
    n_vals <- .get_unique_values(mat_rm)

    if (n_vals == 1) {
      cli::cli_warn(
        "Cannot remove diagonal since there will only be one unique value
         remaining"
      )
    } else {
      mat_in <- mat_rm

      plt_args$na_col <- plt_args$na_col %||% na_color
    }
  }

  # THIS WORKED PREVIOUSLY, BUT NOW DOES NOT WORK
  # using the 'cell_fun' argument seems like a cleaner way of removing upper
  # triangle than adding NAs to the matrix
  #
  # # Remove upper triangle
  # if (rm_upper) {
  #   cell_fun <- function(j, i, x, y, w, h, fill) {
  #     if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #       grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = fill, col = fill))
  #     }
  #   }
  #
  #   plt_args$rect_gp  <- grid::gpar(method = "none")
  #   plt_args$cell_fun <- cell_fun
  # }
  #
  # # Remove diagonal
  # if (rm_diag) {
  #   mat_rm <- .remove_upper_triangle(mat_in, rm_upper = FALSE, rm_diag)
  #   n_vals <- .get_unique_values(mat_rm)
  #
  #   if (n_vals == 1) {
  #     warning(
  #       "Cannot remove diagonal since there ",
  #       "will only be one unique value remaining."
  #     )
  #   } else {
  #     mat_in <- mat_rm
  #
  #     plt_args$na_col <- plt_args$na_col %||% na_color
  #   }
  # }

  # Set plot colors
  # if all values in matrix are the same, only pass one color, otherwise error
  clrs   <- clrs %||% "#619CFF"
  n_vals <- .get_unique_values(mat_in)

  if (n_vals == 1)            clrs <- dplyr::last(clrs)
  else if (length(clrs) == 1) clrs <- c("white", clrs)

  plt_args$col <- plt_args$col %||% clrs

  # Set final heatmap parameters
  # use computed similarities for clustering, so use a distance function that
  # just converts a matrix to a dist object
  # can't directly use as.dist as the function since it accepts too many
  # arguments
  dist_fn <- function(input) stats::as.dist(input)

  lgd_params <- list(
    title_gp      = grid::gpar(fontface = "plain"),
    legend_height = ggplot2::unit(80, "pt"),
    title         = plt_args$name %||% lgd_ttl
  )

  plt_args$matrix <- mat_in

  plt_args$heatmap_legend_param <- plt_args$heatmap_legend_param %||% lgd_params

  plt_args$clustering_distance_rows <-
    plt_args$clustering_distance_rows %||% dist_fn

  plt_args$clustering_distance_columns <-
    plt_args$clustering_distance_columns %||% dist_fn

  # Create heatmap
  res <- lift(ComplexHeatmap::Heatmap)(plt_args)

  res
}

#' Remove upper triangle and/or diagonal from matrix
#'
#' @param mat_in Matrix
#' @param rm_upper If TRUE, upper triangle will be removed from matrix
#' @param rm_diag If TRUE, diagonal will be removed from matrix
#' @noRd
.remove_upper_triangle <- function(mat_in, rm_upper, rm_diag) {

  # Check inputs
  if (!rm_upper && !rm_diag) return(mat_in)

  if (!isSymmetric(mat_in)) {
    cli::cli_abort(
      "Matrix must be symmetrical to remove upper triangle and/or diagonal"
    )
  }

  # Identify values to remove
  nms     <- colnames(mat_in)
  nms_key <- purrr::set_names(nms)

  nms_key <- purrr::imap(nms_key, ~ {
    idx <- grep(paste0("^", .x, "$"), nms)
    v   <- NULL

    if (rm_upper) {
      v <- nms[idx:length(nms)]
      v <- v[v != .y]
    }

    if (rm_diag) v <- c(v, .y)

    v
  })

  for (i in seq_along(nms_key)) {
    mat_in[names(nms_key[i]), nms_key[[i]]] <- NA
  }

  # Remove rows/columns with all NAs
  na_idx <- is.na(mat_in)
  r_idx  <- rowSums(na_idx) != ncol(mat_in)
  c_idx  <- colSums(na_idx) != nrow(mat_in)
  res    <- mat_in[r_idx, c_idx]

  res
}

#' Count number of unique values in matrix
#'
#' @param mat_in Matrix
#' @noRd
.get_unique_values <- function(mat_in) {
  vals <- as.numeric(mat_in)
  vals <- vals[is.finite(vals)]

  n_distinct(vals, na.rm = TRUE)
}


#' Create circos plot
#'
#' @param mat_in Matrix
#' @param clrs Vector of colors for plotting
#' @param na_color Color to use for missing values
#' @param lvls Vector specifying level order
#' @param grps Named vector specifying group for each column/row
#' @param rotate_labels Should labels be rotated to reduce overlapping text
#' @param plt_ttl Plot title
#' @param ... Additional arguments to pass to circlize::chordDiagram()
#' @importFrom scales hue_pal
#' @importFrom graphics title strwidth
#' @noRd
.create_circos <- function(mat_in, clrs = NULL, na_color = "grey90",
                           lvls = NULL, grps = NULL, rotate_labels = FALSE,
                           plt_ttl = NULL, ...) {

  # Check installed packages
  .check_packages("circlize")

  # Check matrix, must have at least one link
  vals <- as.numeric(mat_in)

  if (isSymmetric(mat_in)) {
    vals <- .remove_upper_triangle(mat_in, rm_upper = TRUE, rm_diag = TRUE)
  }

  if (!any(vals > 0, na.rm = TRUE)) {
    cli::cli_abort(
      "To create a circos plot there must be at least one overlap
       between the clusters."
    )
  }

  # Set levels
  all_nms <- union(colnames(mat_in), rownames(mat_in))

  if (!is.null(lvls) && !all(all_nms %in% lvls)) {
    cli::cli_abort(
      "`plot_levels` must include all rows and columns in the matrix"
    )
  }

  if (!is.null(grps)) {
    if (!all(all_nms %in% names(grps))) {
      cli::cli_abort("Groups must be provided for all clusters")
    }

    grps <- grps[lvls]
  }

  # Set plot colors
  # if clrs is NULL, use ggplot2 default palette
  clrs    <- clrs %||% scales::hue_pal()(length(all_nms))
  clr_lst <- .set_circos_cols(mat_in, clrs, na_color)
  clrs    <- clr_lst[[1]]
  c_clrs  <- clr_lst[[2]]

  # Set plot arguments
  plt_args <- list(...)

  plt_args$x          <- mat_in
  plt_args$symmetric  <- plt_args$symmetric %||% TRUE
  plt_args$grid.col   <- clrs
  plt_args$column.col <- c_clrs
  plt_args$order      <- lvls
  plt_args$group      <- plt_args$group %||% grps

  # If rotating labels need to exclude default axis track
  ann_trk <- "grid"

  if (!rotate_labels) ann_trk <- c("name", ann_trk)

  adj_axis <- is.null(plt_args[["annotationTrack"]])

  plt_args[["annotationTrack"]] <- plt_args[["annotationTrack"]] %||% ann_trk

  # Create circos plot
  # circlize::mm_h must be used within chordDiagram
  # strwidth must be used within chordDiagram to avoid
  # 'plot.new() has not been called yet' error
  circlize::circos.clear()

  pre_all  <- is.null(plt_args$preAllocateTracks) && rotate_labels
  track_ht <- is.null(plt_args$annotationTrackHeight)

  circos_fun <- function(...) circlize::chordDiagram(...)

  if (pre_all) {
    if (track_ht) {
      circos_fun <- function(...) {
        circlize::chordDiagram(
          preAllocateTracks = list(
            track.height = max(graphics::strwidth(unlist(dimnames(mat_in))))
          ),
          annotationTrackHeight = circlize::mm_h(c(3, 4)),
          ...
        )
      }
    } else {
      circos_fun <- function(...) {
        circlize::chordDiagram(
          preAllocateTracks = list(
            track.height = max(strwidth(unlist(dimnames(mat_in))))
          ),
          ...
        )
      }
    }
  } else if (track_ht) {
    circos_fun <- function(...) {
      circlize::chordDiagram(
        annotationTrackHeight = circlize::mm_h(c(3, 4)),
        ...
      )
    }
  }

  lift(circos_fun)(plt_args)

  # Add axis track
  if (adj_axis) {
    pan_fun <- function(x, y) {
      circlize::circos.axis(minor.ticks = 0, labels.cex = 0.5)
    }

    circlize::circos.track(track.index = 2, bg.border = NA, panel.fun = pan_fun)
  }

  # Add rotated labels
  if (rotate_labels) {
    circlize::circos.track(
      track.index = 1, bg.border = NA,
      panel.fun = function(x, y) {
        circlize::circos.text(
          circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1],
          circlize::CELL_META$sector.index, facing = "clockwise",
          niceFacing = TRUE, adj = c(-0.3, 0.5), cex = 0.8
        )
      }
    )
  }

  # Add title
  graphics::title(main = plt_ttl, font.main = 1)
}

#' Set colors for circos plot
#'
#' @param mat_in Matrix
#' @param clrs Colors
#' @return List containing a vector of sector colors and a vector of link
#' colors
#' @noRd
.set_circos_cols <- function(mat_in, clrs, na_color = "grey90") {
  c_nms   <- colnames(mat_in)
  r_nms   <- rownames(mat_in)
  r_nms   <- r_nms[!r_nms %in% c_nms]
  all_nms <- union(c_nms, r_nms)
  n_clrs  <- length(clrs)

  if (is.null(names(clrs))) {
    if (n_clrs >= length(all_nms)) {
      clrs   <- clrs[seq_along(all_nms)]
      clrs   <- purrr::set_names(clrs, all_nms)
      c_clrs <- clrs[c_nms]

    } else if (n_clrs >= length(c_nms)) {
      clrs   <- clrs[seq_along(c_nms)]
      c_clrs <- purrr::set_names(clrs, c_nms)
      r_clrs <- purrr::set_names(rep(na_color, length(r_nms)), r_nms)

      clrs   <- c(c_clrs, r_clrs)
      c_clrs <- clrs[c_nms]

    } else if (n_clrs == 1) {
      c_clrs <- clrs

    } else {
      cli::cli_abort(
        "Not enough colors provided ({n_clrs}), provide exactly one
         color, one color for each column in the matrix ({length(c_nms)}),
         or one color for each unique row and column in the matrix
         ({length(all_nms)})"
      )
    }

  } else {
    if (!all(all_nms %in% names(clrs))) {
      cli::cli_abort("A color must be provided for each group being plotted")
    }

    c_clrs <- clrs[c_nms]
  }

  res <- list(clrs, c_clrs)

  res
}

#' Set 'other' groups
#'
#' Label groups that are not among the most frequent as 'other'
#'
#' @param df_in data.frame
#' @param data_col Column containing groups/clusters to rank
#' @param val_col Column containing values for ranking groups/clusters
#' @param method Method to use for ranking groups, possible values include:
#'
#' - A function to use for ranking values in `data_col`, this should take a
#'   single vector as input and will be used to summarize the values in
#'   `val_col` (if `val_col` is not `NULL`)
#' - 'boxplot', rank values in `data_col` based on boxplot stats calculated
#'   for `val_col`
#'
#' @param rev If `TRUE` reverse order of ranked groups so smallest values are
#' shown first
#' @param top Top groups to include, other groups will be labeled as
#' 'other'. Should be integer or character vector.
#' @param other_label Label to use for 'other' groups
#' @noRd
.rank_values <- function(df_in, data_col, val_col = NULL, method = n,
                         rev = FALSE) {

  res <- dplyr::group_by(df_in, !!sym(data_col))

  if (!is.null(val_col)) {
    res <- dplyr::summarize(res, met = method(!!sym(val_col)))

  } else {
    res <- dplyr::summarize(res, met = method())
  }

  if (rev) res <- dplyr::arrange(res, .data$met, !!sym(data_col))
  else     res <- dplyr::arrange(res, dplyr::desc(.data$met), !!sym(data_col))

  res <- dplyr::pull(res, data_col)
  res <- stats::na.omit(res)
  res
}

.rank_boxplot <- function(df_in, data_col, val_col, rev = FALSE, ...) {

  res <- dplyr::group_by(df_in, !!sym(data_col))

  res <- dplyr::summarize(
    res,
    met = list(grDevices::boxplot.stats(
      stats::na.omit(!!sym(val_col)),
      do.conf = FALSE, do.out = FALSE
    )),
    med = median(!!sym(val_col), na.rm = TRUE),
    max = max(!!sym(val_col), na.rm = TRUE),
    .groups = "drop"
  )

  res <- dplyr::mutate(
    res,
    q3 = map_dbl(.data$met, ~ .x$stats[4]),
    q4 = map_dbl(.data$met, ~ .x$stats[5])
  )

  if (rev) {
    res <- dplyr::arrange(
      res,
      .data$med, .data$q3, .data$q4, .data$max, !!sym(data_col)
    )

  } else {
    res <- dplyr::arrange(
      res,
      dplyr::desc(.data$med), dplyr::desc(.data$q3),
      dplyr::desc(.data$q4),  dplyr::desc(.data$max),
      !!sym(data_col)
    )
  }

  res <- dplyr::pull(res, data_col)
  res <- stats::na.omit(res)
  res
}

.set_other_grps <- function(df_in, data_col, val_col = NULL, plot_lvls = NULL,
                            top = NULL, other_label = "other", rev = FALSE,
                            method = n) {

  # Set ranking function
  if (identical(method, "boxplot")) rank_fn <- .rank_boxplot
  else                              rank_fn <- .rank_values

  # If numeric or NULL, do nothing, return df_in
  if (is.null(data_col) || is.numeric(df_in[[data_col]])) return(df_in)

  # Convert logical or factor to character
  dat <- df_in[[data_col]]

  if (is.logical(dat) || is.factor(dat)) {
    df_in <- purrr::modify_at(df_in, data_col, as.character)
  }

  # Rank data_col
  if (is.character(top)) {
    top_grps <- unique(top)
    n_top    <- length(top_grps)

    rnk    <- sort(unique(df_in[[data_col]]))
    rnk    <- c(top_grps, rnk[!rnk %in% top_grps])
    n_grps <- length(rnk)

  } else {
    rnk    <- rank_fn(df_in, data_col, val_col, method = method)
    n_grps <- length(rnk)

    top      <- top[1]
    n_top    <- top %||% ifelse(n_grps > 50, 10, 20)
    top_grps <- head(rnk, n_top)

    # Warn user when top is automatically set
    if (is.null(top) && n_grps > n_top) {
      cli::cli_warn(
        "The top {n_top} groups are shown, the remaining data points are
         labeled as \'{other_label}\'. This behavior can be modified using the
         {.code top} argument."
      )
    }
  }

  res <- df_in

  if (n_top < n_grps && !is.null(other_label)) {
    res <- dplyr::mutate(df_in, !!sym(data_col) := ifelse(
      !!sym(data_col) %in% top_grps | is.na(!!sym(data_col)),
      !!sym(data_col),
      other_label
    ))

    rnk <- c(top_grps, other_label)
  }

  # Always order values in plot_lvls first
  # If all levels are provided, use user provided ordering, do not reverse
  # ACTUALLY MAKES MORE SENSE TO ALWAYS REVERSE
  # remove levels that do not appear in data
  # If user provides NA in plot_levels, keep in order, otherwise always order
  # NA last
  not_in    <- rnk[!rnk %in% plot_lvls]
  plot_lvls <- plot_lvls[plot_lvls %in% c(rnk, NA)]
  lvls      <- c(plot_lvls, not_in)

  if (!NA %in% lvls) lvls <- c(lvls, NA)

  # if (rev && !purrr::is_empty(not_in)) lvls <- rev(lvls)
  if (rev) lvls <- rev(lvls)

  res <- .set_lvls(res, data_col, lvls)

  res
}

#' Add zeros for missing groups
#'
#' This adds zeros for combinations of groups that do not have any cells.
#' This is allows missing groups to be plotted as zeros.
#'
#' @param df_in data.frame or tibble
#' @param dat_col Columns containing data to add zeros
#' @param expand_col Column containing variable to expand for each group
#' @param clst_col Column containing cluster IDs, such as patient IDs,
#' e.g. healthy-1, healthy-2, disease-1, disease-2.
#' @param grp_col Column containing variable to use for grouping cluster IDs,
#' e.g. health and disease.
#' data.frame is expanded so each group has all variables in expand_col and
#' clst_col. However, clst_col is expanded for each group separately since
#' groups are expected to have different values for clst_col.
#' @noRd
.add_missing_zeros <- function(df_in, dat_cols, expand_col, clst_col, grp_col) {

  # If expand_col and grp_col are the same, no need to add zeros since each
  # group contains a single value
  if (identical(expand_col, grp_col)) return(df_in)

  # Save original levels to combinations to maintain ordering
  # convert all factors to character to allow correct merging at the end
  lvls <- purrr::map(df_in, levels)
  lvls <- purrr::discard(lvls, is.null)

  all_cols <- unique(c(expand_col, clst_col, grp_col))

  df_in <- dplyr::mutate(
    df_in,
    dplyr::across(all_of(all_cols), as.character)
  )

  # Get all combinations
  # make expand_col factor so all groups get all values of expand_col
  all <- dplyr::mutate(df_in, !!sym(expand_col) := factor(!!sym(expand_col)))
  all <- split(all, all[[grp_col]])

  all <- purrr::map(
    all, tidyr::expand, !!!syms(all_cols)
  )

  all <- dplyr::bind_rows(all)

  # Convert expand_col back to character to allow for correct merging of NAs
  # which might be present when plotting VDJ data
  all <- dplyr::mutate(
    all, !!sym(expand_col) := as.character(!!sym(expand_col))
  )

  # Add zeros for combinations that are expected but have no cells
  res <- dplyr::right_join(df_in, all, by = all_cols)

  res <- dplyr::mutate(
    res, dplyr::across(all_of(unique(dat_cols)), ~ tidyr::replace_na(.x, 0))
  )

  res <- dplyr::mutate(res, dplyr::across(
    all_of(names(lvls)),
    ~ factor(.x, levels = lvls[[dplyr::cur_column()]], exclude = NULL)
  ))

  res
}

#' Get arguments used solely by provided text function
#'
#' This is important for making sure `label_params` arguments are used for the
#' correct labels
#' i.e. when passing arguments to geom_text, want to exclude any parameters used
#' solely by geom_text_repel
#'
#' @param args_list List of arguments to filter based on `fn`
#' @param fn Function to filter `args_list`
#' @noRd
.get_uniq_text_args <- function(args_list, fn) {

  .get_unique_args <- function(...) {
    args        <- .standardize_aes(list(...))
    args        <- purrr::map(args, unique)
    shared_args <- purrr::reduce(args, dplyr::intersect)
    args        <- purrr::map(args, ~ .x[!.x %in% shared_args])

    args
  }

  # Get args unique for each function
  uniq_args <- .get_unique_args(
    geom_text = c(
      formalArgs(ggplot2::geom_text),
      names(ggplot2::GeomText$default_aes)
    ),
    geom_text_repel = c(
      formalArgs(ggrepel::geom_text_repel),
      names(ggrepel::GeomTextRepel$default_aes)
    )
  )

  if (!fn %in% names(uniq_args)) {
    cli::cli_abort("`fn` must be {.or {names(uniq_args)}}", .internal = TRUE)
  }

  # Get args not used by fn
  uniq_args <- uniq_args[names(uniq_args) != fn]
  uniq_args <- purrr::reduce(uniq_args, c)

  # Remove args from args_list that are used uniquely by other function
  args <- names(args_list)
  args <- args[!args %in% uniq_args]

  res <- args_list[args]

  res
}


#' Set min and max values for column
#'
#' @param df_in Input data.frame
#' @param ft Name of column containing feature values
#' @param mn Minimum quantile
#' @param mx Maximum quantile
#' @return data.frame with modified feature values
#' @noRd
.set_lims <- function(df_in, ft, mn = NULL, mx = NULL) {

  if (!is.numeric(df_in[[ft]])) {
    cli::cli_abort("The {ft} column is not numeric")
  }

  if (is.null(mn) && is.null(mx)) {
    cli::cli_abort("`mn` or `mx` must be provided")
  }

  ft <- sym(ft)

  res <- mutate(df_in, pct = dplyr::percent_rank(!!ft))

  if (!is.null(mn)) {
    mn <- mn / 100

    res <- dplyr::mutate(
      res,
      !!ft := ifelse(.data$pct > mn, !!ft, NA),
      !!ft := ifelse(.data$pct <= mn, min(!!ft, na.rm = TRUE), !!ft)
    )
  }

  if (!is.null(mx)) {
    mx <- mx / 100

    res <- dplyr::mutate(
      res,
      !!ft := ifelse(.data$pct < mx, !!ft, NA),
      !!ft := ifelse(.data$pct >= mx, max(!!ft, na.rm = TRUE), !!ft)
    )
  }

  res <- dplyr::select(res, -"pct")

  res
}


#' Order column levels by another variable
#'
#' @param x Vector to modify
#' @param by Vector of values to sort by
#' @param fun Summary function
#' @param decreasing Should sort be decreasing
#' @return Factor
#' @noRd
.order_lvls <- function(x, by, fun = max, decreasing = TRUE) {

  df <- tibble::tibble(x = as.character(x), by = by)

  l <- dplyr::group_by(df, x)
  l <- dplyr::summarize(l, by = fun(by), .groups = "drop")

  l <- dplyr::arrange(l, by)
  l <- l$x

  if (decreasing) l <- rev(l)

  res <- factor(x, l, exclude = NULL)

  res
}

#' Set column levels
#'
#' @param df_in data.frame
#' @param clmn Column to modify
#' @param lvls Levels
#' @return data.frame with modified clmn levels
#' @noRd
.set_lvls <- function(df_in, clmn, lvls) {

  if (is.null(lvls) || is.null(clmn)) return(df_in)

  dat <- df_in[[clmn]]

  if (!is.character(dat) && !is.factor(dat)) {
    cli::cli_warn(
      "Plot levels can only be modified for characters and factors, levels
       were not modified"
    )

    return(df_in)
  }

  # is.na() will not detect NAs in factor
  if (any(is.na(as.character(dat))) && !any(is.na(lvls))) lvls <- c(NA, lvls)

  u_dat   <- unique(dat)
  missing <- u_dat[!u_dat %in% lvls]

  if (!purrr::is_empty(missing)) {
    cli::cli_abort("Some labels in {clmn} are not in `plot_lvls`: {missing}")
  }

  df_in[clmn] <- factor(df_in[[clmn]], unique(lvls), exclude = NULL)

  df_in
}

#' Set default colors, add levels as names
#'
#' @param df_in data.frame
#' @param data_col Data column
#' @param plot_colors Colors
#' @param plot_lvls Levels
#' @param other_label Label to use for 'other' group
#' @param other_color Color to use for 'other' group
#' @return Vector of colors
#' @noRd
.set_colors <- function(df_in, data_col, plot_colors, plot_lvls = NULL,
                        other_label = "other", other_color = "grey60") {

  if (is.null(data_col)) return(plot_colors)

  # Return default numerical colors
  dat <- df_in[[data_col]]

  if (is.numeric(dat)) {
    plot_colors <- plot_colors %||% c("#132B43", "#56B1F7")

    return(plot_colors)
  }

  set_other_color <- !other_label %in% names(plot_colors)

  # Set default plot_lvls and plot_colors
  if (is.null(plot_lvls)) {
    plot_lvls <- levels(dat) %||% sort(unique(dat))
  }

  clrs <- scales::hue_pal()(length(plot_lvls))

  # Set plot_lvls as names
  lvls <- stats::na.omit(plot_lvls)
  len  <- min(length(lvls), length(clrs))

  clrs <- purrr::set_names(
    clrs[seq_len(len)], lvls[seq_len(len)]
  )

  # If plot_colors provided, replace default colors
  if (!is.null(plot_colors)) {
    if (is.null(names(plot_colors))) {
      plot_colors <- purrr::set_names(
        plot_colors, names(clrs)[seq_along(plot_colors)]
      )
    }

    clrs[names(plot_colors)] <- plot_colors
  }

  if (set_other_color) clrs[[other_label]] <- other_color

  clrs
}

#' Check cluster_col and group_col arguments
#' @noRd
.check_group_cols <- function(cluster_col, group_col, input = NULL,
                              uniq = TRUE) {

  if (!is.null(group_col) && is.null(cluster_col)) {
    cli::cli_abort(
      "`cluster_col` must be provided when `group_col` is specified"
    )
  }

  if (uniq && !is.null(group_col) && identical(group_col, cluster_col)) {
    cli::cli_abort(
      "`cluster_col` and `group_col` must specify different columns"
    )
  }

  if (!is.null(cluster_col) && !is.null(group_col) && !is.null(input)) {
    dat <- .get_meta(input)

    chk <- .check_matching_vals(dat[[cluster_col]], dat[[group_col]])

    if (!chk) {
      cli::cli_abort(
        "There must be a single group label for each cluster,
         i.e. each cluster can only belong to one group,
         check the values in `cluster_col` and `group_col`"
      )
    }
  }
}

.check_matching_vals <- function(x, y) {
  if (length(x) != length(y)) {
    cli::cli_abort("`x` and `y` must be the same length", .internal = TRUE)
  }

  res <- paste0(x, y)

  dplyr::n_distinct(x) == dplyr::n_distinct(res)
}

.get_matching_clmns <- function(df, clmn) {
  dat <- as.list(df)
  clmns <- names(dat)

  if (any(!clmn %in% clmns)) {
    cli::cli_abort("`clmn` is not present in `df`", .internal = TRUE)
  }

  clmns <- clmns[!clmns %in% clmn]
  clmns <- dat[clmns]

  clmn <- dat[clmn]
  clmn <- purrr::reduce(clmn, paste0)

  mtch <- purrr::map_lgl(clmns, ~ .check_matching_vals(clmn, .x))

  names(clmns[mtch])
}


#' Standardize aesthetics
#' e.g. color and colour
#' @noRd
.standardize_aes <- function(aes_list) {
  names(aes_list) <- sub("^color$", "colour", names(aes_list))

  aes_list
}

