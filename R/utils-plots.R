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
#' @param ttl_size Size of axis titles
#' @param txt_size Size of axis text
#' @param ln_size Size of axis lines
#' @param txt_col Color of axis text
#' @param ln_col Color of axis lines
#' @return ggplot theme
#'
#' @examples
#'
#' plot_scatter(vdj_so, data_col = "seurat_clusters") +
#'   djvdj_theme()
#'
#' @export
djvdj_theme <- function(ttl_size = 12, txt_size = 8, ln_size = 0.5,
                        txt_col = "black", ln_col = "grey85") {

  res <- ggplot2::theme(
    strip.background  = ggplot2::element_blank(),
    strip.text        = ggplot2::element_text(size = ttl_size),

    panel.border      = ggplot2::element_rect(
      fill = NA, color = ln_col, linewidth = ln_size
    ),

    panel.background  = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.title      = ggplot2::element_text(size = ttl_size),
    legend.key        = ggplot2::element_blank(),
    legend.text       = ggplot2::element_text(size = txt_size, color = txt_col),
    axis.line         = ggplot2::element_blank(),
    axis.ticks        = ggplot2::element_line(linewidth = ln_size,  color = ln_col),
    axis.text         = ggplot2::element_text(size = txt_size, color = txt_col),
    axis.title        = ggplot2::element_text(size = ttl_size, color = txt_col)
  )

  res
}


#' Create ggplot
#'
#' @param df_in data.frame
#' @param fn ggplot2 function to use for generating plot,
#' e.g. `ggplot2::geom_point`
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param grp Varible to use for `ggplot2::facet_wrap()`
#' @param .color Variable to use for color
#' @param .fill Variable to use for fill
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

  # Check inputs
  if (!is.null(.fill) && !is.null(.color) && !identical(.fill, .color)) {
    cli::cli_abort(
      "If both `.color` and `.fill` are provided,
       they must refer to the same variable", .internal = TRUE
    )
  }

  # Set aesthetics and geom arguments
  scale_fill <- !is.null(.fill)
  scale_clr  <- !is.null(.color)
  num_clr    <- scale_clr && is.numeric(df_in[[.color]])
  num_fill   <- scale_fill && is.numeric(df_in[[.fill]])

  gg_aes  <- ggplot2::aes()
  gg_args <- .standardize_aes(list(...))

  if (is.null(x)) gg_aes$x <- y %||% "x"
  else            gg_aes$x <- sym(x)

  if (!is.null(y)) gg_aes$y      <- sym(y)
  if (scale_fill)  gg_aes$fill   <- sym(.fill)
  if (scale_clr)   gg_aes$colour <- sym(.color)

  if (!scale_fill && !scale_clr) gg_args$colour <- gg_args$colour %||% clrs

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

  res <- .add_n_label(
    res, label_data,
    n_label   = n_label,
    crnr_col  = grp,
    axis_col  = .chk_num(df_in, x),
    lgnd_col  = .chk_num(df_in, .fill %||% .color),
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
#' @param trans Method to use for transforming y-axis
#' @param y_ttl Title for y-axis
#' @param x_ang Angle of x-axis text
#' @param x_hjst Horizontal justification for x-axis text
#' @param n_label Vector indicating where n labels should be added
#' @param ... Additional arguments to pass to `.create_plot()`
#' @return ggplot object
#' @noRd
.create_bars <- function(df_in, x = NULL, y, trans = "identity", y_ttl = y,
                         x_ang = NULL, x_hjst = NULL, n_label = NULL, ...) {

  # Check input values
  if (!is.null(x) && is.numeric(df_in[[x]])) {
    cli::cli_abort("`x` cannot refer to a numeric column")
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

  if (is.null(gg_args$position)) {
    gg_args$position <- ggplot2::position_dodge2(
      preserve = "single", padding = 0
    )
  }

  res <- lift(.create_plot)(gg_args)

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

  if (identical(method, "violin")) pos <- ggplot2::position_dodge()
  else pos <- ggplot2::position_dodge2(preserve = "single")

  gg_args$position <- gg_args$position %||% pos

  # Create boxplots
  # allow user to set point size and color
  if (show_points && identical(method, "boxplot")) gg_args$outlier.color <- NA

  res <- lift(.create_plot)(gg_args)

  # Add additional points
  if (show_points) {
    pt_args       <- list()
    pt_args$alpha <- 1
    pt_args$size  <- point.size %||% 1
    pt_args$color <- point.color %||% gg_args$color

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

  plt_fn <- plt_fns[[method]]

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
  res <- .create_plot(
    df_in, plt_fn,
    mapping = gg_aes,
    x       = x,
    .color  = .color,
    .fill   = .fill,
    trans_x = trans,
    n_label = n_label,
    ...
  )

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
#' @importFrom circlize chordDiagram circos.clear circos.axis circos.track mm_h
#' @importFrom graphics title strwidth
#' @noRd
.create_circos <- function(mat_in, clrs = NULL, na_color = "grey90",
                           lvls = NULL, grps = NULL, rotate_labels = FALSE,
                           plt_ttl = NULL, ...) {

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
  } else {
    if (track_ht) {
      circos_fun <- function(...) {
        circlize::chordDiagram(
          annotationTrackHeight = circlize::mm_h(c(3, 4)),
          ...
        )
      }
    } else {
      circos_fun <- function(...) circlize::chordDiagram(...)
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

  lab_args <- .standardize_aes(lab_args)

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
  just <- 0.5
  char_h_w <- 1.5

  dat <- .format_n_label(df_in)

  dat <- dplyr::mutate(
    dat,
    label = lab_args$label %||% .data$label,
    hjust = nchar(.data$label),
    hjust = 1 + (1 / .data$hjust * (just * char_h_w))
  )

  if (!is.null(lab_args$size)) lab_args$size <- lab_args$size / ggplot2::.pt

  lab_args$mapping     <- ggplot2::aes(label = .data$label, hjust = .data$hjust)
  lab_args$data        <- dat
  lab_args$inherit.aes <- FALSE
  lab_args$x           <- lab_args$x %||% Inf
  lab_args$y           <- lab_args$y %||% Inf
  lab_args$vjust       <- lab_args$vjust %||% 1 + just

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
      ggplot2::scale_x_discrete(labels = dat_labs) +
      ggplot2::theme(axis.text.x = lift(ggplot2::element_text)(lab_args))

  } else if (identical(axis, "y")) {
    res <- gg_in +
      ggplot2::scale_y_discrete(labels = dat_labs) +
      ggplot2::theme(axis.text.y = lift(ggplot2::element_text)(lab_args))

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
    met = list(boxplot.stats(
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
  # remove levels that do not appear in data
  # If user provides NA in plot_levels, keep in order, otherwise always order
  # NA last
  not_in    <- rnk[!rnk %in% plot_lvls]
  plot_lvls <- plot_lvls[plot_lvls %in% c(rnk, NA)]
  lvls      <- c(plot_lvls, not_in)

  if (!NA %in% lvls) lvls <- c(lvls, NA)

  if (rev && !purrr::is_empty(not_in)) lvls <- rev(lvls)

  res <- .set_lvls(res, data_col, lvls)

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
#'
#' e.g. color and colour
#' @noRd
.standardize_aes <- function(aes_list) {
  names(aes_list) <- sub("color", "colour", names(aes_list), fixed = TRUE)

  aes_list
}

