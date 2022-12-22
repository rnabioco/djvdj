#' ggplot2 imports
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_histogram
#' @importFrom ggplot2 geom_density geom_tile geom_boxplot geom_violin geom_col
#' @importFrom ggplot2 scale_x_discrete scale_y_continuous scale_x_continuous
#' @importFrom ggplot2 position_dodge scale_color_manual scale_fill_manual
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
#' plot_features(vdj_so, feature = "seurat_clusters") +
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


#' Create ggplot heatmap
#'
#' @param df_in data.frame
#' @param x Variable to plot on the x-axis
#' @param y Variable to plot on the y-axis
#' @param .fill Variable to use for the fill color
#' @param clrs Vector of colors for plotting
#' @param na_color Color to use for missing values
#' @param plt_ttl Plot title
#' @param trans Method to use for transforming data
#' @param lgd_ttl Legend title
#' @param ang Angle of x-axis text
#' @param hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_gg_heatmap <- function(df_in, x = NULL, y, .fill, clrs = NULL,
                               na_color = NA, plt_ttl = ggplot2::waiver(),
                               trans = "identity", lgd_ttl = .fill, ang = 45,
                               hjst = 1, ...) {

  clrs <- clrs %||% "#619CFF"

  if (length(clrs) == 1) clrs <- c("white", clrs)

  plt_aes <- ggplot2::aes("sample", !!sym(y), fill = !!sym(.fill))

  if (!is.null(x)) plt_aes$x <- sym(x)

  res <- ggplot2::ggplot(df_in, plt_aes) +
    ggplot2::geom_tile(...) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = lgd_ttl)) +
    ggplot2::scale_fill_gradientn(
      colors   = clrs,
      na.value = na_color,
      trans    = trans
    ) +
    djvdj_theme() +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.line   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = ang, hjust = hjst)
    ) +
    labs(title = plt_ttl)

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
      warning(
        "Cannot remove diagonal since there ",
        "will only be one unique value remaining."
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
    stop(
      "Matrix must be symmetrical to remove upper triangle and/or diagonal."
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
    stop(
      "To create a circos plot there must be at least one overlap ",
      "between the clusters."
    )
  }

  # Set levels
  all_nms <- union(colnames(mat_in), rownames(mat_in))

  if (!is.null(lvls) && !all(all_nms %in% lvls)) {
    stop("plot_levels must include all rows and columns in the matrix")
  }

  if (!is.null(grps)) {
    if (!all(all_nms %in% names(grps))) {
      stop("Groups must be provided for all clusters.")
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
      stop(
        "Not enough colors provided (", n_clrs, "), provide exactly one ",
        "color, one color for each column in the matrix (", length(c_nms),
        "), or one color for each unique row and column in the matrix (",
        length(all_nms), ")."
      )
    }

  } else {
    if (!all(all_nms %in% names(clrs))) {
      stop("A color must be provided for each group being plotted")
    }

    c_clrs <- clrs[c_nms]
  }

  res <- list(clrs, c_clrs)

  res
}


#' Create ggplot bar graph
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param trans Method to use for transforming data
#' @param y_exp Specification to pass to [ggplot2::expansion()]
#' @param y_ttl Title for y-axis
#' @param ang Angle of x-axis text
#' @param hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_bars <- function(df_in, x, y, .fill = NULL, clrs = NULL,
                         trans = "identity", y_ttl = y, y_exp = c(0.05, 0.05),
                         ang = 45, hjst = 1, ...) {

  # Set aesthetics and geom_col arguments
  gg_aes <- aes(!!sym(x), !!sym(y))

  gg_args <- list(...)

  if (!is.null(.fill)) gg_aes$fill  <- sym(.fill)
  else                 gg_args$fill <- clrs

  if (is.null(gg_args$position)) {
    gg_args$position <- ggplot2::position_dodge(preserve = "single")
  }

  # Create bar graph
  y_args <- list(trans = trans, expand = ggplot2::expansion(y_exp))

  res <- ggplot2::ggplot(df_in, gg_aes) +
    lift(geom_col)(gg_args) +
    lift(ggplot2::scale_y_continuous)(y_args)

  # Set plot colors
  if (!is.null(.fill) && !is.null(clrs)) {
    res <- res +
      ggplot2::scale_fill_manual(values = clrs)
  }

  # Set theme
  res <- res +
    ggplot2::labs(y = y_ttl) +
    djvdj_theme() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(angle = ang, hjust = hjst)
    )

  res
}


#' Create ggplot boxplot
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param grp Variable to use for facet_wrap
#' @param .color Variable to use for line color
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param method Method to use for plotting, either 'boxplot' or 'violin'
#' @param trans Transformation to use for plotting data, e.g. 'log10', refer
#' to [ggplot2::continuous_scale()] for more options.
#' @param y_exp Specification to pass to [ggplot2::expansion()]
#' @param nrow Number of rows for facets
#' @param scales scales specification for facet_wrap
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_boxes <- function(df_in, x = NULL, y, grp = NULL, .color = NULL,
                          .fill = NULL, clrs = NULL, method = "boxplot",
                          trans = "identity", y_exp = c(0.05, 0.05),
                          nrow = NULL, scales = "fixed", ...) {

  # Check input
  typs <- c("boxplot", "violin")

  if (!method %in% typs) {
    stop("method must be one of ", paste0(typs, collapse = ", "))
  }

  # Set aesthetics
  plt_aes <- ggplot2::aes(y, !!sym(y))

  if (!is.null(x))      plt_aes$x      <- sym(x)
  if (!is.null(.fill))  plt_aes$fill   <- sym(.fill)
  if (!is.null(.color)) plt_aes$colour <- sym(.color)

  # Adjust theme and transform axis
  res <- ggplot2::ggplot(df_in, plt_aes) +
    djvdj_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x    = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous(
      trans = trans, expand = ggplot2::expansion(y_exp)
    )

  if (!is.null(grp)) {
    res <- res +
      ggplot2::facet_wrap(
        stats::as.formula(paste("~", grp)),
        nrow = nrow,
        scales = scales
      )
  }

  if (!is.null(clrs)) {
    res <- res +
      ggplot2::scale_color_manual(values = clrs) +
      ggplot2::scale_fill_manual(values = clrs)
  }

  # Return violins
  if (identical(method, "violin")) {
    res <- res +
      ggplot2::geom_violin(...) +
      ggplot2::stat_summary(
        geom = "point", fun = stats::median, color = "black"
      )

    return(res)
  }

  # Return boxes
  res <- res +
    ggplot2::geom_boxplot(...)

  res
}


#' Create ggplot histogram
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param grp Variable to use for facet_wrap
#' @param .color Variable to use for line color
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param method Method to use for plotting, either 'histogram' or 'density'
#' @param units Units to use for y-axis when plotting histogram. Use
#' 'frequency' to show the raw number of values or 'percent' to show
#' the percentage of total values.
#' @param trans Transformation to use for plotting data, e.g. 'log10', refer
#' to [ggplot2::continuous_scale()] for more options.
#' @param y_exp Specification to pass to [ggplot2::expansion()]
#' @param nrow Number of rows for facets
#' @param scales scales specification for facet_wrap
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_hist <- function(df_in, x, grp, .color = NULL, .fill = NULL,
                         clrs = NULL, method = "histogram",
                         units = "frequency", y_ttl = units, trans = "identity",
                         y_exp = c(0.05, 0.1), nrow = NULL, scales = "fixed",
                         ...) {

  # Check inputs
  typs <- c("histogram", "density")

  if (!method %in% typs) {
    stop("method must be one of ", paste0(typs, collapse = ", "))
  }

  axs <- c("frequency", "percent")

  if (!units %in% axs) {
    stop("units must be one of ", paste0(axs, collapse = ", "))
  }

  # Set aesthetics
  plt_aes <- ggplot2::aes()

  # Only plot percent for histogram
  if (identical(units, "percent") && identical(method, "histogram")) {
    plt_aes <- ggplot2::aes(y = ggplot2::after_stat((count / max(count)) * 100))
  }

  plt_aes$x <- sym(x)

  if (!is.null(.fill))  plt_aes$fill   <- sym(.fill)
  if (!is.null(.color)) plt_aes$colour <- sym(.color)

  # Adjust theme and transform axis
  res <- ggplot2::ggplot(df_in, plt_aes) +
    djvdj_theme() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::scale_x_continuous(trans = trans) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(y_exp))

  if (!is.null(grp)) {
    res <- res +
      ggplot2::facet_wrap(
        stats::as.formula(paste("~", grp)),
        nrow = nrow,
        scales = scales
      )
  }

  if (!is.null(clrs)) {
    res <- res +
      ggplot2::scale_color_manual(values = clrs) +
      ggplot2::scale_fill_manual(values = clrs)
  }

  # Return density plot
  if (identical(method, "density")) {
    res <- res +
      ggplot2::geom_density(...)

    return(res)
  }

  # Return histogram
  res <- res +
    ggplot2::geom_histogram(...) +
    ggplot2::labs(y = y_ttl)

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


#' Add n label to plot
#'
#' @param gg_in ggplot2 object
#' @param n n to use for label
#' @param lab_args named list with aesthetic parameters to used for modifying
#' n label
#' @noRd
.add_n_label <- function(gg_in, lab_args, n = NULL) {
  lab_args$geom  <- "text"
  lab_args$x     <- lab_args$x %||% Inf
  lab_args$y     <- lab_args$y %||% Inf
  lab_args$hjust <- lab_args$hjust %||% 1.2
  lab_args$vjust <- lab_args$vjust %||% 1.5

  lab_args$label <- lab_args$label %||%
    paste0("n = ", scales::label_comma()(n))

  res <- gg_in +
    lift(ggplot2::annotate)(lab_args)

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
    stop(ft, " is not numeric")
  }

  if (is.null(mn) && is.null(mx)) {
    stop("mn or mx must be provided")
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
    warning(
      "Plot levels were not modified, ",
      "levels are only modified for characters and factors"
    )
    return(df_in)
  }

  if (any(is.na(dat)) && !any(is.na(lvls))) lvls <- c(NA, lvls)

  u_dat     <- unique(dat)
  missing   <- u_dat[!u_dat %in% lvls]
  n_missing <- length(missing)

  if (n_missing > 0) {
    sfx <- ""

    if (n_missing > 10) {
      sfx <- "..."
      missing <- missing[seq_len(10)]
    }

    missing <- paste0(missing, collapse = ", ")
    missing <- paste0(missing, sfx)

    stop("Not all labels in ", clmn, " are included in plot_lvls: ", missing)
  }

  df_in[clmn] <- factor(df_in[[clmn]], unique(lvls), exclude = NULL)

  df_in
}


#' Check cluster_col and group_col arguments
#' @noRd
.chk_group_cols <- function(cluster_col, group_col, input = NULL) {
  if (!is.null(group_col) && is.null(cluster_col)) {
    stop("cluster_col must be provided when group_col is specified.")
  }

  if (!is.null(group_col) && identical(group_col, cluster_col)) {
    stop("group_col and cluster_col must specify different columns.")
  }

  if (!is.null(cluster_col) && !is.null(group_col) && !is.null(input)) {
    dat <- .get_meta(input)

    chk <- .chk_matching_vals(dat[[cluster_col]], dat[[group_col]])

    if (!chk) {
      stop(
        "There must be a single group label for each cluster, ",
        "i.e. each cluster can only belong to one group. ",
        "Check the values in group_col and cluster_col."
      )
    }
  }
}

.chk_matching_vals <- function(x, y) {
  if (length(x) != length(y)) stop("x and y must be the same length.")

  res <- paste0(x, y)

  dplyr::n_distinct(x) == dplyr::n_distinct(res)
}

.get_matching_clmns <- function(df, clmn) {
  dat <- as.list(df)

  clmns <- names(dat)
  clmns <- clmns[!clmns %in% clmn]
  clmns <- dat[clmns]

  clmn <- dat[clmn]
  clmn <- purrr::reduce(clmn, paste0)

  mtch <- purrr::map_lgl(clmns, ~ .chk_matching_vals(clmn, .x))

  names(clmns[mtch])
}

