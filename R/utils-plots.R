#' ggplot2 imports
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_histogram geom_density geom_tile geom_boxplot geom_violin geom_col
#' @importFrom ggplot2 scale_x_discrete scale_y_continuous scale_x_log10 scale_y_log10 position_dodge
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_color_gradientn scale_fill_gradientn
#' @importFrom ggplot2 stat_summary facet_wrap guides guide_legend labs
#' @importFrom ggplot2 theme element_blank element_text element_line
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
#' @export
djvdj_theme <- function(ttl_size = 12, txt_size = 8, ln_size = 0.5, txt_col = "black",
                        ln_col = "grey85") {

  res <- ggplot2::theme(
    strip.background  = ggplot2::element_blank(),
    strip.text        = ggplot2::element_text(size = ttl_size),

    panel.border      = ggplot2::element_rect(fill = NA, color = ln_col, size = ln_size),

    panel.background  = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.title      = ggplot2::element_text(size = ttl_size),
    legend.key        = ggplot2::element_blank(),
    legend.text       = ggplot2::element_text(size = txt_size, color = txt_col),

    axis.line         = ggplot2::element_blank(),
    # axis.line         = ggplot2::element_line(size = ln_size,  color = ln_col),

    axis.ticks        = ggplot2::element_line(size = ln_size,  color = ln_col),
    axis.text         = ggplot2::element_text(size = txt_size, color = txt_col),
    axis.title        = ggplot2::element_text(size = ttl_size, color = txt_col)
  )

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
#' @param ttl Legend title
#' @param ang Angle of x-axis text
#' @param hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_heatmap <- function(df_in, x, y, .fill, clrs = NULL, na_color = "white",
                            ttl = .fill, ang = 45, hjst = 1, ...) {

  if (is.null(clrs)) {
    clrs <- "#6A51A3"
  }

  if (length(clrs) == 1) {
    clrs <- c("grey90", clrs)
  }

  if (!is.null(x)) {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes(!!sym(x), !!sym(y), fill = !!sym(.fill))
    )

  } else {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes("sample", !!sym(y), fill = !!sym(.fill))
    )
  }

  res <- res +
    ggplot2::geom_tile(...) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = ttl)) +
    ggplot2::scale_fill_gradientn(colors = clrs, na.value = na_color) +
    djvdj_theme() +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.line   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = ang, hjust = hjst)
    )

  res
}


#' Create ggplot bar graph
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param y_ttl Title for y-axis
#' @param ang Angle of x-axis text
#' @param hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_bars <- function(df_in, x, y, .fill, clrs = NULL, y_ttl = y, ang = 45,
                         hjst = 1, ...) {

  # Bar position to use for plot
  bar_pos <- ggplot2::position_dodge(preserve = "single")

  # Reverse bar order
  lvls  <- rev(levels(pull(df_in, x)))
  df_in <- .set_lvls(df_in, x, lvls)

  # Color fill by variable
  if (!is.null(.fill)) {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes(!!sym(x), !!sym(y), fill = !!sym(.fill))
    ) +
    ggplot2::geom_col(..., position = bar_pos)

    if (!is.null(clrs)) {
      res <- res +
        ggplot2::scale_fill_manual(values = clrs)
    }

  # Use single fill color
  } else {
    res <- ggplot2::ggplot(
      df_in,
      ggplot2::aes(!!sym(x), !!sym(y))
    )

    if (!is.null(clrs)) {
      res <- res +
        ggplot2::geom_col(..., fill = clrs, position = bar_pos)

    } else {
      res <- res +
        ggplot2::geom_col(..., position = bar_pos)
    }
  }

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
#' @param .color Variable to use for line color
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param type Type of plot to create, either 'boxplot' or 'violin'
#' @param log_trans If TRUE, log10 transform the y-axis
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_boxes <- function(df_in, x = NULL, y, .color = NULL, .fill = NULL, clrs = NULL,
                          type = "boxplot", log_trans = FALSE, ...) {

  # Check input
  typs <- c("boxplot", "violin")

  if (!type %in% typs) {
    stop("type must be one of ", paste0(typs, collapse = ", "))
  }

  # Set aesthetics
  plt_aes <- ggplot2::aes(y, !!sym(y))

  if (!is.null(x)) {
    plt_aes$x <- sym(x)
  }

  if (!is.null(.fill)) {
    plt_aes$fill <- sym(.fill)
  }

  if (!is.null(.color)) {
    plt_aes$colour <- sym(.color)
  }

  # Adjust theme
  res <- ggplot2::ggplot(df_in, plt_aes) +
    djvdj_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x    = ggplot2::element_blank()
    )

  if (!is.null(clrs)) {
    res <- res +
      ggplot2::scale_color_manual(values = clrs) +
      ggplot2::scale_fill_manual(values = clrs)
  }

  # Log10 tranform y-axis
  if (log_trans) {
    res <- res +
      ggplot2::scale_y_log10()
  }

  # Return violins
  if (identical(type, "violin")) {
    res <- res +
      ggplot2::geom_violin(...) +
      ggplot2::stat_summary(geom = "point", fun = stats::median, color = "black")

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
#' @param .color Variable to use for line color
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param type Type of plot to create, either 'histogram' or 'density'
#' @param yaxis Units to use for y-axis when plotting histogram. Use
#' 'frequency' to show the raw number of values or 'percent' to show
#' the percentage of total values.
#' @param log_trans If TRUE, log10 transform the x-axis
#' @param ... Additional arguments to pass to ggplot2, e.g. color, fill, size,
#' linetype, etc.
#' @return ggplot object
#' @noRd
.create_hist <- function(df_in, x, .color = NULL, .fill = NULL, clrs = NULL, type = "histogram",
                         yaxis = "frequency", log_trans = FALSE, ...) {

  # Check inputs
  typs <- c("histogram", "density")

  if (!type %in% typs) {
    stop("type must be one of ", paste0(typs, collapse = ", "))
  }

  axs <- c("frequency", "percent")

  if (!yaxis %in% axs) {
    stop("yaxis must be one of ", paste0(axs, collapse = ", "))
  }

  # Set aesthetics
  plt_aes <- ggplot2::aes()

  # Only plot percent for histogram
  if (identical(yaxis, "percent") && identical(type, "histogram")) {
    plt_aes <- ggplot2::aes(y = .data$..count.. / sum(.data$..count..) * 100)
  }

  plt_aes$x <- sym(x)

  if (!is.null(.fill)) {
    plt_aes$fill <- sym(.fill)
  }

  if (!is.null(.color)) {
    plt_aes$colour <- sym(.color)
  }

  # Adjust theme
  res <- ggplot2::ggplot(df_in, plt_aes) +
    djvdj_theme() +
    ggplot2::theme(legend.position = "top")

  if (!is.null(clrs)) {
    res <- res +
      ggplot2::scale_color_manual(values = clrs) +
      ggplot2::scale_fill_manual(values = clrs)
  }

  # Log10 transform axis
  if (log_trans) {
    res <- res +
      ggplot2::scale_x_log10()
  }

  # Return density plot
  if (identical(type, "density")) {
    res <- res +
      ggplot2::geom_density(...)

    return(res)
  }

  # Return histogram
  res <- res +
    ggplot2::geom_histogram(...) +
    ggplot2::labs(y = yaxis)

  res
}


#' Add list of aes params to ggplot object
#'
#' @param gg_in ggplot object
#' @param add_aes List of aes params to add or override
#' @param lay Layer number to modify
#' @return ggplot object
#' @noRd
.add_aes <- function(gg_in, add_aes, lay) {

  # Need to use colour instead of color
  aes_names      <- names(add_aes)
  names(add_aes) <- replace(aes_names, aes_names == "color", "colour")

  # Add aes params
  curr_aes <- gg_in$layers[[lay]]$aes_params

  curr_aes[names(add_aes)]       <- add_aes
  gg_in$layers[[lay]]$aes_params <- curr_aes

  gg_in
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

  res <- dplyr::select(res, -.data$pct)

  res
}


#' Set column levels
#'
#' @param df_in data.frame
#' @param clmn Column to modify
#' @param lvls Levels
#' @return data.frame
#' @noRd
.set_lvls <- function(df_in, clmn, lvls) {

  if (!is.null(lvls) && !is.null(clmn)) {
    dat <- pull(df_in, clmn)

    if (!is.character(dat) && !is.factor(dat)) {
      warning("Plot levels were not modified, levels are only modified for characters and factors")

      return(df_in)
    }

    if (!all(pull(df_in, clmn) %in% lvls)) {
      stop("Not all labels in ", clmn, " are included in plot_levels")
    }

    df_in <- dplyr::mutate(
      df_in,
      !!sym(clmn) := factor(!!sym(clmn), levels = unique(lvls))
    )
  }

  df_in
}

