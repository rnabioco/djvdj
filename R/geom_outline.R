#' Highlight points
#'
#' @param mapping Set of aesthetic mappings created by [aes()] or
#'    [aes_()]. If specified and `inherit.aes = TRUE` (the
#'    default), it is combined with the default mapping at the top level of the
#'    plot. You must supply `mapping` if there is no plot mapping.
#' @param data The data to be displayed in this layer. There are three
#'    options:
#'
#'    If `NULL`, the default, the data is inherited from the plot
#'    data as specified in the call to [ggplot()].
#'
#'    A `data.frame`, or other object, will override the plot
#'    data. All objects will be fortified to produce a data frame. See
#'    [fortify()] for which variables will be created.
#'
#'    A `function` will be called with a single argument,
#'    the plot data. The return value must be a `data.frame`, and
#'    will be used as the layer data. A `function` can be created
#'    from a `formula` (e.g. `~ head(.x, 10)`).
#' @param geom The geometric object to use display the data
#' @param stat The statistical transformation to use on the data for this
#'    layer, as a string.
#' @param position Position adjustment, either as a string, or the result of
#'    a call to a position adjustment function.
#' @param ... Other arguments passed on to [layer()]. These are often
#'     aesthetics, used to set an aesthetic to a fixed value, like
#'     `colour = "red"` or `size = 3`. They may also be parameters to the
#'     paired geom/stat.
#' @param na.rm If `FALSE`, the default, missing values are removed with a
#'     warning. If `TRUE`, missing values are silently removed.
#' @param show.legend logical. Should this layer be included in the legends?
#'    `NA`, the default, includes if any aesthetics are mapped.
#'    `FALSE` never includes, and `TRUE` always includes.
#'    It can also be a named logical vector to finely select the aesthetics to
#'    display.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics,
#'    rather than combining with them. This is most useful for helper functions
#'    that define both data and aesthetics and shouldn't inherit behaviour from
#'    the default plot specification, e.g. [borders()].
#' @param outline_size Outline size
#' @param outline_color Outline color
#' @param outline_position Outline position
#'    * "all", a separate outline is drawn for each group plotted
#'    * "bottom", a single outline is drawn around all groups
#'    * A vector of group names can be passed to only draw outlines around
#'      certain groups
#' @export
geom_outline <- function(mapping = NULL, data = NULL, geom = "point", stat = "identity", position = "identity",
                         ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, outline_size = 1,
                         outline_color = "black", outline_position = "all") {

  geoms          <- list(GeomPoint, GeomLine, "point", "line")
  geom_supported <- any(purrr::map_lgl(geoms, identical, geom))

  if (!geom_supported) {
    stop("geom_outline currently only supports point and line geoms")
  }

  structure(
    list(
      mapping          = mapping,
      data             = data,
      stat             = stat,
      position         = position,
      na.rm            = na.rm,
      show.legend      = show.legend,
      inherit.aes      = inherit.aes,
      geom             = geom,
      outline_size     = outline_size,
      outline_color    = outline_color,
      outline_position = outline_position,
      other_params     = list(...)
    ),
    class = "geom_outline"
  )
}

#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.geom_outline <- function(object, plot, object_name) {

  # Identify data to use for grouping layers
  plt_map <- plot$mapping
  plt_dat <- plot$data

  grps <- dplyr::case_when(
    "group" %in% names(plt_map)  ~ "group",
    "colour" %in% names(plt_map) ~ "colour"
  )

  if (!is.na(grps)) {
    grps    <- as_name(plt_map[[grps]])
    grp_dat <- dplyr::pull(plt_dat, grps)
  }

  # Modify plot groups based on outline_position
  out_pos <- object$outline_position

  if (out_pos == "bottom" || is.na(grps)) {
    plt_dat <- dplyr::mutate(plt_dat, .GROUPS = "BOTTOM")

  } else if (out_pos == "all") {
    plt_dat <- dplyr::mutate(plot$data, .GROUPS = !!sym(grps))

  } else {
    if (!all(out_pos %in% grp_dat)) {
      stop("outline_position must be 'all', 'bottom', or a label(s) present in the column passed to group or color within aes()")
    }

    plt_dat <- dplyr::mutate(
      plt_dat,
      .GROUPS = ifelse(!!sym(grps) %in% out_pos, as.character(!!sym(grps)), "OTHER"),
      .GROUPS = factor(.GROUPS, c("OTHER", out_pos))
    )
  }

  # Pull plot levels
  lvls <- pull(plt_dat, .GROUPS)
  lvls <- unique(lvls)

  if (is.factor(lvls)) {
    lvls <- levels(lvls)
  }

  # Plot parameters
  na <- list(na.rm = object$na.rm)

  out_prms <- append(
    list(
      color = object$outline_color,
      size  = object$outline_size
    ),
    na
  )

  other_prms <- append(na, object$other_params)

  # Create layers
  lyrs <- purrr::map(lvls, ~ {
    lyr_dat <- filter(plt_dat, .GROUPS == .x)

    list(
      layer(
        data        = lyr_dat,
        geom        = object$geom,
        mapping     = object$mapping,
        stat        = object$stat,
        position    = object$position,
        show.legend = object$show.legend,
        inherit.aes = object$inherit.aes,
        params      = out_prms
      ),
      layer(
        data        = lyr_dat,
        geom        = object$geom,
        mapping     = object$mapping,
        stat        = object$stat,
        position    = object$position,
        show.legend = object$show.legend,
        inherit.aes = object$inherit.aes,
        params      = other_prms
      )
    )
  })

  lyrs <- purrr::flatten(lyrs)

  # Add layers to plot
  res <- plot + lyrs

  res
}

