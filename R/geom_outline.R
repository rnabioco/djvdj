#' Highlight points
#'
#' @param outline_size Outline size
#' @param outline_color Outline color
#' @param outline_position Outline position
#' @export
geom_outline <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity",
                         ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, geom = "point",
                         outline_size = 1, outline_color = "black", outline_position = "bottom") {

  geoms          <- list(GeomPoint, GeomLine, "geom", "line")
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

  plt_dat <- plot$data
  plt_map <- plot$mapping

  # Pull plot levels
  clr_col <- as_name(plt_map$colour)

  lvls <- pull(plt_dat, clr_col)
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
    lyr_dat <- filter(plt_dat, !!sym(clr_col) == .x)

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


tiny_so@meta.data %>%
  mutate(orig.ident = forcats::fct_relevel(orig.ident, c("avid_2", "avid_1"))) %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, color = orig.ident)) +
  geom_outline(
    outline_size = 11,
    size         = 10,
    geom         = "line"
  )
