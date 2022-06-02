
#' Create circos plot showing overlap of clonotypes between clusters
#'
#' @param input input data
#' @param cluster_col cluster column
#' @param clonotype_col clonotype column
#' @export
plot_circos <- function(input, cluster_col, clonotype_col = "clonotype_id",
                        plot_colors = NULL, plot_lvls = names(plot_colors),
                        width = 0.2) {

  # Format input data
  meta <- .get_meta(input)
  vdj  <- dplyr::filter(meta, !is.na(!!sym(clonotype_col)))

  vdj <- dplyr::select(
    vdj,
    all_of(c(CELL_COL, clonotype_col, cluster_col))
  )

  # Remove very small groups
  frac_lim   <- 0.03
  clsts_n    <- table(vdj[[cluster_col]])
  clst_frac  <- clsts_n / sum(clsts_n)
  tiny_clsts <- names(clst_frac[clst_frac < frac_lim])

  if (!is_empty(tiny_clsts)) {
    warning(
      "Some clusters are too small to plot and were removed: ",
      paste0(tiny_clsts, collapse = ", ")
    )
  }

  vdj <- dplyr::filter(vdj, !(!!sym(cluster_col) %in% tiny_clsts))

  # Set plot levels
  # run .set_lvls twice so levels are still set whe plot_lvls is NULL
  vdj <- .set_lvls(vdj, cluster_col, plot_lvls)

  clsts <- unique(vdj[[cluster_col]])
  clsts <- as.character(sort(clsts))

  vdj <- .set_lvls(vdj, cluster_col, clsts)

  vdj   <- dplyr::arrange(
    vdj,
    !!sym(cluster_col),
    desc(dplyr::n_distinct(!!sym(clonotype_col))),
    !!sym(clonotype_col)
  )

  # Sort and number rows for each cluster
  # get coordinate for labels
  vdj <- dplyr::mutate(vdj, id = row_number())

  vdj <- dplyr::group_by(vdj, !!sym(cluster_col))
  vdj <- dplyr::mutate(vdj, coord = median(id))
  vdj <- dplyr::ungroup(vdj)

  lab_coords <- dplyr::distinct(vdj, !!sym(cluster_col), coord)
  lab_coords <- purrr::set_names(lab_coords$coord, lab_coords[[cluster_col]])

  # Set plot_colors
  if (!is.null(plot_colors)) {
    if (!is.null(names(plot_colors))) {
      plot_colors <- plot_colors[clsts]

    } else {
      plot_colors <- plot_colors[seq_along(clsts)]
    }

  } else {
    plot_colors <- NA
  }

  # Get ID ranges for each clonotype
  ids <- dplyr::select(vdj, -all_of(CELL_COL))
  ids <- dplyr::group_by(ids, !!!syms(c(clonotype_col, cluster_col)))

  ids <- dplyr::summarize(
    ids,
    id = list(unique(range(.data$id))),
    n  = n(),
    .groups = "drop"
  )

  # Cluster combinations to plot
  combs <- utils::combn(clsts, 2, simplify = FALSE)

  # Create circos plot
  circlize::circos.clear()
  circlize::circos.par("track.height" = width)
  circlize::circos.initialize(vdj[[cluster_col]], vdj$id)

  circlize::circos.track(
    ylim   = c(0, 1),
    bg.col = plot_colors
  )

  # JUST USE LEGEND INSTEAD
  # circlize::circos.text(
  #   x      = unname(lab_coords),
  #   y      = rep(2, length(clsts)),
  #   labels = names(lab_coords),
  #   facing = "downward"
  # )

  purrr::walk(combs, ~ {
    d <- dplyr::group_by(ids, !!sym(clonotype_col))
    d <- dplyr::filter(d, all(.x %in% !!sym(cluster_col)))
    d <- dplyr::filter(d, !!sym(cluster_col) %in% .x)

    d <- split(d, d[[clonotype_col]])  # THERE MUST BE A BETTER WAY

    purrr::walk(d, ~ {
      lnk_clr <- dplyr::slice_max(.x, n, n = 1, with_ties = FALSE)
      lnk_clr <- plot_colors[lnk_clr[[cluster_col]]]

      circlize::circos.link(
        .x[[cluster_col]][1], .x$id[[1]],
        .x[[cluster_col]][2], .x$id[[2]],
        col = unname(lnk_clr)

      )
    })
  })
}






# sect <- table
#
# circos.clear()
#
# circos.clear()
# circos.par("track.height" = 0.1)
# circos.initialize(df$sectors, x = df$x)
#
# # add track
# circos.initialize(
#   sectors = dat$orig.ident,
#   x = dat$clonotype_id
# )
#
# circos.track(ylim = c(0, 1))
#
# # Example from vignette
# set.seed(999)
#
# n = 1000
#
# df = data.frame(
#   sectors = sample(letters[1:8], n, replace = TRUE),
#   x = rnorm(n),
#   y = runif(n)
# )
#
# # initialize
# circos.clear()
# circos.par("track.height" = 0.1)
# circos.initialize(df$sectors, x = df$x)
#
# # add track
# circos.track(
#   df$sectors,
#   y = df$y
#   # panel.fun = function(x, y) {
#   #   circos.text(
#   #     CELL_META$xcenter,
#   #     CELL_META$cell.ylim[2] + mm_y(5),
#   #     CELL_META$sector.index
#   #   )
#   #   circos.axis(labels.cex = 0.6)
#   # }
# )
#
# # col = rep(c("#FF0000", "#00FF00"), 4)
# # circos.trackPoints(df$sectors, df$x, df$y, col = col, pch = 16, cex = 0.5)
# # circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1)
#
# # add links
# # circos.link("a", 0, "b", 0, h = 0.4)
# # circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red", border = "blue", h = 0.2)
# circos.link(
#   "e", c(-2, 1), "g", c(-1,1),
#   # "e", 0, "g", c(-1,1),
#   col = "green",
#   lwd = 2
#   # lty = 2
# )

