#' Create two dimensional scatter plot
#' 
#' @param obj_in Seurat object or data.frame containing data for plotting
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param feature Variable to use for coloring points
#' @param data_slot Slot in the Seurat object to pull data from
#' @param split_id Variable to use for splitting plot
#' @return ggplot object
plot_features <- function(obj_in, x = "UMAP_1", y = "UMAP_2", feature, data_slot = "data", 
                          split_id = NULL,
                          
                          pt_size = 0.25, 
                          pt_outline = NULL, 
                          plot_cols = NULL,
                          feat_levels = NULL, 
                          split_levels = NULL, 
                          
                          min_pct = NULL,
                          max_pct = NULL,
                          
                          calc_cor = F, 
                          lm_line = F, 
                          lab_size = 3.7, 
                          lab_pos = c(0.8, 0.9),
                          
                          na_value = "#f0f0f0", 
                          ...) {
  
  # Format imput data
  counts <- obj_in
  
  if ("Seurat" %in% class(obj_in)) {
    vars <- c(x, y, feature)
    
    if (!is.null(split_id)) {
      vars <- c(vars, split_id)
    }
    
    counts <- obj_in %>%
      FetchData(vars = unique(vars), slot = data_slot) %>%
      as_tibble(rownames = "cell_ids")
  }
  
  # Rename features
  if (!is.null(names(feature))) {
    counts <- counts %>%
      rename(!!!syms(feature))
    
    feature <- names(feature)
  }
  
  if (!is.null(names(x))) {
    counts <- counts %>%
      rename(!!!syms(x))
    
    x <- names(x)
  }
  
  if (!is.null(names(y))) {
    counts <- counts %>%
      rename(!!!syms(y))
    
    y <- names(y)
  }
  
  # Set min and max values for feature
  if (!is.null(min_pct) || !is.null(max_pct)) {
    counts <- counts %>%
      mutate(
        pct_rank = percent_rank(!!sym(feature)),
        max_val  = ifelse(pct_rank > max_pct, !!sym(feature), NA),
        max_val  = min(max_val, na.rm = T),
        min_val  = ifelse(pct_rank < min_pct, !!sym(feature), NA),
        min_val  = max(min_val, na.rm = T),
        !!sym(feature) := if_else(!!sym(feature) > max_val, max_val, !!sym(feature)),
        !!sym(feature) := if_else(!!sym(feature) < min_val, min_val, !!sym(feature))
      )
  }
  
  # Set feature order
  if (!is.null(feat_levels)) {
    counts <- counts %>%
      mutate(!!sym(feature) := fct_relevel(!!sym(feature), feat_levels))
  }
  
  # Set facet order
  if (!is.null(split_id) && length(split_id) == 1) {
    counts <- counts %>%
      mutate(split_id = !!sym(split_id))
    
    if (!is.null(split_levels)) {
      counts <- counts %>%
        mutate(split_id = fct_relevel(split_id, split_levels))
    }
  }
  
  # Calculate correlation
  if (calc_cor) {
    if (!is.null(split_id)) {
      counts <- counts %>%
        group_by(split_id)
    }
    
    counts <- counts %>%
      mutate(
        r       = tidy(cor.test(!!sym(x), !!sym(y)))$estimate,
        r       = round(r, digits = 2),
        pval    = tidy(cor.test(!!sym(x), !!sym(y)))$p.value,
        cor_lab = str_c("r = ", r, ", p = ", format(pval, digits = 2))
      )
    
    if (lab_pos != "strip") {
      counts <- counts %>%
        mutate(
          min_x = min(!!sym(x)),
          max_x = max(!!sym(x)),
          min_y = min(!!sym(y)),
          max_y = max(!!sym(y)),
          lab_x = (max_x - min_x) * lab_pos[1] + min_x,
          lab_y = (max_y - min_y) * lab_pos[2] + min_y
        )
    }
  }
  
  # Create scatter plot
  # To add outline for each cluster create separate layers
  res <- counts %>%
    arrange(!!sym(feature))
  
  if (!is.null(pt_outline)) {
    
    if (!is.numeric(counts[[feature]])) {
      res <- res %>%
        ggplot(aes(!!sym(x), !!sym(y), color = !!sym(feature), fill = !!sym(feature)))
      
      feats <- counts[[feature]] %>%
        unique()
      
      if (!is.null(feat_levels)) {
        feats <- feat_levels[feat_levels %in% feats]
      }
      
      for (feat in feats) {
        if (is.na(feat)) {
          feat <- "NA"
        }
        
        f_counts <- counts %>%
          mutate(!!sym(feature) := replace_na(!!sym(feature), "NA")) %>%
          filter(!!sym(feature) == feat)
        
        res <- res +
          geom_point(
            data        = f_counts,
            mapping     = aes(fill = !!sym(feature)),
            size        = pt_outline,
            color       = "black",
            show.legend = F
          ) +
          geom_point(data = f_counts, size = pt_size)
      }
      
    } else {
      res <- res %>%
        ggplot(aes(!!sym(x), !!sym(y), color = !!sym(feature))) +
        geom_point(aes(fill = !!sym(feature)), size = pt_outline, color = "black", show.legend = F) +
        geom_point(size = pt_size)
    }
    
  } else {
    res <- res %>%
      ggplot(aes(!!sym(x), !!sym(y), color = !!sym(feature))) +
      geom_point(size = pt_size)
  }
  
  # Add regression line
  if (lm_line) {
    res <- res +
      geom_smooth(
        method   = "lm",
        se       = F,
        color    = "black",
        size     = 0.5,
        linetype = 2
      )
  }
  
  # Add correlation coefficient label
  if (calc_cor && lab_pos != "strip") {
    res <- res +
      geom_text(
        mapping       = aes(lab_x, lab_y, label = cor_lab),
        color         = "black",
        size          = lab_size,
        check_overlap = T, 
        show.legend   = F
      )
  }
  
  # Set feature colors
  if (!is.null(plot_cols)) {
    if (is.numeric(counts[[feature]])) {
      res <- res +
        scale_color_gradientn(colors = plot_cols, na.value = na_value)
      
    } else {
      res <- res +
        scale_color_manual(values = plot_cols, na.value = na_value) +
        scale_fill_manual(values = plot_cols, na.value = na_value)
    }
  }
  
  # Split plot into facets
  cor_labeller <- function(labels) {
    labels %>%
      map(~ {
        cor_labs <- counts %>%
          ungroup() %>%
          select(!!sym(feature), cor_lab) %>%
          unique()
        
        cor_labs <- set_names(cor_labs$cor_lab, cor_labs[[feature]])
        
        str_c(.x, "\n", cor_labs[.x])
      })
  }
  
  if (!is.null(split_id)) {
    if (length(split_id) == 1) {
      
      if (calc_cor && lab_pos == "strip") {
        my_labs <- cor_labeller
        
      } else {
        my_labs <- "label_value"
      }
      
      res <- res +
        facet_wrap(~ split_id, labeller = my_labs, ...)
      
    } else if (length(split_id) == 2) {
      eq <- str_c(split_id[1], " ~ ", split_id[2])
      
      res <- res +
        facet_grid(as.formula(eq), ...)
    }
  }
  
  res
}


#' Plot clonotype abundance
#' 
#' @param sobj_in Seurat object containing VDJ data
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param abundance_col meta.data column containing pre-calculated abundances
#' @param plot_colors Character vector containing colors to use for plotting
#' @param plot_levels Character vector containing levels to use for ordering
#' @param ... Additional arguments to pass to geom_line
#' @return Seurat object with calculated clonotype abundance
#' @export
#' 
plot_abundance <- function(sobj_in, clonotype_col = "clonotype_id", cluster_col = NULL,
                           abundance_col = NULL, plot_colors = NULL, plot_levels = NULL, ...) {
  
  # Calculate clonotype abundance
  if (is.null(abundance_col)) {
    sobj_in <- djvdj::calc_abundance(
      sobj_in, 
      clonotype_col = clonotype_col, 
      cluster_col   = cluster_col,
      prefix        = "."
    )
    
    abundance_col <- ".clone_abund"
  }
  
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::filter(meta_df, !is.na(!!dplyr::sym(clonotype_col)))
  
  # Rank by abundance
  if (!is.null(cluster_col)) {
    if (!is.null(plot_levels)) {
      meta_df <- dplyr::mutate(
        meta_df, 
        !!dplyr::sym(cluster_col) := factor(
          !!dplyr::sym(cluster_col), 
          levels = plot_levels
        )
      )
    }
    
    meta_df <- dplyr::group_by(meta_df, !!dplyr::sym(cluster_col))
  }
  
  meta_df <- dplyr::mutate(
    meta_df,
    rank = dplyr::row_number(!!dplyr::sym(abundance_col)),
    rank = (rank - (max(rank) + 1)) * -1
  )
  
  # Plot abundance vs rank
  gg <- ggplot2::ggplot(meta_df, aes(rank, !!dplyr::sym(abundance_col))) +
    labs(y = "abundance")
  
  if (is.null(cluster_col)) {
    gg <- gg +
      ggplot2::geom_line(...)
    
  } else {
    gg <- gg +
      ggplot2::geom_line(aes(color = !!dplyr::sym(cluster_col)), ...)
  }
  
  if (!is.null(plot_colors)) {
    gg <- gg +
      scale_color_manual(values = plot_colors)
  }
  
  gg
}






