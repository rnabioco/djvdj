# Test data
df_1 <- vdj_so@meta.data

df_2 <- vdj_so@meta.data |>
  as_tibble(rownames = ".cell_id")

# Check all calc_gene_usage arguments
arg_lst <- list(
  input       = list(vdj_so, vdj_sce, df_1, df_2),
  data_cols   = list("v_gene", "j_gene", c("v_gene", "j_gene")),
  cluster_col = list(NULL, "seurat_clusters"),
  chain       = list(NULL, "IGH", "IGL", "IGK"),
  chain_col   = "chains"
)

test_all_args(
  arg_lst = arg_lst,
  .fn     = calc_gene_usage,
  desc    = "calc_usage args",
  chk     = expr(expect_s3_class(.res, "tbl"))
)

# Check calc_gene_usage calculations
.check_gene_usage <- function(input, genes, chain = NULL) {

  # calc_gene_usage results
  res <- calc_gene_usage(
    input,
    data_cols = genes,
    chain = chain
  ) |>
    rowwise() |>
    mutate(nm = stringr::str_c(!!!syms(genes))) |>
    arrange(desc(freq), nm)

  freq <- purrr::set_names(res$freq, res$nm)
  pct  <- purrr::set_names(res$pct, res$nm)

  # data for manual calc
  dat <- vdj_so@meta.data |>
    filter(!is.na(clonotype_id))

  gns <- genes

  if (!is.null(chain)) {
    gns <- c("chains", genes)
  }

  x <- gns |>
    purrr::set_names() |>
    map(~ {
      dat[[.x]] |>
        map(stringr::str_split, ";") |>
        purrr::flatten()
    })

  # Get total n cells
  n_cells <- map_int(x, length) |>
    unique()

  stopifnot(length(n_cells) == 1)

  # Filter chains
  if (!is.null(chain)) {
    chns <- x$chains

    chns <- chns |>
      map(~ .x %in% chain)

    x <- x[names(x) != "chains"]

    x <- x |>
      map(~ purrr::map2(.x, chns, ~ .x[.y])) |>
      map(~ map(.x, ~ {
        if (purrr::is_empty(.x)) .x <- "None"
        .x
      }))
  }

  # Calc frequency
  x <- x |>
    purrr::reduce(~ {
      stopifnot(length(.x) == length(.y))
      purrr::map2(.x, .y, ~ stringr::str_c(.x, .y))
    }) |>
    map(unique) |>
    purrr::reduce(c) |>
    table() |>
    base::sort(decreasing = TRUE)

  # Check results
  x <- purrr::set_names(
    as.integer(x),
    names(x)
  )

  expect_identical(x, freq)
  expect_identical((x / n_cells) * 100, pct)
}

test_that("calc_gene_usage check calcs", {

  vdj_genes <- c(
    "v_gene", "d_gene",
    "j_gene", "c_gene"
  )

  vdj_genes <- seq_along(vdj_genes) |>
    map(combn, x = vdj_genes, simplify = FALSE) |>
    purrr::flatten()

  ins <- list(vdj_so, vdj_sce, vdj_so@meta.data)

  purrr::walk(ins, ~ {
    vdj_genes |>
      purrr::walk(.check_gene_usage, input = .x)
  })

  purrr::walk(ins, ~ {
    vdj_genes |>
      purrr::walk(.check_gene_usage, input = .x, chain = "IGH")
  })

  purrr::walk(ins, ~ {
    vdj_genes |>
      purrr::walk(.check_gene_usage, input = .x, chain = c("IGH", "IGK"))
  })
})

# # Bad vdj gene
# test_that("calc_gene_usage bad gene", {
#
#   ins <- list(vdj_so, vdj_sce, vdj_so@meta.data)
#
#   purrr::walk(ins, ~ {
#     expect_error(
#       calc_gene_usage(.x, "BAD"),
#       "Column .+ doesn't exist"
#     )
#   })
# })

