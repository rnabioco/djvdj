library(testthat)
library(dplyr)
library(tidyr)
library(purrr)
library(djvdj)

# Helper to test all combinations of provided arguments
test_all_args <- function(arg_lst, .fn, ttl, chk_fn, chk_arg = NULL) {
  arg_lst <- arg_lst %>%
    expand.grid(stringsAsFactors = FALSE)

  n <- 1

  pwalk(arg_lst, ~ {
    test_that(paste(ttl, n), {

      if (!is.null(chk_arg)) {
        return(chk_fn(.fn(...), chk_arg))
      }

      chk_fn(.fn(...))
    })

    n <<- n + 1
  })
}

test_check("djvdj")
