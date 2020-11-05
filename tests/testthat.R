library(testthat)
library(dplyr)
library(tidyr)
library(purrr)
library(djvdj)

# Helper to test all combinations of provided arguments
test_all_args <- function(lst, .fn, ttl, chk_fn, chk_arg = NULL) {
  lst <- lst %>%
    expand.grid(stringsAsFactors = FALSE)

  n <- 1

  pwalk(lst, ~ {
    test_that(paste(ttl, n), {

      res <- .fn(...)

      if (!is.null(chk_arg)) {
        return(chk_fn(res, chk_arg))
      }

      chk_fn(res)
    })

    n <<- n + 1
  })
}

test_check("djvdj")
