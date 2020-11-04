library(testthat)
library(dplyr)
library(purrr)
library(djvdj)

# Helper to test all combinations of provided arguments
test_all_args <- function(lst, .fn, ttl, chk_fn = NULL, chk_arg = NULL) {
  lst <- lst %>%
    expand.grid() %>%
    mutate(across(where(is.factor), as.character))

  n <- 1

  pwalk(lst, ~ {
    test_that(paste(ttl, n), {
      res <- .fn(...)

      if (!is.null(chk_fn)) {
        if (!is.null(chk_arg)) {
          return(chk_fn(res, chk_arg))
        }

        chk_fn(res)
      }
    })

    n <<- n + 1
  })
}

test_check("djvdj")
