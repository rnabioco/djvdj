library(testthat)
library(dplyr)
library(tidyr)
library(purrr)
library(djvdj)

# Helper to test all combinations of provided arguments
test_all_args <- function(arg_lst, .fn, ttl, chk) {
  arg_lst <- arg_lst %>%
    expand.grid(stringsAsFactors = FALSE)

  n <- 1

  pwalk(arg_lst, ~ {
    test_that(paste(ttl, n), {

      if (is.call(chk)) {
        .res <- .fn(...)

        return(eval(chk))
      }

      chk(.fn(...))
    })

    n <<- n + 1
  })
}

test_check("djvdj")
