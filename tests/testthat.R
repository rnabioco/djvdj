library(testthat)
library(djvdj)

check_args <- function(.fn, ..., n = NULL, print_n = FALSE) {
  if (print_n && !is.null(n)) {
    print(n)
  }

  .fn(...)
}

test_check("djvdj")
