library(testthat)
library(djvdj)

check_args <- function(.fn, ..., n, print_n = FALSE) {
  if (print_n) {
    print(n)
  }

  .fn(...)
}

test_check("djvdj")
