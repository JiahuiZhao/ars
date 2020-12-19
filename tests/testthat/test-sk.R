
library(pracma)

# Test sample_sk function in sk.R file.

test_that("Return error if the number of abscisscae is not one less than the length of Z, the number of intervals,
          the number of intervals", {

  g <- dnorm
  h <- function(x) log(g(x))

  expect_error(sample_sk(Tk=c(-1,0,1), z=c(-Inf, 0, Inf), h, exp_uks))
  expect_error(sample_sk(Tk=c(-1, 0, 1, 2), z=c(-Inf, 0, 4, Inf), h, exp_uks))
  expect_equal(class(sample_sk(Tk=c(-1, 1), z=c(-Inf, 0, Inf), h, exp_uks)), "numeric")
})
