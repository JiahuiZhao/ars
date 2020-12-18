# Test functions in the zk.R file.

test_that("Test that the number of intersection points is
          always just one more than Tk", {

  g <- dnorm
  h <- function(x) log(g(x))
  Tk <- c(-1,1)

  expect_equal(length(calc_z(c(-Inf, Inf), Tk, h)), length(Tk)+1)
})
