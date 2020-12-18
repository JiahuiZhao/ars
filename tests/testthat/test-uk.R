# Test functions in the uk.R file.

# Test calc_uk function.
test_that("Check class and length of calc_uk outputs", {
  print("Checking classes of results returned by calc_uk")

  g <- dnorm
  h <- function(x) log(g(x))
  Tk <- c(-2,0,2)

  expect_equal(class(calc_uk(h, Tk, F)), "list")
  expect_equal(class((calc_uk(h, Tk, F))[[1]]), "function")
  expect_equal(class((calc_uk(h, Tk, F))[[2]]), "function")
  expect_equal(length(calc_uk(h, Tk, F)), length(Tk))
})

# Test get_uk_x function.
test_that("Throw error if domain specification of z is not valid for the given x.", {

  print("Expect error if domain doesn't include the given x value")
  Tk = c(-2,0,2)
  h = dnorm
  uks = calc_uk(h, Tk, F)

  expect_error(get_uk_x(x=-15, z=c(-10,0,10), uks))
  expect_error(get_uk_x(x=15, z=c(-10,0,10), uks))
  expect_equal(class(get_uk_x(x=5, z=c(-10,0,10), uks)), "numeric")
})

