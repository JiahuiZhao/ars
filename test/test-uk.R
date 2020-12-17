# Test functions in the uk.R file.

# Testing calc_uk function.
test_that("Check class and length of calc_uk outputs", {
  print("Checking classes of results returned by calc_uk")
  
  g <- dnorm
  h <- function(x) log(g(x))
  Tk <- c(-1,1)
  
  expect_equal(class(calc_uk(h, Tk, F)), "list")
  expect_equal(class((calc_uk(h, Tk, F))[[1]]), "function")
  expect_equal(class((calc_uk(h, Tk, F))[[2]]), "function")
  expect_equal(length(calc_uk(h, Tk, F)), length(Tk))
})

# Testing get_uk_x function.
test_that("Throw error if domain specification of z is not valid for the given x.", {
  expect_error(get_uk_x(x=-15, z=c(-10,0,10), uks))
  expect_error(get_uk_x(x=15, z=c(-10,0,10), uks))
})

