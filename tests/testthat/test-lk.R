# Test functions in the lk.R file.

# Test calc_lk function.
test_that("Check class and length of calc_lk outputs", {
  print("Checking the class and length of results returned by calc_lk")
  
  g <- dnorm
  h <- function(x) log(g(x))
  Tk <- c(-1,1)
  
  expect_equal(class(calc_lk(h, Tk)), "list")
  expect_equal(class((calc_lk(h, Tk))[[1]]), "function")
  expect_equal(length(calc_lk(h, Tk)), length(Tk)-1)
})

# Test get_lk_x function
test_that("Check if the sampled x star value is outside the abscissae returns -Inf", {
  
  lks <- calc_lk(h, Tk)
  
  expect_equal(get_lk_x(x=-10, Tk=c(-1,1), lks),-Inf)
  expect_false(get_lk_x(x=0, Tk=c(-1,1), lks) == -Inf)
  expect_equal(class(get_lk_x(x=0, Tk=c(-1,1), lks)), "numeric")
})
