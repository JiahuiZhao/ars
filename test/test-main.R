# Test fderiv function.
test_that("Check cases for different distributions.", {
  print("The derivatives are being properly calculated for the input distributions.")
  expect_gt(fderiv(dnorm, -1e-10), 0)
  expect_lt(fderiv(dnorm, 1e-10), 0)
  expect_equal(fderiv(dt, 0, df=1), 0)
})

# Test calc_init_vals function.

test_that("Check that the derivatives for initial points are correct for the unbounded cases", {
  
  norm_dist <- function(x){dnorm(x, mean=5)}
  t_dist <- function(x){dt(x, df = 4)}
  gamma_dist <- function(x){dgamma(x, shape = 5, rate=3)}
  f_dist <- function(x){df(x, df1 = 10, df2 = 3)}
  trunc_norm_dist <- function(x){dtruncnorm(x, a=-Inf, b=0)}
  beta_dist <- function(x){dbeta(x, shape1=2, shape2=3)}
  
  print("Checking negative infinity to infinity case")
  expect_gt(fderiv(norm_dist, calc_init_vals(norm_dist, bounds = c(-Inf, Inf))[1]), 0)
  expect_lt(fderiv(norm_dist, calc_init_vals(norm_dist, bounds = c(-Inf, Inf))[2]), 0)
  expect_gt(fderiv(t_dist, calc_init_vals(t_dist, bounds = c(-Inf, Inf))[1]), 0)
  expect_lt(fderiv(t_dist, calc_init_vals(t_dist, bounds = c(-Inf, Inf))[2]), 0)
  
  print("Checking finite lower bound to infinity case")
  expect_lt(fderiv(gamma_dist, calc_init_vals(gamma_dist, bounds = c(0, Inf))[2]), 0)
  expect_lt(fderiv(f_dist, calc_init_vals(f_dist, bounds = c(0, Inf))[2]), 0)
  
  print("Checking negative infinity to finite upper bound case")
  expect_gt(fderiv(trunc_norm_dist, calc_init_vals(trunc_norm_dist, bounds = c(-Inf, 0))[1]), 0)
  
  print("Checking finite lower bound to finite upper bound case")
  expect_length(calc_init_vals(beta_dist, bounds = c(0,1)), 2)
})



