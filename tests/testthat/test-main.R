
library(ars)

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
  t_dist <- function(x){dt(x, df = 2)}
  gamma_dist <- function(x){dgamma(x, shape = 2, rate=3)}
  f_dist <- function(x){df(x, df1 = 10, df2 = 3)}
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
  expect_gt(fderiv(dnorm, calc_init_vals(dnorm, bounds = c(-Inf, 0))[1]), 0)

  print("Checking finite lower bound to finite upper bound case")
  expect_length(calc_init_vals(beta_dist, bounds = c(0,1)), 2)
})


test_that("Check Length Consistency of z", {
  print("Z is calculated properly.")

  g = dnorm
  h = function(x) log(g(x))
  Tk = c(-1,1)
  bounds = c(-2, 2)
  z = calc_z(bounds, Tk, h)

  expect_equal(length(z), 1+length(Tk))
})

test_that("Distributions are correct", {

  ptrunc <- function(x, spec, a = -Inf, b = Inf, ...)
  {
    tt <- x
    aa <- rep(a, length(x))
    bb <- rep(b, length(x))
    G <- get(paste("p", spec, sep = ""), mode = "function")
    tt <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
    tt <- tt - G(aa, ...)
    tt <- tt/(G(bb, ...) - G(aa, ...))
    return(tt)
  }

  # the ptrunc function above is adapted from
  # https://github.com/boyinggong/Adaptive-rejection-sampling/blob/master/ars/tests/testthat/test-main.R

  print("test normal and truncated normal with mean 0 variance 1")
  expect_true(ks.test(ars(dnorm, c(-Inf,Inf), n = 1000), pnorm)$p.value > 0.01)
  expect_true(ks.test(ars(dnorm, c(-2,2), n = 1000), function(x) {ptrunc(x, "norm", a = -2, b = 2)})$p.value > 0.01)
  expect_true(ks.test(ars(dnorm, c(1,Inf), n = 1000), function(x){ptrunc(x, "norm", a = 1, b = Inf)})$p.value > 0.01)

  print("test gamma distribution with shape = 5 and rate = 3")
  gamma_dist <- function(x){dgamma(x, shape = 5, rate=3)}
  expect_true(ks.test(ars(g=gamma_dist, n=1000, initial = NULL, bounds = c(0, Inf)),
                      function(x) {ptrunc(x, "gamma", a = 0, b = Inf, shape = 5, rate=3)})$p.value > 0.01)
})
