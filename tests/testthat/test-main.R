library(pracma)
library(truncnorm)

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


# get_uk_x
# get_lk_x

#################  z  #########################

test_that("Check Length Consistency of z", {
  print("Z is calculated properly.")
  expect_equal(length(z), 1+length(Tk))
})

test_that("General Functions are correct", {

  print("test normal and truncated normal with mean 0 variance 1")
  expect_true(ks.test(ars(dnorm, c(-Inf,Inf), n = 1000), pnorm)$p.value > 0.05)
  expect_true(ks.test(ars(dnorm, c(-1,1), n = 1000), function(x) {ptrunc(x, "norm", a = -1, b = 1)})$p.value > 0.05)
  expect_true(ks.test(ars(dnorm, c(-Inf,1), n = 1000), function(x){ptrunc(x, "norm", a = -Inf, b = 1)})$p.value > 0.05)
  expect_true(ks.test(ars(dnorm, c(1,Inf), n = 1000), function(x){ptrunc(x, "norm", a = 1, b = Inf)})$p.value > 0.05)

  print("test Chi square and truncated Chi square with df > 2")
  expect_true(ks.test(ars(dchisq, c(0,Inf), n = 1000, df = 3), function(x) {ptrunc(x, "chisq", a = 0, b = Inf, df=3)})$p.value > 0.05)
  expect_true(ks.test(ars(dchisq, c(1,4), n = 1000, df = 3), function(x) {ptrunc(x, "chisq", a = 1, b = 4, df=3)})$p.value > 0.05)

  print("test Beta and truncated Beta with df1 > 1 and df2 > 1")
  expect_true(ks.test(ars(dbeta, c(0,1), n = 1000, shape1 = 2, shape2 = 2), function(x) {pbeta(x, 2, 2)})$p.value > 0.05)
  expect_true(ks.test(ars(dbeta, c(0,1/2), n = 1000, shape1 = 2, shape2 = 2), function(x) {ptrunc(x, "beta", a = 0, b = 1/2, 2, 2)})$p.value > 0.05)
  expect_true(ks.test(ars(dbeta, c(1/4,3/4), n = 1000, shape1 = 2, shape2 = 2), function(x) {ptrunc(x, "beta", a = 1/4, b = 3/4, 2, 2)})$p.value > 0.05)
})
