# Test sample_sk function in sk.R file.

test_that("Return NA if length(Tk) is greater than length(z), the number of intervals", {
  expect_true(is.na(sample_sk(Tk=c(-1,0, 1), z=c(-Inf, 0, Inf), h, exp_uks)))
  expect_false(is.na(sample_sk(Tk=c(-1, 1), z=c(-Inf, 0, Inf), h, exp_uks)))
})
  