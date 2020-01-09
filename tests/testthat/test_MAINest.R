library(uGMAR)
context("functions in MAINest")

test_that("get_minval works correctly", {
  expect_equal(get_minval(VIX), -9999)
})
