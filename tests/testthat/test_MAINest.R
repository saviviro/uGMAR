library(uGMAR)
context("functions in MAINest")

# Note that fitGSMAR is not covered at all
test_that("get_minval works correctly", {
  expect_equal(get_minval(VIX), -9999)
  expect_equal(get_minval(rep(0, 1000)), -9999)
  expect_equal(get_minval(rep(0, 1001)), -99999)
})

params12 <- c(0.18, 0.93, 0.01, 0.86, 0.68, 0.02, 0.88)
gmar12 <- GSMAR(data=logVIX, p=1, M=2, params=params12, model="GMAR")

params11t <- c(0.18, 0.93, 0.86, 4)
stmar11 <- GSMAR(data=logVIX, p=1, M=1, params=params11t, model="StMAR")

params22gs <- c(2.44, 0.06, 0.01, 0.017, 0.15, 0.93, 0.01, 0.002, 0.04, 6)
gstmar22 <- GSMAR(data=logVIX, p=2, M=c(1, 1), params=params22gs, model="G-StMAR")


test_that("get_minval works correctly", {
  gmar12it <- iterate_more(gmar12, calc_std_errors=FALSE)
  stmar11it <- iterate_more(stmar11, calc_std_errors=FALSE)
  gstmar22it <- iterate_more(gstmar22, calc_std_errors=FALSE)

  expect_equal(gmar12it$params, c(0.182979876, 0.926503609, 0.002144505, 0.855628637, 0.682666342, 0.019045056, 0.883703986), tolerance=1e-4)
  expect_equal(stmar11it$params, c(0.144624407, 0.940474160, 0.007601943, 2.627557356), tolerance=1e-4)
  expect_equal(gstmar22it$params, c(3.639637340, 0.234833577, -0.656612282, 0.007346265, 0.208349151, 0.917842454, -0.002236792,
                                    0.003150485, 0.061418991, 6.008131952), tolerance=1e-4)
})

