context("Wald and LR tests")
library(uGMAR)

## A(M)(p)_(p)(M)(d)

# Structural GMVAR(2, 2), d=2 model identified with sign-constraints:
params12 <- c(1.7041595, 0.8492518, 0.3006342, 4.1194618, 0.7271892, 1.9817267, 0.6331898)
mod12 <- GSMAR(simudata, p=1, M=2, params=params12, model="GMAR")

A1 <- rbind(c(0, 1, 0, 0, -1, 0, 0),
           c(1, 0, 0, -1, 0, 0, 0))
c1 <- c(0.1, -2.4)
wald1 <- Wald_test(mod12, A1, c1)

test_that("Wald_test works correctly", {
  expect_equal(wald1$df, 2)
  expect_equal(wald1$test_stat, 1.32194, tolerance=1e-1)
  expect_equal(wald1$p_value, 0.5163504, tolerance=1e-1)
})


# The same model as above but with the AR parameters restricted to be the
# same in both regimes.
params12r <- c(2.1849744, 2.9342138, 0.8052687, 0.2948068, 1.9019095, 0.6437754)
mod12r <- GSMAR(simudata, p=1, M=2, params=params12r, model="GMAR", restricted=TRUE)

lr1 <- LR_test(mod12, mod12r)

test_that("LR_test works correctly", {
  expect_equal(lr1$df, 1)
  expect_equal(lr1$test_stat, 1.378559, tolerance=1e-4)
  expect_equal(lr1$p_value, 0.2403468, tolerance=1e-4)
})
