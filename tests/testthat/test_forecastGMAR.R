library(uGMAR)
context("forecastGMAR")

params12 <- c(1.0, 0.9, 0.25, 4.5, 0.7, 3.0, 0.8)
params23 <- c(2.7, 0.8, -0.06, 0.3, 3.5, 0.8, -0.07, 2.6, 7.2, 0.3, -0.01, 0.1, 0.6, 0.25)
params12r <- c(1.4, 1.8, 0.9, 0.3, 3.3, 0.8)
params23r <- c(1.7, 1.9, 2.1, 0.8, -0.05, 0.3, 0.7, 4.5, 0.7, 0.2)
R1 <- matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
R3 <- matrix(c(0.5, 0.5), ncol=1)
R4 <- diag(1, ncol=2, nrow=2)
params22c <- c(1, 0.1, -0.1, 1, 2, 0.2, 2, 0.8, 11, 12) # R4, R3, StMAR
params32cr <- c(1, 2, 0.3, -0.3, 1, 2, 0.6) # R1

set.seed(1); pred12 <- forecastGMAR(VIX, 1, 2, params12, nsteps=1, nsimu=50, printRes=FALSE, plotRes=FALSE, useMean=TRUE)
set.seed(2); pred23 <- forecastGMAR(VIX, 2, 3, params23, nsteps=3, nsimu=50, conflevel=c(0.99, 0.90, 0.60),
                                    printRes=FALSE, plotRes=FALSE)
set.seed(3); pred12r <- forecastGMAR(VIX, 1, 2, params12r, nsteps=3, nsimu=50, conflevel=c(0.999, 0.001),
                                     restricted=TRUE, printRes=FALSE, plotRes=FALSE)
set.seed(4); pred23r <- forecastGMAR(VIX, 2, 3, params23r, nsteps=1, nsimu=50, conflevel=c(0.5),
                                     restricted=TRUE, printRes=FALSE, plotRes=FALSE, useMean=TRUE)
set.seed(5); pred22c <- forecastGMAR(VIX, 2, 2, params22c, StMAR=TRUE, constraints=TRUE, R=list(R4, R3), nsteps=1, nsimu=20,
                                    conflevel=c(0.99, 0.90, 0.60), printRes=FALSE, plotRes=FALSE)
set.seed(6); pred32cr <- forecastGMAR(VIX, 3, 2, params32cr, restricted=TRUE, constraints=TRUE, R=R1, nsteps=1, nsimu=10,
                                      conflevel=c(0.90), printRes=FALSE, plotRes=FALSE, useMean=TRUE)

test_that("forecastGMAR gives correct best prediction", {
  expect_equal(pred12$MEAN, 11.15113, tolerance=1e-3)
  expect_equal(pred23$`0.5`[3], 10.77946, tolerance=1e-3)
  expect_equal(pred12r$`0.5`[2], 11.94374, tolerance=1e-3)
  expect_equal(pred23r$MEAN, 10.41577, tolerance=1e-3)
  expect_equal(pred22c$`0.5`, 5.702248, tolerance=1e-3)
  expect_equal(pred32cr$MEAN, 2.112541, tolerance=1e-3)
})

test_that("forecastGMAR gives correct confidence intervals", {
  expect_equal(pred12$`0.025`, 10.1512, tolerance=1e-3)
  expect_equal(pred12$`0.9`, 11.69211, tolerance=1e-3)
  expect_equal(pred23$`0.005`[3], 9.144938, tolerance=1e-3)
  expect_equal(pred23$`0.8`[2], 11.56006, tolerance=1e-3)
  expect_equal(pred12r$`5e-04`[3], 8.716014, tolerance=1e-3)
  expect_equal(pred12r$`0.5005`[1], 11.58921, tolerance=1e-3)
  expect_equal(pred23r$`0.75`, 11.14477, tolerance=1e-3)
  expect_equal(pred22c$`0.2`, 2.971152, tolerance=1e-3)
  expect_equal(pred32cr$`0.95`, 3.680976, tolerance=1e-3)
})

params12gs <- c(3.98, 0.68, 0.36, 0.70, 0.94, 11.75, 0.25, 2.03)
params23gs <- c(2.0, 0.83, 0.01, 0.36, 1.14, 0.90, 0.01, 0.06, 4.23, 0.72, 0.01, 3.85, 0.6, 0.20, 3.3)
pred12gs <- forecastGMAR(VIX, 1, c(1, 1), params12gs, GStMAR=TRUE, oneStepCondMean=TRUE, plotRes=FALSE, printRes=FALSE)
pred23gs <- forecastGMAR(VIX, 2, c(2, 1), params23gs, GStMAR=TRUE, oneStepCondMean=TRUE, plotRes=FALSE, printRes=FALSE)

test_that("forecastGMAR oneStepCondMean works correctly", {
  expect_equal(pred12gs, 11.3083, tolerance=1e-4)
  expect_equal(pred23gs, 11.425, tolerance=1e-4)
})
