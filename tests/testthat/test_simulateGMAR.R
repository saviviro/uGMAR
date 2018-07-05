library(uGMAR)
context("simulate GMAR")

params13gs <- c(0, 0.9, 1, -1, 0.5, 0.8, 1, -0.5, 0.5, 0.3, 0.4, 5, 7) # M1=1, M2=2
params12 <- c(0.8, 0.5, 0.5, 0.3, 0.7, 0.1, 0.6, 10, 12) # StMAR
params23 <- c(1, 0.1, 0.1, 0.1, 2, -0.2, -0.2, 0.2, 3, 0.3, 0.3, 0.3, 0.6, 0.3)
set.seed(1); data12 <- simulateGMAR(1, 2, params12, StMAR=TRUE, nsimu=3)
set.seed(2); data23 <- simulateGMAR(2, 3, StMAR=FALSE, params23, nsimu=1)
set.seed(3); data23init <- simulateGMAR(2, 3, StMAR=FALSE, params23, nsimu=2, initvalues=c(1, 1.2))
set.seed(1); data13gs <- simulateGMAR(1, c(1, 2), params13gs, GStMAR=TRUE, nsimu=3)
sample12 <- data12$sample
component12 <- data12$component
sample23 <- data23$sample
component23 <- data23$component
weights23 <- data23$mixingWeights
sample23init <- data23init$sample

test_that("simulateGMAR simulates correctly from non-restricted process", {
  expect_equal(data13gs$sample[1], -0.3833485, tolerance=1e-6)
  expect_equal(data13gs$mixingWeights[3,3], 0.00969763, tolerance=1e-6)
  expect_equal(sample12[3], 2.655959, tolerance=1e-6)
  expect_equal(component12[3], 1)
  expect_equal(sample23, 1.017869, tolerance=1e-6)
  expect_equal(component23, 2)
  expect_equal(weights23[,2], 0.2067762, tolerance=1e-6)
  expect_equal(sample23init[2], 1.351309, tolerance=1e-6)
})

params12r <- c(1, 2, -0.9, 1, 2, 0.7)
params23r <- c(-0.1, -0.2, -0.3, -0.2, 0.5, 0.1, 0.2, 0.3, 0.5, 0.3, 3, 14, 50) # StMAR
params22gsr <- c(1, 1.2, 0.3, 0.3, 1, 1.2, 0.5, 6)
set.seed(1); data12r <- simulateGMAR(1, 2, params12r, StMAR=FALSE, nsimu=1, restricted=TRUE)
set.seed(2); data12rinit <- simulateGMAR(1, 2, params12r, StMAR=FALSE, nsimu=2, initvalues=c(1), restricted=TRUE)
set.seed(3); data23r <- simulateGMAR(2, 3, params23r, StMAR=TRUE, nsimu=2, restricted=TRUE)
set.seed(1); data22gsr <- simulateGMAR(2, c(1, 1), params22gsr, GStMAR=TRUE, restricted=TRUE, nsimu=2)
sample12r <- data12r$sample
component12r <- data12r$component
weights12r <- data12r$mixingWeights
sample12init <- data12rinit$sample
sample23r <- data23r$sample
component23r <- data23r$component

test_that("simulateGMAR simulates correctly from restricted process", {
  expect_equal(data22gsr$sample[2], 2.583507, tolerance=1e-6)
  expect_equal(data22gsr$component[1], 2)
  expect_equal(sample12r, 1.018146, tolerance=1e-6)
  expect_equal(component12r, 2)
  expect_equal(weights12r[,1],  0.7716839, tolerance=1e-6)
  expect_equal(sample12init[1], 0.6312408, tolerance=1e-6)
  expect_equal(sample23r[2], -0.5006049, tolerance=1e-6)
  expect_equal(component23r[2], 2)
})

R1 <- matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
R3 <- matrix(c(0.5, 0.5), ncol=1)
R4 <- diag(1, ncol=2, nrow=2)
params22c <- c(1, 0.1, -0.1, 1, 2, 0.2, 2, 0.8, 11, 12) # R4, R3, StMAR
params32cr <- c(1, 2, 0.3, -0.3, 1, 2, 0.6) # R1
params12gsrc <- c(0, 0, -0.99, 1, 1.2, 0.3, 10) # M1=1, M2=1, R3
set.seed(1); data22c <- simulateGMAR(2, 2, params22c, StMAR=TRUE, nsimu=1, restricted=FALSE, constraints=TRUE, R=list(R4, R3))
set.seed(2); data32cr <- simulateGMAR(3, 2, params32cr, nsimu=2, restricted=TRUE, constraints=TRUE, R=R1)
set.seed(3); data12gsrc <- simulateGMAR(2, c(1, 1), params12gsrc, GStMAR=TRUE, nsimu=2, restricted=TRUE, constraints=TRUE, R=R3)
sample22c <- data22c$sample
weights22c <- data22c$mixingWeights
sample32cr <- data32cr$sample
component32cr <- data32cr$component

test_that("simulateGMAR simulates correctly from constrained process", {
  expect_equal(data12gsrc$sample[2], -1.732949, tolerance=1e-6)
  expect_equal(data12gsrc$component[1], 2)
  expect_equal(data12gsrc$mixingWeights[,1], c(0.3122889, 0.3084366), tolerance=1e-6)
  expect_equal(sample22c,  0.1681538, tolerance=1e-6)
  expect_equal(weights22c[,1], 0.8956446, tolerance=1e-6)
  expect_equal(sample32cr[2], 0.7406432, tolerance=1e-6)
  expect_equal(component32cr[2], 1)
})
