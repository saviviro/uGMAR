library(uGMAR)
context("parameter reforms")

params11 <- c(1, 0.9, 1, 10) # StMAR
params12 <- c(0.8, 0.5, 0.5, 2, -1, 0.1, 0.6)
params12_2 <- c(2, -1, 0.1, 0.8, 0.5, 0.5, 0.4, 12, 30) #StMAR
params22 <- c(0.2, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.49)
params13 <- c(0.1, 0.99, 0.1, 0.2, -0.99, 0.2, 0.3, 0.01, 0.3, 0.5, 0.5)
params23 <- c(0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 0.5, 0.2, 0.3, 0.3, 0.3, 0.3, 0.8, 0.05, 11, 12, 13) # StMAR
params23_2 <- c(0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 0.5, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.05, 11, 12, 13) # StMAR

params12r <- c(0.1, 0.1, 1, 0.1, 0.1, 0.6, 11, 12) # StMAR
params22r <- c(0.1, 0.2, 0.99, 0.01, 0.1, 0.2, 0.05)
params23r <- c(0.1, 0.3, 0.4, -0.4, 0.3, 1, 2, 3, 0.5, 0.1, 100, 112, 130) # StMAR
params13r <- c(1, 2, 3, 0.99999, 1, 2, 3, 0.2, 0.15)
params23r2 <- c(0.1, 0.2, 0.3, 0.9, 0.2, 0.1, 0.2, 0.3, 0.3, 0.5)

ref11 <- reformParameters(1, 1, params=params11, StMAR=TRUE)
ref12 <- reformParameters(1, 2, params=params12)
ref12_2 <- reformParameters(1, 2, params=params12_2, StMAR=TRUE)
ref13 <- reformParameters(1, 3, params=params13)
ref23 <- reformParameters(2, 3, params=params23, StMAR=TRUE)

ref12r <- reformParameters(1, 2, params=params12r, StMAR=TRUE, restricted=TRUE)
ref23r <- reformParameters(2, 3, params=params23r, StMAR=TRUE, restricted=TRUE)
ref13r <- reformParameters(1, 3, params=params13r, restricted=TRUE)

test_that("reformParameters works correctly", {
  expect_equal(ref11$pars[2, 1], 0.9)
  expect_equal(ref11$alphas, 1)
  expect_equal(ref11$dfs[1], 10)
  expect_equal(ref12$pars[3, 2], 0.1)
  expect_equal(ref12$alphas[2], 0.4)
  expect_equal(ref12_2$pars[1, 1], 2)
  expect_equal(ref12_2$alphas[2], 0.6)
  expect_equal(ref13$pars[2, 3], 0.01)
  expect_equal(ref13$alphas[1], 0.5)
  expect_equal(ref23$pars[2, 2], 0.5)
  expect_equal(ref23$pars[3, 3], 0.3)
  expect_equal(ref23$dfs[3], 13)

  expect_equal(ref12r$params[5], 1)
  expect_equal(ref12r$pars[3, 1], 0.1)
  expect_equal(ref12r$dfs[2], 12)
  expect_equal(ref23r$params[4], 1)
  expect_equal(ref23r$pars[3, 3], 0.3)
  expect_equal(ref23r$alphas[3], 0.4)
  expect_equal(ref23r$dfs[3], 130)
  expect_equal(ref13r$params[3], 1)
  expect_equal(ref13r$params[8], 0.99999)
  expect_equal(ref13r$pars[2, 3], 0.99999)
  expect_equal(ref13r$pars[1, 3], 3)
  expect_equal(ref13r$pars[1, 2], 2)
  expect_equal(ref13r$alphas[1], 0.2)
  expect_equal(ref13r$alphas[3], 0.65)
})

R1 <- matrix(c(1, 0, 0, 0, 0, 1), ncol=2)
R2 <- diag(1, ncol=3, nrow=3)
R3 <- matrix(c(0.5, 0.5), ncol=1)
R4 <- diag(1, ncol=2, nrow=2)

params32c <- c(1, 0.1, -0.1, 1, 2, 0.2, -0.2, 2, 0.6, 11, 12)
refc32c <- reformConstrainedPars(3, 2, params=params32c, StMAR=TRUE, R=list(R1, R1))

params33c <- c(1, 0.1, 0.1, 0.1, 1, 2, 0.2, 0.2, 0.2, 2, 3, 0.3, -0.3, 3, 0.5, 0.4)
refc33c <- reformConstrainedPars(3, 3, params=params33c, R=list(R2, R2, R1))

params21c <- c(1, 1, 1, 3)
refc21c <- reformConstrainedPars(2, 1, params=params21c, StMAR=TRUE, R=list(R3))

params22c <- c(1, 0.1, -0.1, 1, 2, 0.2, 2, 0.8, 11, 12)
refc22c <- reformConstrainedPars(2, 2, params=params22c, StMAR=TRUE, R=list(R4, R3))

params21cr <- c(1, 1, 1)
refc21cr <- reformConstrainedPars(2, 1, params=params21cr, restricted=TRUE, R=R3)

params22cr <- c(1, 2, 0.8, 1, 2, 0.7, 11, 12)
refc22cr <- reformConstrainedPars(2, 2, params=params22cr, StMAR=TRUE, restricted=TRUE, R=R3)

params32cr <- c(1, 2, 0.3, -0.3, 1, 2, 0.6)
refc32cr <- reformConstrainedPars(3, 2, params=params32cr, restricted=TRUE, R=R1)

test_that("reformConstrainedPars works correctly", {
  expect_equal(refc32c[1], 1)
  expect_equal(refc32c[3], 0)
  expect_equal(refc32c[8], 0)
  expect_equal(refc32c[10], 2)
  expect_equal(refc32c[13], 12)
  expect_equal(refc33c[4], 0.1)
  expect_equal(refc33c[6], 2)
  expect_equal(refc33c[13], 0)
  expect_equal(refc33c[14], -0.3)
  expect_equal(refc33c[17], 0.4)
  expect_equal(refc21c[2], 0.5)
  expect_equal(refc21c[3], 0.5)
  expect_equal(refc22c[3], -0.1)
  expect_equal(refc22c[7], 0.1)
  expect_equal(refc22c[11], 12)

  expect_equal(refc21cr[2], 0.5)
  expect_equal(refc21cr[3], 0.5)
  expect_equal(refc22cr[3], 0.4)
  expect_equal(refc22cr[7], 0.7)
  expect_equal(refc22cr[9], 12)
  expect_equal(refc32cr[1], 1)
  expect_equal(refc32cr[3], 0.3)
  expect_equal(refc32cr[4], 0)
  expect_equal(refc32cr[5], -0.3)
  expect_equal(refc32cr[7], 2)
  expect_equal(refc32cr[8], 0.6)
})

test_that("sortComponents sorts correctly", {
 expect_equal(sortComponents(1, 1, params11, StMAR=TRUE), params11)
 expect_equal(sortComponents(1, 2, params12), params12)
 expect_equal(sortComponents(1, 2, params12_2, StMAR=TRUE), c(0.8, 0.5, 0.5, 2, -1, 0.1, 0.6, 30, 12))
 expect_equal(sortComponents(2, 2, params22), c(0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.2, 0.2, 0.51))
 expect_equal(sortComponents(1, 3, params13), params13)
 expect_equal(sortComponents(2, 3, params23, StMAR=TRUE), c(0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3, 0.2, 0.5, 0.5, 0.2, 0.8, 0.15, 11, 13, 12))
 expect_equal(sortComponents(2, 3, params23_2, StMAR=TRUE), c(0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1, 0.2, 0.5, 0.5, 0.2, 0.65, 0.3, 13, 11, 12))

 expect_equal(sortComponents(1, 2, params12r, restricted=TRUE, StMAR=TRUE), params12r)
 expect_equal(sortComponents(2, 2, params22r, restricted=TRUE), c(0.2, 0.1, 0.99, 0.01, 0.2, 0.1, 0.95))
 expect_equal(sortComponents(2, 3, params23r, restricted=TRUE, StMAR=TRUE), c(0.1, 0.4, 0.3, -0.4, 0.3, 1, 3, 2, 0.5, 0.4, 100, 130, 112))
 expect_equal(sortComponents(1, 3, params13r, restricted=TRUE), c(3, 1, 2, 0.99999, 3, 1, 2, 0.65, 0.2))
 expect_equal(sortComponents(2, 3, params23r2, restricted=TRUE), c(0.2, 0.1, 0.3, 0.9, 0.2, 0.2, 0.1, 0.3, 0.5, 0.3))
})
