library(uGMAR)
context("get_test_Omega")

## The results slightly differ depending on whether numerical integration (without the package "gsl")
# or hypergeometric function is used to calculate the quantile residuals. Also, the limited precission
# of floating point numbers (default precission used) might accumulate to differences big enough for
# the tests to fail with too small tolerance (this particularly happens if the parameters contain
# large degree of freedom parameters).

# get_test_Omega calls other functions to deal with constraints, so they are not tested separately
params11t <- c(0.9, 0.92, 1.01, 2.89)
params12 <- c(1.7, 0.85, 0.3, 4.12, 0.73, 1.98, 0.63)
params23t <- c(2.54, 0.99, -0.21, 0.36, 5.13, 0.9, -0.24, 1.88, 7.93, 0.2, 0.03, 0.1, 0.53, 0.36, 9, 10, 11)
params12r <- c(2.18, 2.93, 0.81, 0.29, 1.9, 0.64)
params12gs <- c(4.13, 0.73, 1.98, 1.7, 0.85, 0.3, 0.37, 9) # M1=1, M2=1
params13gsr <- c(4.8, 3.31, 3.74, 0.69, 2, 0.19, 0.41, 0.34, 0.3, 9) # M1=2, M2=1

g_norm <- function(r) cbind(r^2-1, r^3, r^4-3)
omega_n11t <- get_test_Omega(simudata, 1, 1, params11t, model="StMAR", g=g_norm, dim_g=3)
omega_n12 <- get_test_Omega(simudata, 1, 2, params12, g=g_norm, dim_g=3)
omega_n12r <- get_test_Omega(simudata, 1, 2, params12r, restricted=TRUE, g=g_norm, dim_g=3)
omega_n23t <- get_test_Omega(simudata, 2, 3, params23t, model="StMAR", g=g_norm, dim_g=3)
omega_n23gs <- get_test_Omega(simudata, 1, c(1, 1), params12gs, model="G-StMAR", g=g_norm, dim_g=3)
omega_n13gsr <- get_test_Omega(simudata, 1, c(2, 1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_norm, dim_g=3)

test_that("get_test_Omega works for normality tests", {
  expect_equal(omega_n11t[2, ], c(-0.0538172, 4.1840459, 0.2731159), tolerance=1e-3)
  expect_equal(omega_n12[1, ], c(0.2006648, -0.1384332, 1.2561043), tolerance=1e-3)
  expect_equal(omega_n12[2, 1], -0.1384332, tolerance=1e-3)
  expect_equal(omega_n12r[, 3], c(1.1721352, -0.9816178, 14.5631772), tolerance=1e-3)
  expect_equal(omega_n12r[2, ], c(-0.1497232, 4.1424036, -0.9816178), tolerance=1e-3)
  expect_equal(omega_n23t[3, ], c(4.0933073, 0.2691352, 33.7264024), tolerance=1e-3)
  expect_equal(omega_n23t[1, ], c(0.6997594, -0.3466295, 4.0933073), tolerance=1e-3)
  expect_equal(omega_n23gs[2, ], c(-0.4480529, 6.3617625, -0.5927894), tolerance=1e-3)
  expect_equal(omega_n13gsr[1, ], c(0.4162614, -0.3976433, 2.4052323), tolerance=1e-3)
})

g0_ac <- function(r, lag) {
  sapply((1 + lag):length(r), function(i1) sapply(1:lag, function(i2) r[i1]*r[i1 - i2]))
}
g_ac1 <- function(r) {
  lag <- 1
  if(lag > 1) {
    return(t(g0_ac(r, lag)))
  } else {
    return(as.matrix(g0_ac(r, lag)))
  }
}
g_ac3 <- function(r) {
  lag <- 3
  if(lag > 1) {
    return(t(g0_ac(r, lag)))
  } else {
    return(as.matrix(g0_ac(r, lag)))
  }
}

omega_ac1_11t <- get_test_Omega(simudata, 1, 1, params11t, model="StMAR", g=g_ac1, dim_g=1)
omega_ac1_12 <- get_test_Omega(simudata, 1, 2, params12, model="GMAR", g=g_ac1, dim_g=1)
omega_ac1_12r <- get_test_Omega(simudata, 1, 2, params12r, restricted=TRUE, g=g_ac1, dim_g=1)
omega_ac1_23t <- get_test_Omega(simudata, 2, 3, params23t, model="StMAR", g=g_ac1, dim_g=1)

omega_ac3_11t <- get_test_Omega(simudata, 1, 1, params11t, model="StMAR", g=g_ac3, dim_g=3)
omega_ac3_12 <- get_test_Omega(simudata, 1, 2, params12, g=g_ac3, dim_g=3)
omega_ac3_12r <- get_test_Omega(simudata, 1, 2, params12r, restricted=TRUE, g=g_ac3, dim_g=3)
omega_ac3_23t <- get_test_Omega(simudata, 2, 3, params23t, model="StMAR", g=g_ac3, dim_g=3)

omega_ac1_12gs <- get_test_Omega(simudata, 1, c(1,1), params12gs, model="G-StMAR", g=g_ac1, dim_g=1)
omega_ac1_13gsr <- get_test_Omega(simudata, 1, c(2,1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_ac1, dim_g=1)
omega_ac3_12gs <- get_test_Omega(simudata, 1, c(1,1), params12gs, model="G-StMAR", g=g_ac3, dim_g=3)
omega_ac3_13gsr <- get_test_Omega(simudata, 1, c(2,1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_ac3, dim_g=3)

test_that("get_test_Omega works for autocorrelation test", {
  expect_equal(omega_ac1_12gs[1, 1], 0.9821655, tolerance=1e-3)
  expect_equal(omega_ac1_13gsr[1, 1], 0.9637771, tolerance=1e-3)
  expect_equal(omega_ac3_12gs[2, ], c(-0.17748730, 0.90975252, 0.02570752), tolerance=1e-3)
  expect_equal(omega_ac3_13gsr[3, ], c(-0.32905727, 0.02084744, 0.97008598), tolerance=1e-3)
  expect_equal(omega_ac1_11t[1, 1], 1.12772, tolerance=1e-3)
  expect_equal(omega_ac1_12[1, 1], 0.9295376, tolerance=1e-3)
  expect_equal(omega_ac1_12r[1, 1], 1.020508, tolerance=1e-3)
  expect_equal(omega_ac1_23t[1, 1], 0.434482, tolerance=1e-3)
  expect_equal(omega_ac3_11t[, 3], c(-0.42201697, 0.05511734, 1.17249224), tolerance=1e-3)
  expect_equal(omega_ac3_12[3, ], c(-0.3139409, 0.0465679, 0.8805554), tolerance=1e-3)
  expect_equal(omega_ac3_12[1, ], c(0.9340736, -0.1687184, -0.3139409), tolerance=1e-3)
  expect_equal(omega_ac3_12r[, 2], c(-0.18435431, 0.86086297, 0.03502557), tolerance=1e-3)
  expect_equal(omega_ac3_23t[, 1], c(0.42946504, -0.18449283, 0.02649484), tolerance=1e-3)
  expect_equal(omega_ac3_23t[2, ], c(-0.1844928, 0.7946713, -0.2194248), tolerance=1e-3)
})

g0_ch <- function(r, lag) {
  sapply((1 + lag):length(r), function(i1) sapply(1:lag, function(i2) (r[i1]^2 - 1)*r[i1 - i2]^2))
}
g_ch1 <- function(r) {
  lag <- 1
  if(lag > 1) {
    return(t(g0_ch(r, lag)))
  } else {
    return(as.matrix(g0_ch(r, lag)))
  }
}
g_ch3 <- function(r) {
  lag <- 3
  if(lag > 1) {
    return(t(g0_ch(r, lag)))
  } else {
    return(as.matrix(g0_ch(r, lag)))
  }
}

omega_ch1_11t <- get_test_Omega(simudata, 1, 1, params11t, model="StMAR", g=g_ch1, dim_g=1)
omega_ch1_12 <- get_test_Omega(simudata, 1, 2, params12, g=g_ch1, dim_g=1)
omega_ch1_12r <- get_test_Omega(simudata, 1, 2, params12r, restricted=TRUE, g=g_ch1, dim_g=1)
omega_ch1_23t <- get_test_Omega(simudata, 2, 3, params23t, model="StMAR", g=g_ch1, dim_g=1)

omega_ch3_11t <- get_test_Omega(simudata, 1, 1, params11t, model="StMAR", g=g_ch3, dim_g=3)
omega_ch3_12 <- get_test_Omega(simudata, 1, 2, params12, g=g_ch3, dim_g=3)
omega_ch3_12r <- get_test_Omega(simudata, 1, 2, params12r, restricted=TRUE, g=g_ch3, dim_g=3)
omega_ch3_23t <- get_test_Omega(simudata, 2, 3, params23t, model="StMAR", g=g_ch3, dim_g=3)

omega_ch1_12gs <- get_test_Omega(simudata, 1, c(1,1), params12gs, model="G-StMAR", g=g_ch1, dim_g=1)
omega_ch1_13gsr <- get_test_Omega(simudata, 1, c(2,1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_ch1, dim_g=1)
omega_ch3_12gs <- get_test_Omega(simudata, 1, c(1,1), params12gs, model="G-StMAR", g=g_ch3, dim_g=3)
omega_ch3_13gsr <- get_test_Omega(simudata, 1, c(2,1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_ch3, dim_g=3)

test_that("get_test_Omega works for conditional heteroskedasticity tests", {
  expect_equal(omega_ch1_12gs[1, 1], 4.823228, tolerance=1e-3)
  expect_equal(omega_ch1_13gsr[1, 1], 5.350241, tolerance=1e-3)
  expect_equal(omega_ch3_13gsr[3, ], c(1.5098061, -0.2062663, 4.1129781), tolerance=1e-3)
  expect_equal(omega_ch3_12gs[, 1], c(4.861730, -0.730305, 1.025121), tolerance=1e-3)
  expect_equal(omega_ch1_11t[1, 1], 4.234744, tolerance=1e-3)
  expect_equal(omega_ch1_12[1, 1], 4.960322, tolerance=1e-3)
  expect_equal(omega_ch1_12r[1, 1], 5.983109, tolerance=1e-3)
  expect_equal(omega_ch1_23t[1, 1], 3.638563, tolerance=1e-3)
  expect_equal(omega_ch3_11t[, 1], c(4.2835988, -0.1581301, 1.2795732), tolerance=1e-3)
  expect_equal(omega_ch3_12[2, ], c(-0.6596341, 2.3962239, -0.1877080), tolerance=1e-3)
  expect_equal(omega_ch3_12[, 3], c(1.219144, -0.187708, 4.560058), tolerance=1e-3)
  expect_equal(omega_ch3_12r[1, 2], -0.9075125, tolerance=1e-3)
  expect_equal(omega_ch3_23t[3, ], c(0.15414291, -0.05188368, 6.49308987), tolerance=1e-3)
  expect_equal(omega_ch3_23t[, 2], c(-0.95669286, 3.53875813, -0.05188368), tolerance=1e-3)
})
