library(uGMAR)
context("getOmega")

## The results slightly differ depending on whether numerical integration (without the package "gsl")
# or hypergeometric function is used to calculate the quantile residuals. Also, the limited precission
# of floating point numbers (default used) might accumulate to differences big enough for the tests to fail
# (this particularly happens if the parameters contain overly large degree of freedom parameters). As the
# tests are moreover computationally demanding, they are commented out and used for development only.
test_that("quantileResiduals works", {
  expect_equal(quantileResiduals_int(VIX, 1, 1, c(-2, 0.8, 1, 12), model="StMAR")[13], 0.9494652, tolerance=1e-3)
})

# # getOmega calls other functions to deal with constraints, so they are not tested separately
# params11t <- c(-2, 0.8, 1, 12)
# params12 <- c(1.1, 0.9, 0.29, 4.5, 0.7, 3.2, 0.8)
# params23t <- c(1.8, 0.9, -0.06, 0.4, 7.2, 0.3, -0.009, 0.1, 3.1, 0.83, -0.05, 3.7, 0.7, 0.2, 11, 339, 198)
# params12r <- c(1, 4, 0.8, 0.3, 3, 0.8)
# params12gs <- c(1.2, 0.8, 0.6, 1.3, 0.6, 1.1, 0.6, 3) # M1=1, M2=1
# params13gsr <- c(1.3, 2.2, 1.4, 0.8, 2.4, 4.6, 0.4, 0.25, 0.15, 20) # M1=2, M2=1
#
# g_norm <- function(r) cbind(r^2-1, r^3, r^4-3)
# omega_n11t <- getOmega(VIX, 1, 1, params11t, model="StMAR", g=g_norm, dim_g=3)
# omega_n12 <- getOmega(VIX, 1, 2, params12, g=g_norm, dim_g=3)
# omega_n12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_norm, dim_g=3)
# omega_n23t <- getOmega(VIX, 2, 3, params23t, model="StMAR", g=g_norm, dim_g=3)
# omega_n23gs <- getOmega(VIX, 1, c(1, 1), params12gs, model="G-StMAR", g=g_norm, dim_g=3)
# omega_n13gsr <- getOmega(VIX, 1, c(2, 1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_norm, dim_g=3)
#
# test_that("getOmega works for normality tests", {
#   expect_equal(omega_n11t[2, ], c(29.36079, 52.05743, 87.65027), tolerance=1e-4)
#   expect_equal(omega_n12[1, 1], 0.409619, tolerance=1e-4)
#   expect_equal(omega_n12[2, 1], 0.1073828, tolerance=1e-4)
#   expect_equal(omega_n12r[, 3], c(32.55617, -50.83875, 89.97590), tolerance=1e-4)
#   expect_equal(omega_n12r[2, 3], -50.83875, tolerance=1e-4)
#   expect_equal(omega_n23t[3, ], c(9.476058, 31.364721, 145.473749), tolerance=1e-4)
#   expect_equal(omega_n23t[3, 2], 31.36472, tolerance=1e-4)
#   expect_equal(omega_n23gs[2, 3], 94.66165, tolerance=1e-4)
#   expect_equal(omega_n13gsr[1, ], c(1.5743312, 0.6217055, 2.8404580), tolerance=1e-4)
# })
#
# g0_ac <- function(r, lag) {
#   sapply((1+lag):length(r), function(i1) sapply(1:lag, function(i2) r[i1]*r[i1-i2]))
# }
# g_ac1 <- function(r) {
#   lag <- 1
#   if(lag>1) {
#     return(t(g0_ac(r, lag)))
#   } else {
#     return(as.matrix(g0_ac(r, lag)))
#   }
# }
# g_ac3 <- function(r) {
#   lag <- 3
#   if(lag>1) {
#     return(t(g0_ac(r, lag)))
#   } else {
#     return(as.matrix(g0_ac(r, lag)))
#   }
# }
#
# omega_ac1_11t <- getOmega(VIX, 1, 1, params11t, model="StMAR", g=g_ac1, dim_g=1)
# omega_ac1_12 <- getOmega(VIX, 1, 2, params12, model="GMAR", g=g_ac1, dim_g=1)
# omega_ac1_12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_ac1, dim_g=1)
# omega_ac1_23t <- getOmega(VIX, 2, 3, params23t, model="StMAR", g=g_ac1, dim_g=1)
#
# omega_ac3_11t <- getOmega(VIX, 1, 1, params11t, model="StMAR", g=g_ac3, dim_g=3)
# omega_ac3_12 <- getOmega(VIX, 1, 2, params12, g=g_ac3, dim_g=3)
# omega_ac3_12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_ac3, dim_g=3)
# omega_ac3_23t <- getOmega(VIX, 2, 3, params23t, model="StMAR", g=g_ac3, dim_g=3)
#
# omega_ac1_12gs <- getOmega(VIX, 1, c(1,1), params12gs, model="G-StMAR", g=g_ac1, dim_g=1)
# omega_ac1_13gsr <- getOmega(VIX, 1, c(2,1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_ac1, dim_g=1)
# omega_ac3_12gs <- getOmega(VIX, 1, c(1,1), params12gs, model="G-StMAR", g=g_ac3, dim_g=3)
# omega_ac3_13gsr <- getOmega(VIX, 1, c(2,1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_ac3, dim_g=3)
#
# test_that("getOmega works for autocorrelation test", {
#   expect_equal(omega_ac1_12gs[1, 1], 31.9099, tolerance=1e-4)
#   expect_equal(omega_ac1_13gsr[1, 1], 0.8864289, tolerance=1e-4)
#   expect_equal(omega_ac3_12gs[2, ], c(31.72712, 32.00536, 32.01475), tolerance=1e-4)
#   expect_equal(omega_ac3_13gsr[3, ], c(0.8803784, 0.8692675, 0.9629880), tolerance=1e-4)
#   expect_equal(omega_ac1_11t[1, 1], 16.76182, tolerance=1e-4)
#   expect_equal(omega_ac1_12[1, 1], 0.533819, tolerance=1e-4)
#   expect_equal(omega_ac1_12r[1, 1], 15.07792, tolerance=1e-4)
#   expect_equal(omega_ac1_23t[1, 1], 0.8719067, tolerance=1e-4)
#   expect_equal(omega_ac3_11t[, 3], c(16.68747, 16.67081, 16.77604), tolerance=1e-4)
#   expect_equal(omega_ac3_12[3, 2], -0.008611741, tolerance=1e-4)
#   expect_equal(omega_ac3_12[1, 3], -0.2138224, tolerance=1e-4)
#   expect_equal(omega_ac3_12r[, 2], c(14.78865, 15.01056, 14.53314), tolerance=1e-4)
#   expect_equal(omega_ac3_23t[1, 1], 0.8769805, tolerance=1e-4)
#   expect_equal(omega_ac3_23t[2, ], c(0.1645746, 0.7979296, 0.0148655), tolerance=1e-4)
# })
#
# g0_ch <- function(r, lag) {
#   sapply((1+lag):length(r), function(i1) sapply(1:lag, function(i2) (r[i1]^2-1)*r[i1-i2]^2))
# }
# g_ch1 <- function(r) {
#   lag <- 1
#   if(lag>1) {
#     return(t(g0_ch(r, lag)))
#   } else {
#     return(as.matrix(g0_ch(r, lag)))
#   }
# }
# g_ch3 <- function(r) {
#   lag <- 3
#   if(lag>1) {
#     return(t(g0_ch(r, lag)))
#   } else {
#     return(as.matrix(g0_ch(r, lag)))
#   }
# }
#
# omega_ch1_11t <- getOmega(VIX, 1, 1, params11t, model="StMAR", g=g_ch1, dim_g=1)
# omega_ch1_12 <- getOmega(VIX, 1, 2, params12, g=g_ch1, dim_g=1)
# omega_ch1_12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_ch1, dim_g=1)
# omega_ch1_23t <- getOmega(VIX, 2, 3, params23t, model="StMAR", g=g_ch1, dim_g=1)
#
# omega_ch3_11t <- getOmega(VIX, 1, 1, params11t, model="StMAR", g=g_ch3, dim_g=3)
# omega_ch3_12 <- getOmega(VIX, 1, 2, params12, g=g_ch3, dim_g=3)
# omega_ch3_12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_ch3, dim_g=3)
# omega_ch3_23t <- getOmega(VIX, 2, 3, params23t, model="StMAR", g=g_ch3, dim_g=3)
#
# omega_ch1_12gs <- getOmega(VIX, 1, c(1,1), params12gs, model="G-StMAR", g=g_ch1, dim_g=1)
# omega_ch1_13gsr <- getOmega(VIX, 1, c(2,1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_ch1, dim_g=1)
# omega_ch3_12gs <- getOmega(VIX, 1, c(1,1), params12gs, model="G-StMAR", g=g_ch3, dim_g=3)
# omega_ch3_13gsr <- getOmega(VIX, 1, c(2,1), params13gsr, model="G-StMAR", restricted=TRUE, g=g_ch3, dim_g=3)
#
# test_that("getOmega works for conditional heteroskedasticity tests", {
#   expect_equal(omega_ch1_12gs[1, 1], 16.83635, tolerance=1e-4)
#   expect_equal(omega_ch1_13gsr[1, 1], 1.034198, tolerance=1e-4)
#   expect_equal(omega_ch3_13gsr[3, ], c(0.7978984, 0.5466978, 1.0220235), tolerance=1e-4)
#   expect_equal(omega_ch3_12gs[, 1], c(17.09112, 17.61228, 17.30636), tolerance=1e-4)
#   expect_equal(omega_ch1_11t[1, 1], 43.6595, tolerance=1e-4)
#   expect_equal(omega_ch1_12[1, 1], 2.152529, tolerance=1e-4)
#   expect_equal(omega_ch1_12r[1, 1], 33.58398, tolerance=1e-4)
#   expect_equal(omega_ch1_23t[1, 1], 10.93848, tolerance=1e-4)
#   expect_equal(omega_ch3_11t[1, 1], 43.99142, tolerance=1e-4)
#   expect_equal(omega_ch3_12[2, ], c(0.1145639, 3.1456096, 0.2640641), tolerance=1e-4)
#   expect_equal(omega_ch3_12[, 3], c(0.4926400, 0.2640641, 3.1821835), tolerance=1e-4)
#   expect_equal(omega_ch3_12r[1, 2], 33.95629, tolerance=1e-4)
#   expect_equal(omega_ch3_23t[3, ], c(-1.54953644, -0.05031768, 4.42308303), tolerance=1e-4)
#   expect_equal(omega_ch3_23t[, 2], c(1.93888644, 5.14543685, -0.05031768), tolerance=1e-4)
# })
