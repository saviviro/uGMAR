library(uGMAR)
context("getOmega")

# Omegas are very sensitive to even very small changes in the quantile residuals. The quantile residuals are
# calculated only with numerical integration if package "gsl" is not available, which will cause these very
# small differences (in 1e-6 or 1e-8 tolerance for eample) that will then accumulate to big differences in Omega.
# This is why these tests are commented out (and used only for development).
test_that("quantileResiduals works", {
  expect_equal(quantileResiduals_int(VIX, 1, 1, c(-2, 0.8, 1, 12), StMAR=TRUE)[13], 0.9494652, tolerance=1e-3)
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
# omega_n11t <- getOmega(VIX, 1, 1, params11t, StMAR=TRUE, g=g_norm, dim_g=3)
# omega_n12 <- getOmega(VIX, 1, 2, params12, g=g_norm, dim_g=3)
# omega_n12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_norm, dim_g=3)
# omega_n23t <- getOmega(VIX, 2, 3, params23t, StMAR=TRUE, g=g_norm, dim_g=3)
# omega_n23gs <- getOmega(VIX, 1, c(1, 1), params12gs, GStMAR=TRUE, g=g_norm, dim_g=3)
# omega_n13gsr <- getOmega(VIX, 1, c(2, 1), params13gsr, GStMAR=TRUE, restricted=TRUE, g=g_norm, dim_g=3)
#
# test_that("getOmega works for normality tests", {
#   expect_equal(omega_n11t[2, 2], 52.05756, tolerance=1e-3)
#   expect_equal(omega_n12[1,1], 0.4096198, tolerance=1e-3)
#   expect_equal(omega_n12[2,1], 0.1073858, tolerance=1e-3)
#   expect_equal(omega_n12r[1,3], 32.55611, tolerance=1e-3)
#   expect_equal(omega_n12r[2,3], -50.83869, tolerance=1e-3)
#   expect_equal(omega_n23t[3,3], 272641.4, tolerance=1e-1)
#   expect_equal(omega_n23t[3,2], 54713.53, tolerance=1e-1)
#   expect_equal(omega_n23gs[2, 3], 94.66124, tolerance=1e-3)
#   expect_equal(omega_n13gsr[2, 2], 1.192931, tolerance=1e-3)
# })
#
# g0_ac <- function(r, lag) {
#   sapply((1+lag):length(r), function(i1) sapply(1:lag, function(i2) r[i1]*r[i1-i2]))
# }
# g_ac1 <- function(r) {
#   lag = 1
#   if(lag>1) {
#     return(t(g0_ac(r, lag)))
#   } else {
#     return(as.matrix(g0_ac(r, lag)))
#   }
# }
# g_ac3 <- function(r) {
#   lag = 3
#   if(lag>1) {
#     return(t(g0_ac(r, lag)))
#   } else {
#     return(as.matrix(g0_ac(r, lag)))
#   }
# }
#
# omega_ac1_11t <- getOmega(VIX, 1, 1, params11t, StMAR=TRUE, g=g_ac1, dim_g=1)
# omega_ac1_12 <- getOmega(VIX, 1, 2, params12, g=g_ac1, dim_g=1)
# omega_ac1_12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_ac1, dim_g=1)
# omega_ac1_23t <- getOmega(VIX, 2, 3, params23t, StMAR=TRUE, g=g_ac1, dim_g=1)
#
# omega_ac3_11t <- getOmega(VIX, 1, 1, params11t, StMAR=TRUE, g=g_ac3, dim_g=3)
# omega_ac3_12 <- getOmega(VIX, 1, 2, params12, g=g_ac3, dim_g=3)
# omega_ac3_12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_ac3, dim_g=3)
# omega_ac3_23t <- getOmega(VIX, 2, 3, params23t, StMAR=TRUE, g=g_ac3, dim_g=3)
#
# omega_ac1_12gs <- getOmega(VIX, 1, c(1,1), params12gs, GStMAR=TRUE, g=g_ac1, dim_g=1)
# omega_ac1_13gsr <- getOmega(VIX, 1, c(2,1), params13gsr, GStMAR=TRUE, restricted=TRUE, g=g_ac1, dim_g=1)
# omega_ac3_12gs <- getOmega(VIX, 1, c(1,1), params12gs, GStMAR=TRUE, g=g_ac3, dim_g=3)
# omega_ac3_13gsr <- getOmega(VIX, 1, c(2,1), params13gsr, GStMAR=TRUE, restricted=TRUE, g=g_ac3, dim_g=3)
#
# test_that("getOmega works for autocorrelation test", {
#   expect_equal(omega_ac1_12gs[1,1], 31.90999, tolerance=1e-3)
#   expect_equal(omega_ac1_13gsr[1,1], 0.8864362, tolerance=1e-3)
#   expect_equal(omega_ac3_12gs[2,3], 31.99145, tolerance=1e-3)
#   expect_equal(omega_ac3_13gsr[3,1], 0.889686, tolerance=1e-3)
#   expect_equal(omega_ac1_11t[1,1], 16.75613, tolerance=1e-3)
#   expect_equal(omega_ac1_12[1,1], 0.5338139, tolerance=1e-3)
#   expect_equal(omega_ac1_12r[1,1], 15.08044, tolerance=1e-3)
#   expect_equal(omega_ac1_23t[1,1], 1.378263, tolerance=1e-3)
#   expect_equal(omega_ac3_11t[3,3], 16.79347, tolerance=1e-3)
#   expect_equal(omega_ac3_12[3,2], -0.01035503, tolerance=1e-3)
#   expect_equal(omega_ac3_12[1,3], -0.221124, tolerance=1e-3)
#   expect_equal(omega_ac3_12r[2,2], 17.24854, tolerance=1e-3)
#   expect_equal(omega_ac3_23t[1,1], 1.393981, tolerance=1e-3)
#   expect_equal(omega_ac3_23t[2,1], -1.86544, tolerance=1e-3)
# })
#
# g0_ch <- function(r, lag) {
#   sapply((1+lag):length(r), function(i1) sapply(1:lag, function(i2) (r[i1]^2-1)*r[i1-i2]^2))
# }
# g_ch1 <- function(r) {
#   lag = 1
#   if(lag>1) {
#     return(t(g0_ch(r, lag)))
#   } else {
#     return(as.matrix(g0_ch(r, lag)))
#   }
# }
# g_ch3 <- function(r) {
#   lag = 3
#   if(lag>1) {
#     return(t(g0_ch(r, lag)))
#   } else {
#     return(as.matrix(g0_ch(r, lag)))
#   }
# }
#
# omega_ch1_11t <- getOmega(VIX, 1, 1, params11t, StMAR=TRUE, g=g_ch1, dim_g=1)
# omega_ch1_12 <- getOmega(VIX, 1, 2, params12, g=g_ch1, dim_g=1)
# omega_ch1_12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_ch1, dim_g=1)
# omega_ch1_23t <- getOmega(VIX, 2, 3, params23t, StMAR=TRUE, g=g_ch1, dim_g=1)
#
# omega_ch3_11t <- getOmega(VIX, 1, 1, params11t, StMAR=TRUE, g=g_ch3, dim_g=3)
# omega_ch3_12 <- getOmega(VIX, 1, 2, params12, g=g_ch3, dim_g=3)
# omega_ch3_12r <- getOmega(VIX, 1, 2, params12r, restricted=TRUE, g=g_ch3, dim_g=3)
# omega_ch3_23t <- getOmega(VIX, 2, 3, params23t, StMAR=TRUE, g=g_ch3, dim_g=3)
#
# omega_ch1_12gs <- getOmega(VIX, 1, c(1,1), params12gs, GStMAR=TRUE, g=g_ch1, dim_g=1)
# omega_ch1_13gsr <- getOmega(VIX, 1, c(2,1), params13gsr, GStMAR=TRUE, restricted=TRUE, g=g_ch1, dim_g=1)
# omega_ch3_12gs <- getOmega(VIX, 1, c(1,1), params12gs, GStMAR=TRUE, g=g_ch3, dim_g=3)
# omega_ch3_13gsr <- getOmega(VIX, 1, c(2,1), params13gsr, GStMAR=TRUE, restricted=TRUE, g=g_ch3, dim_g=3)
#
# test_that("getOmega works for conditional heteroskedasticity tests", {
#   expect_equal(omega_ch1_12gs[1,1], 16.83659, tolerance=1e-3)
#   expect_equal(omega_ch1_13gsr[1,1], 1.034139, tolerance=1e-3)
#   expect_equal(omega_ch3_13gsr[3,3], 1.043636, tolerance=1e-3)
#   expect_equal(omega_ch3_12gs[2,1], 17.52431, tolerance=1e-3)
#   expect_equal(omega_ch1_11t[1,1], 43.66495, tolerance=1e-3)
#   expect_equal(omega_ch1_12[1,1], 2.149916, tolerance=1e-3)
#   expect_equal(omega_ch1_12r[1,1], 33.61512, tolerance=1e-3)
#   expect_equal(omega_ch1_23t[1,1], 26.12893, tolerance=1e-3)
#   expect_equal(omega_ch3_11t[1,1], 44.10323, tolerance=1e-3)
#   expect_equal(omega_ch3_12[2,2], 3.105598, tolerance=1e-3)
#   expect_equal(omega_ch3_12[1,3], 0.4683494, tolerance=1e-3)
#   expect_equal(omega_ch3_12r[1,2], 39.1227, tolerance=1e-3)
#   expect_equal(omega_ch3_23t[3,1], 35.43488, tolerance=1e-3)
#   expect_equal(omega_ch3_23t[2,3], -149.1611, tolerance=1e-3)
# })
