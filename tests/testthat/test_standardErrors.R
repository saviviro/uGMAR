library(uGMAR)
context("standard errors")

params11t <- c(0.162640526, 0.933158192, 0.005438104, 3.032187928)
params12 <- c(0.18280320, 0.92657724, 0.00214549, 0.85737192, 0.68206229, 0.01900195, 0.88340618)
params12gs <- c(0.273978569, 0.889069634, 0.002023354, 0.617278126, 0.771494763, 0.017264129, 0.800791275, 3.498366389)

test_that("standardErrors work", {
  expect_equal(standardErrors(logVIX, p=1, M=1, params=params11t, model="StMAR"),
               c(0.063922594, 0.025761040, 0.002954027, 0.885619265), tolerance=1e-2)
  expect_equal(standardErrors(logVIX, p=1, M=2, params=params12, model="GMAR"),
               c(0.069584459, 0.027736375, 0.000252121, 0.390350711, 0.142669772, 0.005994098, 0.058184980), tolerance=1e-2)
  expect_equal(standardErrors(logVIX, p=1, M=c(1, 1), params=params12gs, model="G-StMAR"),
               c(0.0994685509, 0.0404645265, 0.0002611683, 0.3657080324, 0.1336568285, 0.0153480868, 0.0983687861, 2.9642334081), tolerance=1e-2)
})


