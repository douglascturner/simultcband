test_that("output structure", {
  expect_is(cband(yobs~s(x), alpha=.05, predgen = pgen_gam, data=data.frame(x = 1:100, yobs = (1:100)^2)), "data.frame")
})

test_that("B size for Bonferroni test", {
  expect_error(cband(yobs~x, alpha=.1, predgen = pgen_loess, data=data.frame(x = 1:100, yobs = (1:100)^2), method = "Sidak",B=500),"Number of bootstrap replicates is too low for accurate estimation of bands at this alpha level.")
})

test_that('weight error', {
  expect_error(cband(yobs~WEIGHT_FOR_BOOT, predgen = pgen_loess, data=data.frame(WEIGHT_FOR_BOOT = 1:100, yobs = (1:100)^2)), "Variable name of 'WEIGHT_FOR_BOOT' not permitted since used internally.")
  expect_error(cband(yobs~weight, predgen = pgen_loess, data=data.frame(weight = 1:100, yobs = (1:100)^2)), "Variable name of 'weight' not permitted since the 'boot' package would overwrite it.")
})

