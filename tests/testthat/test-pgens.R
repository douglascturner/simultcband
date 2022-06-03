# These tests are just meant to check for mechanical errors, such as an
# unspecified dependency, or change in interface from an underlying package
# update.
#
# These tests are similar in code to output-tests. In the future we could
# consider convenience functions to centralize the shared code.


ssize <- 100
dat <- data.frame(x = rnorm(ssize), y = rnorm(ssize))
dat_kde <- data.frame(y = rnorm(ssize))

# todo: once we implement an abstract layer for a formula to
# work with multiple predictors, loop through predictors with
# just one formula to simplify things.

# formulas for the tests
f1 <- y ~ x
f_kde <- ~y
f_lp <- y ~ locfit::lp(x)
f_s <- y ~ s(x)

# todo: we should have a way to get a list of all predictors
# and all that support derivative prediction.
predictors_l <- c(
                  pgen_gam = pgen_gam,
                  pgen_lm = pgen_lm,
                  pgen_rq = pgen_rq,
                  pgen_loess = pgen_loess,
                  pgen_glm = pgen_glm,
                  pgen_npreg = pgen_npreg,
                  pgen_loc = pgen_loc
)

predictors_deriv_l <- c(
                  pgen_gam_deriv = pgen_gam_deriv,
                  pgen_spline = pgen_spline
)

for (i in seq_along(predictors_l)) {
  test_that(paste("Estimator", names(predictors_l)[[i]], "does not give error"), {
      expect_error(cband(formula = f1, alpha = 0.1, predgen = predictors_l[[i]], data = dat, B = 10),
                   NA)
  })
}

for (i in seq_along(predictors_deriv_l)) {
  test_that(paste("Deriv predictor", names(predictors_deriv_l)[[i]], "does not give error"), {
      # pgen_spline does not estimate SE so need to set robust_var = TRUE.
      expect_error(cband(formula = f1, alpha = 0.1, predgen = predictors_deriv_l[[i]], data = dat, B = 10, deriv_order = 1, robust_var = TRUE),
                   NA)
  })
}

test_that(paste("pgen_kde works with f_kde"), {
expect_error(cband(formula = f_kde, alpha = 0.1, predgen = pgen_kde, data = dat_kde, B = 10),
                   NA)
})

test_that(paste("pgen_loc works with f_lp"), {
expect_error(cband(formula = f_lp, alpha = 0.1, predgen = pgen_loc, data = dat, B = 10),
                   NA)
})

test_that(paste("pgen_gam works with f_s"), {
expect_error(cband(formula = f_s, alpha = 0.1, predgen = pgen_gam, data = dat, B = 10),
                   NA)
})

# bad because "sefit" should be "se.fit"
pgen_bad <- function(data_orig, formula, ...) {
  ret_fn <- function(data, newdata, se.fit, use_boot_weights = FALSE) {
    badret <- list(xyzfit = 3, xyzsefit = 5)
    return(badret)
  }
  return(ret_fn)
}
test_that("incorrectly formatted predictors give errors", {
expect_error(cband(formula = f1, alpha = 0.1, predgen = pgen_bad, data = dat, B = 10),
                   "The predictor must return a list with elements named 'fit' and 'se.fit'.")
})


## TODO {after public}: reference new issue after moving over old issue #24
##
## TODO: figure out whether there is something we should do to address this.
## if we use dat (that has an auxiliary column), it fails with the following error:
## <<length of bandwidth vector does not match number of columns of 'tdat'>>
#expect_error(cband(formula = f_kde, alpha = 0.1, predgen = pgen_kde, data = dat, B = 10),
#                   NA)
