# These integration tests make it easy and quick to detect if output changes
# for various combinations of options and predictors. They rely heavily on
# set.seed() and thus anything that changes the seed will invalidate these
# tests and will require tests with a large number of bootstraps to ensure that
# there is no regression.
#
# These tests are similar in code to test-predictors.R. In the future we could
# consider convenience functions to centralize the shared code.

library(simultcband)

set.seed(1)

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
  # not all predictors support weighted_boot, so we only vary robust_var here. We test weighted_boot below.
  print(cband(formula = f1, alpha = 0.1, predgen = predictors_l[[i]], data = dat, B = 10, robust_var = FALSE, weighted_boot = FALSE))
  print(cband(formula = f1, alpha = 0.1, predgen = predictors_l[[i]], data = dat, B = 10, robust_var = TRUE, weighted_boot = FALSE))
}

for (i in seq_along(predictors_deriv_l)) {
  # pgen_spline does not estimate SE so need to set robust_var = TRUE.
  print(cband(formula = f1, alpha = 0.1, predgen = predictors_deriv_l[[i]], data = dat, B = 10, deriv_order = 1, robust_var = TRUE))
}

print(cband(formula = f_kde, alpha = 0.1, predgen = pgen_kde, data = dat_kde, B = 10))

print(cband(formula = f_lp, alpha = 0.1, predgen = pgen_loc, data = dat, B = 10))

print(cband(formula = f_s, alpha = 0.1, predgen = pgen_gam, data = dat, B = 10, weighted_boot = FALSE, robust_var = FALSE))
print(cband(formula = f_s, alpha = 0.1, predgen = pgen_gam, data = dat, B = 10, weighted_boot = FALSE, robust_var = TRUE))
print(cband(formula = f_s, alpha = 0.1, predgen = pgen_gam, data = dat, B = 10, weighted_boot = TRUE, robust_var = FALSE))
print(cband(formula = f_s, alpha = 0.1, predgen = pgen_gam, data = dat, B = 10, weighted_boot = TRUE, robust_var = TRUE))

print(cband(formula = f1, alpha = 0.1, predgen = pgen_lm, data = dat, B = 10, weighted_boot = FALSE, robust_var = FALSE))
print(cband(formula = f1, alpha = 0.1, predgen = pgen_lm, data = dat, B = 10, weighted_boot = FALSE, robust_var = TRUE))
print(cband(formula = f1, alpha = 0.1, predgen = pgen_lm, data = dat, B = 10, weighted_boot = TRUE, robust_var = FALSE))
print(cband(formula = f1, alpha = 0.1, predgen = pgen_lm, data = dat, B = 10, weighted_boot = TRUE, robust_var = TRUE))
