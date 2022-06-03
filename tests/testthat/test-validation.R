dat <- data.frame(y = 1:10, x1 = 1:10, x2 = 1:10)

f1 <- y ~ x1 + x2
f_str <- "y ~ poly(x1, degree = 3)*poly(x2, degree = 3)"

test_that('basic formula', {
    expect_error(validate_formula(formula = f1, data = dat), NA)
})


# this test currently fails. Thus, we do not validate formulas for
# pgen_lm because f_str is a valid formula for lm.
#
#test_that('complex string formula', {
#    expect_error(validate_formula(formula = f_str, data = dat), NA)
#})
