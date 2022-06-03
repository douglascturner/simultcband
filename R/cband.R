# This file has the core functions that use pgens to estimate simultaneous
# confidence bands.


# this value is hardcoded in pgens.R so if change here, change there also.
weight_varname <- "WEIGHT_FOR_BOOT"
utils::globalVariables(weight_varname)


get_prediction <- function(predictor, data, newdata, se.fit, use_boot_weights = FALSE) {
  # This wrapper around the predictor is to centralize cleaning and sanity
  # checking of the result of the predictor, which can vary because the predictor
  # can be written by the user, so we try to be flexible in what we accept.
  # The cast2vec() call in this function converts from various formats to a
  # standard format that we return.

  if (use_boot_weights) {
    if (use_boot_weights && !weight_varname %in% names(data)) {
      # boot::boot() calls the the predictor once without setting weights for the
      # predictions from the original data. Instead of setting flat weights, we
      # could just call the predictor with "use_boot_weights = FALSE", but
      # the user might not expect the non-weighted variant to be used when
      # using a weighted bootstrap.
      data[, weight_varname] <- rep(1, nrow(data))
    }
    pred_ <- predictor(data = data, newdata = newdata, se.fit = se.fit, use_boot_weights = TRUE)
  } else {
    pred_ <- predictor(data = data, newdata = newdata, se.fit = se.fit)
  }

  if (se.fit) {
    if(!is.list(pred_) || !all(c("fit", "se.fit") %in% names(pred_))) {
      stop("The predictor must return a list with elements ",
           "named 'fit' and 'se.fit'.")
    }
  }

  if (se.fit) {
    pred_clean <- lapply(X = pred_, FUN = cast2vec)
  } else {
    pred_clean <- cast2vec(pred_)
  }
  return(pred_clean)
}


#' @importFrom stats rexp
data.rg.wbo <- function(data, mle) {
  bdata <- data
  n<-dim(data)[1]
  multipliers  <- rexp(n)
  bdata[, weight_varname] <- multipliers / sum(multipliers)
  return(bdata)
}


boot.stat.weight <- function(data, predictor, fmla, newdata, yhat_orig, robust_var) {
  pred1_ <- get_prediction(predictor = predictor, data = data,
                          newdata = newdata, se.fit = !robust_var, use_boot_weights = TRUE)
  if (robust_var) {
    return(pred1_)
  } else {
    cb <- max(abs(pred1_$fit - yhat_orig)/pred1_$se.fit)
    return(cb)
  }
}

boot.stat = function(x, b, predictor, newdata, yhat_orig, robust_var) {
  pred1_ <- get_prediction(predictor = predictor, data = x[b, ],
                         newdata = newdata, se.fit = !robust_var)
  if (robust_var) {
    # If robust variance, we will estimate the variance after the bootstraps
    # are done, so here we just return the predicted values. That is, we don't
    # actually calculate the real test statistic for this case. Here the
    # "statistic" is the vector of predictions.
    return(pred1_)
  } else {
    if (anyNA(pred1_$se.fit)) {
      stop("predictor must provide standard errors when 'robust_var' is FALSE")
    }
    cb <- max(abs(pred1_$fit - yhat_orig)/pred1_$se.fit)
    return(cb)
  }
}


boot.stat.fitted = function(x, b, predictor, newdata) {
  pred1_ <- get_prediction(predictor = predictor, data = x[b, ],
                         newdata = newdata, se.fit = FALSE)
  return(pred1_)
}


#' @importFrom boot boot
#' @importFrom stats quantile qnorm
emp_boot <- function(formula_, predictor, newdata, alpha = 0.05, data, B,
                     robust_var, weighted_boot, est_name = est_name) {
  # returns a critical value.

  # Construct boot vector and (possibily) matrix of fitted values
  pred_ <- get_prediction(predictor = predictor, data = data,
                        newdata = newdata, se.fit = TRUE)
  if (weighted_boot) {
    boot_m = boot(data = data, statistic = boot.stat.weight, sim = "parametric",
                  ran.gen = data.rg.wbo,
                  # MLE=NULL as no parameters passed to ran.gen
                  mle = NULL,
                  predictor = predictor, fmla = formula_, newdata = newdata,
                  yhat_orig = pred_$fit, robust_var = robust_var, R = B)
    data = data[, names(data) != weight_varname]
    boot_m = boot_m$t
  } else {
    boot_m <- boot(data, boot.stat, B, predictor = predictor, newdata = newdata,
                   yhat_orig = pred_$fit, robust_var = robust_var)$t
  }
  # boot_m is a matrix where the first row is the d_b vector if
  # robust_var is FALSE. If robust_var is TRUE, then it is a N x B matrix.
  # Rows are points in the grid and columns are bootstrap replicates.

  #Compute se
  if (robust_var) {
    fitted_values_m <- boot_m
    centered_fitval_m <- fitted_values_m - matrix(pred_$fit,
                                                  nrow = nrow(fitted_values_m),
                                                  ncol = ncol(fitted_values_m),
                                                  byrow = TRUE)

    # TODO: add this paper to references in the R help file.
    # robust variance estimator
    # see Remark 3.2 in https://arxiv.org/pdf/0904.0951.pdf
    # by default, use interquartile spread, but as the article mentions:
    #   "Other choices of quantile spreads are also possible."
    q1 <- 0.75
    quants <- apply(centered_fitval_m, 2, quantile, c(q1, 1 - q1), na.rm = TRUE)
    se_numerator <- quants[1, ] - quants[2, ]
    se_denominator <- qnorm(q1) - qnorm(1 - q1)
    se <- se_numerator/se_denominator

    # "bvec" = "boot vector" (i.e., vector of values from bootstraps)
    d_bvec <- apply(abs(centered_fitval_m /matrix(se, nrow = nrow(fitted_values_m), ncol = ncol(fitted_values_m), byrow = TRUE)), 1, max)
    # critical value
    crit <- quantile(d_bvec, 1 - alpha)
  } else {
    # Compute critical values.
    crit <- quantile(x = boot_m, probs = 1 - alpha)
    se <- pred_$se.fit
  }
  return(list(crit = crit, se = se))
}


get_fitted_val <-function(formula_, predictor,
                              newdata, yhat_orig, data) {
  boot_ind <- sample(nrow(data), replace = TRUE)
  data <- data[boot_ind, ]
  pred1_ <- get_prediction(predictor = predictor, data = data,
                         newdata = newdata, se.fit = FALSE)
  return(pred1_)
}


#' @importFrom boot boot
#' @importFrom stats quantile
alpha_correction <- function(formula_, predictor,
                             newdata, alpha, data, method, B) {
  N <- nrow(newdata)

  if (method == "Pointwise") {
    alpha_corrected = alpha
  } else if (method == "Bonferroni") {
    alpha_corrected = alpha/(N)
  } else if (method == "Sidak") {
    alpha_corrected = 1 - (1 - alpha)^(1/N)
  }

  if (B*alpha_corrected < 5) {
    # This threshold of 5 is arbitrary.
    stop("Number of bootstrap replicates is too low for accurate estimation of bands at this alpha level.")
  }

  fitted_values_m <- boot(data, boot.stat.fitted, B, predictor = predictor, newdata = newdata)$t
  fitted_values_m <- t(fitted_values_m) #rows are grid points, columns are bootstrap replicates
  bound_values_m <- apply(fitted_values_m, MARGIN = 1, FUN = function(point) {
    return(quantile(point, probs = c(alpha_corrected/(2), 1 - alpha_corrected/(2) ) ) )
  })
  return(t(bound_values_m))
}


#' Create a simultaneous confidence band
#'
#' @param formula The formula giving the specification.
#' @param alpha The significance level of the confidence band.
#' @param predgen The predictor generator. Can be the function or
#'        character name of the function. See ?predictor_generators for more
#'        information.
#' @param weighted_boot If TRUE, a multiplier bootstrap is used with normalized
#'        i.i.d. standard exponential multipliers. If FALSE, a non-parametric
#'        bootstrap is used. See Kosorok (2008) page 20.
#' @param newdata The grid at which to give predictions.
#' @param data The data.
#' @param method The confidence band method.
#' @param B Number of bootstrap replications.
#' @param robust_var if TRUE, use robust variance estimation.
#' @param ... dots are passed to predgen.
#'
#' @return A data frame with columns "lower", "est", and "upper".
#'
#' @examples
#' x <- rnorm(1000)
#' y <- rnorm(1000)
#' dframe <- data.frame(x = x, y = y)
#' # choose small 'B' to make quick; in practice B should be much larger.
#' cband_ex1 <- cband(formula = y ~ s(x), data = dframe, B = 10)
#' # for pgen_loc, we must use the verbose locfit::lp(x) and not just lp(x)
#' cband_ex2 <- cband(formula = y ~ locfit::lp(x), data = dframe,
#'                    predgen = pgen_loc, B = 10)
#'
#' @export
cband <- function(formula, alpha = 0.05, predgen = pgen_gam, newdata,
                  data, method = "Empirical bootstrap", B = 500,
                  robust_var = FALSE, weighted_boot = FALSE, ...) {
  if (missing(newdata)) {
    # optop (optimization opportunity):
    # could alternatively extract only the x variables from data
    # and pass those to newdata.
    # TODO: test whether this^ would speed things up with a very
    #       wide data frame (e.g., thousands of variables). It might
    #       make a difference for some pgens and not others.
    newdata <- data
  }
  if (weight_varname %in% colnames(data)) {
    stop("Variable name of '", weight_varname, "' not permitted since used internally.")
  }
  if ("weight" %in% colnames(data)) {
    stop("Variable name of 'weight' not permitted since the 'boot' package would overwrite it.")
  }
  if (is.character(predgen)) {
    est_name <- predgen
  } else {
    est_name <- toString(substitute(predgen))
  }

  if (!is.numeric(alpha))
    stop("The argument 'alpha' must be numeric.")
  if (!identical(length(alpha), 1L))
    stop("The argument 'alpha' must be of length 1.")
  if (alpha <= 0 || alpha >= 1)
    stop("The argument 'alpha' must be in (0, 1).")

  # this allows the arg to be either the function name (string) or the
  # function itself.
  predgen <- match.fun(predgen)

  predictor <- predgen(data_orig = data, formula = formula, ...)

  # optop: this is a duplicate calculation, since we also calculate inside
  # emp_boot(). We could alternatively pass pred_ as an arg to emp_boot(), or
  # have emp_boot() handle it and get it returned; or can have get_prediction
  # perform argument caching.
  pred_ <- get_prediction(predictor = predictor, data = data,
                        newdata = newdata, se.fit = TRUE)

  if (method == "Empirical bootstrap") {
    res <- emp_boot(formula_ = formula, predictor = predictor,
                               newdata = newdata, alpha = alpha, data = data,
                               B = B, robust_var = robust_var,
                               weighted_boot = weighted_boot, est_name = est_name)
    ret <- data.frame(lower = pred_$fit - res$crit*res$se,
                      est   = pred_$fit,
                      upper = pred_$fit + res$crit*res$se)
    return(ret)
  } else if (method %in% c("Bonferroni", "Sidak", "Pointwise")) {
    bound <- alpha_correction(formula_ = formula, predictor = predictor,
                newdata = newdata, alpha = alpha, data = data, method = method,
                B = B)
    ret <- data.frame(lower = bound[, 1],
                      est   = pred_$fit,
                      upper = bound[, 2])
    return(ret)
  } else {
    stop("unknown 'method'")
  }
}
