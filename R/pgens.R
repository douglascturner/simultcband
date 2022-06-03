# This file has definitions for predictor generators (pgens) that we ship.

# See ?predictor_generators for user documentation. We document some
# additional internal details below.

# Decisions regarding the pgen design: this design would be inefficient if
# we wanted to call the predictor with the same data but with more than one
# newdata. We don't have a situation in mind where we would do this.

# To use many of the pgens effectively, it is necessary to use functions
# in the formula. For example, s() with pgen_gam, or lp() with
# pgen_loc(). In the case of pgen_gam(), the user can use s()
# directly and s() does not need to be imported in this file's definition of
# pgen_gam(). The case of pgen_loc() is more tricky, and currently we
# require the user to either use the more verbose locfit::lp() in the formula,
# or to load the package "locfit" themselves before using pgen_loc().

# Miscellaneous interface discussion:
#
# Ideally we would not need to explicitly accept "weights" arg, and would just
# pass it through '...', but we already use '...' from a higher level.
# We could alternatively merge both (i.e., manually change '...').

# Source code guide:
# optop: optimimization opportunity. There is potential for making the code
#        run faster. Often, these are micro-optimizations, i.e., not expected
#        to make a big difference for most use cases.


#' Estimation using mgcv::gam
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#' @param ... dots are passed to gam.
#'
#' @importFrom mgcv gam
#' @export
pgen_gam <- function(data_orig, formula, ...) {
  # the formula is required to be a gam formula

  # optop: this gam() call is unnecessary if length(sp_) == 0.
  #        We could instead check from "formula" arg.
  gam_ <- gam(formula = formula, data = data_orig)

  # We don't currently allow specifying one and not the other, e.g.:
  # gam(formula = y ~ s(lnl, sp = 1) + s(lnk), data = data)
  # This could be easily implemented, but is currently erroneously detected as
  # a bivariate smooth.
  sp_ <- gam_$sp

  # formula_fixed is what we will change to make sure sp is fixed
  formula_fixed <- formula
  if (length(sp_) == 0) {
    # We don't need to do anything:
    # formula_fixed will use the same as formula.
  } else if (length(sp_) == 1) {
    # This case is for, e.g., y ~ s(x1, x2)
    # we check that sure user did not already fix sp
    if (is.null(formula_fixed[[3]]$sp)) {
      formula_fixed[[3]]$sp <- sp_
    }
  } else if (length(sp_) == 2) {
    # This case is for, e.g., y ~ s(x1) + s(x2)
    formula_fixed[[3]][[2]]$sp <- sp_[[1]]
    formula_fixed[[3]][[3]]$sp <- sp_[[2]]
  } else {
    stop("We are not prepared for this formula type. ",
         "Please make a feature request.")
  }

  ret_fn <- function(data, newdata, se.fit, use_boot_weights = FALSE) {
    if (use_boot_weights) {
      gam_evaluated <- gam(formula = formula_fixed, data = data, weights = WEIGHT_FOR_BOOT, ...)
    } else {
      gam_evaluated <- gam(formula = formula_fixed, data = data, ...)
    }
    pred_ <- predict(gam_evaluated, newdata = newdata, se.fit = se.fit)
    return(pred_)
  }
  return(ret_fn)
}


#' Estimation using stats::lm
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#' @param ... dots are passed to lm.
#'
#' @export
#' @importFrom stats lm predict
pgen_lm <- function(data_orig, formula, ...) {
  # we do not validate because we incorrectly give errors in some cases. See the
  # commented out test here:
  #   test-validation.R
  #validate_formula(formula = formula, data = data_orig)

  ret_fn <- function(data, newdata, se.fit, use_boot_weights = FALSE) {
    if (use_boot_weights) {
      lm_ <- lm(formula = formula, data = data, weights = WEIGHT_FOR_BOOT, ...)
    } else {
      lm_ <- lm(formula = formula, data = data, ...)
    }
    pred_ <- predict(lm_, newdata = newdata, se.fit = se.fit)
    return(pred_)
  }
  return(ret_fn)
}


# We don't export this function to avoid confusion.
# It's only for testing that we get the same as relying on predict(), as in
# pgen_lm; and as a simple template for future pgens.
#' @importFrom stats lm coef as.formula model.matrix vcov
pgen_lm_nopredict <- function(data_orig, formula, ...) {
  # formula can sometimes be given as a character vector (of length 1).
  # For example:
  #   "y ~ poly(x1, degree = 3)*poly(x2, degree = 3)"
  # We now convert it to a formula because we use formula methods below
  # (e.g., we rely on formula[[3]] returning the RHS).
  formula <- as.formula(formula)

  ret_fn <- function(data, newdata, se.fit) {
    lm_ <- lm(formula = formula, data = data, ...)
    # rhs_f will be something like "~x1 + x2" we do not want the lhs term in
    # the formula, because if it is not in newdata (which is often the case if)
    # model.matrix will give an error.
    rhs_f <- as.formula(paste0("~ ", deparse(formula[[3]])))
    newdata_mm <- model.matrix(rhs_f, data = newdata)
    fit_ <- newdata_mm %*% coef(lm_)
    if (se.fit) {
      vcov_ <- vcov(lm_)
      # Uses an optimized way to calculate only the diagonal of
      # the var matrix. See equivalent (and more readable) below.
      se_ <- sqrt(rowSums((newdata_mm %*% vcov_) * newdata_mm))
      # TODO: create an assert option for extra checks.
      #       also do this assert in pgen_rq which uses
      #       the same optimization.
      # v_xb_new <- newdata_mm %*% vcov_ %*% t(newdata_mm)
      # se2_ <- sqrt(diag(v_xb_new))
      # stopifnot(isTRUE(all.equal(se_, se2_)))
      pred_ <- list(fit = fit_, se.fit = se_)
    } else {
      pred_ <- fit_
    }
    return(pred_)
  }
  return(ret_fn)
}


#' Estimation quantreg::rq
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#' @param ... dots are passed to rq.
#'
#' @importFrom quantreg rq
#' @importFrom stats coef as.formula model.matrix
#' @export
pgen_rq <- function(data_orig, formula, ...) {
  formula <- as.formula(formula)

  ret_fn <- function(data, newdata, se.fit) {
    rq_ <- rq(formula = formula, data = data, ...)
    rhs_f <- as.formula(paste0("~ ", deparse(formula[[3]])))
    newdata_mm <- model.matrix(rhs_f, data = newdata)
    fit_ <- newdata_mm %*% coef(rq_)
    if (se.fit) {
      # There are five possible ways to get SEs. First, we list
      # the options that give the fewest errors and warnings.
      #
      ## no warnings or errors
      # "ker"
      ## many warnings of "solution may be non-unique", but no errors
      # "iid"
      ## non-positives
      # "boot"
      ## non-positives
      # "nid"
      ## returns stats::cov.
      # "BLB"
      se_method <- "ker"

      # TODO: document this. Also give explicit examples of warnings.
      #       Maybe capture the warnings and use a list of "expected"
      #       warnings so we can still give error on unexpected warnings.
      #
      # for some methods, we need to allow warnings
      reset_warn_to_2 <- FALSE
      if (getOption("warn") == 2 && se_method == "iid") {
        reset_warn_to_2 <- TRUE
        options(warn = 0)
      }
      summary_ <- summary(rq_, covariance = TRUE, se = se_method)
      if (reset_warn_to_2) {
        options(warn = 2)
      }

      vcov_ <- summary_$cov
      # Uses an optimized way to calculate only the diagonal of
      # the var matrix.
      se_ <- sqrt(rowSums((newdata_mm %*% vcov_) * newdata_mm))
      pred_ <- list(fit = fit_, se.fit = se_)
    } else {
      pred_ <- fit_
    }
    return(pred_)
  }
  return(ret_fn)
}


#' Estimation using stats::loess
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#' @param ... dots are passed to loess.
#'
#' @export
#' @importFrom stats loess loess.control
pgen_loess <- function(data_orig, formula, ...) {
  validate_formula(formula = formula, data = data_orig)
  ret_fn <- function(data, newdata, se.fit) {
    loess_ <- loess(formula = formula, data = data, control = loess.control(surface = "direct"), ...)
    pred_ <- predict(loess_, newdata = newdata, se = se.fit)
    return(pred_)
  }
  return(ret_fn)
}


#' Estimation using stats::glm
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#' @param ... dots are passed to glm.
#'
#' @export
#' @importFrom stats predict glm
pgen_glm <- function(data_orig, formula, ...) {
  ret_fn <- function(data, newdata, se.fit) {
    glm_ <- glm(formula = formula, data = data, ...)
    pred_ <- predict(glm_, newdata = newdata, se.fit = se.fit, type = "response")
    return(pred_)
  }
  return(ret_fn)
}


#' Estimation using np::npudens
#'
#' @param data_orig The user's original data passed to the client function.
#'        The data cannot contain columns of variables not in the formula.
#' @param formula The user's original formula passed to the client function.
#'        Must be of the form '~x'.
#' @param ... dots are passed to npudens.
#'
#' @importFrom np npudensbw npudens se
#' @importFrom stats fitted
#' @export
pgen_kde <- function(data_orig, formula, ...) {
  bw <- npudensbw(formula = formula, data = data_orig)$bw
  var_names <- all.vars(formula[[-1]])
  ret_fn <- function(data, newdata, se.fit) {
    # Column labels are dropped in bootstrap indexing.
    if (!(is.data.frame(data))) {
      data <- data.frame(data)
      colnames(data) <- var_names
    }
    kde_ <- npudens(formula = formula, tdat = data, edat = newdata, bws = bw, ...)
    pred_ <-list(fit = fitted(kde_), se.fit = se(kde_))
    return(pred_)
  }
  return(ret_fn)
}


#' Estimation using np::npregbw
#'
#' Use npreg's local linear (regtype = "ll") estimation.
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#' @param regtype_ use local linear ("ll") or local constant, i.e., Nadaraya-Watson.
#' @param ... dots are passed to npregbw.
#'
#' @importFrom np npregbw npreg
#' @export
pgen_npreg <- function(data_orig, formula, regtype_ = "ll", ...) {
  #fix bandwidth based on data_orig
  bw_ <- npregbw(formula, data = data_orig, regtype = regtype_, ...)$bw
  f <- formula(formula)
  vars <- all.vars(f)
  ret_fn <- function(data, newdata, se.fit, bw = bw_) {
    #Make explanatory dataframe
    expldata <- data[, vars[2:length(vars)]]
    # Make response dataframe
    respdata <- data[, vars[1]]
    # Make new explanatory dataframe from newdata
    newexpldata <- newdata[, vars[2:length(vars)]]
    # Fit and Predict
    np_ <- npreg(bws = bw_, txdat = expldata, tydat = respdata, regtype = regtype_, ...)
    pred_ <- predict(np_, data.frame(newexpldata), se.fit = se.fit)

    return(pred_)
  }
  return(ret_fn)
}


#' Estimation of derivative using mgcv::gam
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#' @param ... dots are passed to gam.
#' @param deriv_order order of derivative.
#'
#' @importFrom mgcv gam
#' @importFrom stats coef terms.formula
#' @export
pgen_gam_deriv <- function(data_orig, formula, deriv_order, ...) {
  # see also comments in pgen_gam().

  # the formula is required to be a gam formula
  nterms <- ncol(attr(terms.formula(formula), "factors"))
  if (nterms != 1) {
    stop("Only 1-d bands supported for this pgen.")
  }
  if (deriv_order > 2) {
    stop("This order of derivative not supported for this pgen.")
  }

  gam_ <- gam(formula = formula, data = data_orig)
  sp_ <- gam_$sp

  # formula_fixed is what we will change to make sure sp is fixed
  formula_fixed <- formula
  if (length(sp_) == 0) {
    # We don't need to do anything:
    # formula_fixed will use the same as formula.
  } else if (length(sp_) == 1) {
    # This case is for, e.g., y ~ s(x1, x2)
    # we check that sure user did not already fix sp
    if (is.null(formula_fixed[[3]]$sp)) {
      formula_fixed[[3]]$sp <- sp_
    }
  } else if (length(sp_) == 2) {
    # This case is for, e.g., y ~ s(x1) + s(x2)
    formula_fixed[[3]][[2]]$sp <- sp_[[1]]
    formula_fixed[[3]][[3]]$sp <- sp_[[2]]
  } else {
    stop("We are not prepared for this formula type. ",
         "Please make a feature request.")
  }
  ret_fn <- function(data, newdata, se.fit) {
    gam_evaluated <- gam(formula = formula_fixed, data = data, ...)
    x.mesh <- as.vector(newdata) #evaluation mesh
    X0 <- predict(gam_evaluated, newdata, type = "lpmatrix") #values of linear predictor
    eps <- 1e-7 # finite difference interval
    shiftednewdata <- newdata + eps # shift the evaluation mesh
    X1 <- predict(gam_evaluated, shiftednewdata, type = "lpmatrix")
    if (deriv_order == 1) {
      Xp <- (X1-X0)/eps
      est <- Xp%*%coef(gam_evaluated) #predicted first derivative
      if (se.fit) {
        se <- rowSums(Xp%*%gam_evaluated$Vp*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
        pred_ <- data.frame(fit = est, se.fit = se)
      } else {
        pred_ <- data.frame(fit = est)
      }
    } else if (deriv_order == 2) {
      backshiftednewdata <- newdata - eps #back shift data
      X_1 <- predict(gam_evaluated, backshiftednewdata, type = "lpmatrix") #predicts on back shifted
      Xpp <- (X1 + X_1 - 2*X0)  / eps^2
      est <- Xpp%*%coef(gam_evaluated) #predicted second derivative
      if (se.fit) {
        se <- rowSums(Xpp%*%gam_evaluated$Vp*Xpp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
        pred_ <- data.frame(fit = est, se.fit = se)
      } else {
        pred_ <- data.frame(fit = est)
      }
    } else {
      stop("'deriv_order' must be '1' or '2'. The value '", deriv_order,
           "' is not valid.")
    }

    return(pred_)
  }
  return(ret_fn)
}


#' Estimation using locfit::locfit
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#'        The formula must use the verbose 'locfit::lp(x)', or the user must
#'        run "library(locfit)" before using this pgen.
#' @param ... dots are passed to locfit.
#' @param deriv_order order of derivative.
#'
#' @importFrom locfit locfit
#' @export
pgen_loc <- function(data_orig, formula, deriv_order = 0, ...) {
  if (length(deriv_order) > 0) {
    if (deriv_order > 2) {
      stop("This order of derivative not supported for this pgen.")
    }
    if (deriv_order == 2) {
      deriv_order <- c(1, 1)
    } else {
      if (deriv_order == 0) {
        deriv_order = numeric(0)
      }
    }
  }
  ret_fn <- function(data, newdata, se.fit) {
    lp <-locfit(formula, data = data, deriv = deriv_order, ...)
    pred_<- predict(lp, newdata = newdata, se.fit = se.fit)
    return(pred_)
  }
  return(ret_fn)
}


#' Estimation using stats::smooth.spline
#'
#' @param data_orig The user's original data passed to the client function.
#' @param formula The user's original formula passed to the client function.
#' @param ... dots are passed to spline.
#' @param deriv_order order of derivative.
#'
#' @importFrom stats smooth.spline
#' @export
pgen_spline <- function(data_orig, formula, deriv_order, ...) {
  if (deriv_order == 0) {
    # level estimation has poor performance in some of the cases we tested, so
    # until we figure out why we only enable predictors of derivative.
    # TODO {sims} reference simulations where we show the poor performance.
    stop("Spline predictor only supported for deriv_order > 0.")
  }
  f <- formula(formula)
  vars <- all.vars(f)
  ret_fn <- function(data, newdata, se.fit) {
    # Make new explanatory dataframe from newdata
    new_x <- newdata[, vars[2:length(vars)]]
    # x-variable
    x <- data[, vars[2:length(vars)]]
    # y-variable
    y <-  data[, vars[1]]
    # Fit and Predict
    sm_ <- smooth.spline(x, y, ...)
    fit <- predict(sm_, x = new_x, deriv = deriv_order)$y
    if (se.fit) {
      pred_ <- data.frame(fit = fit, se.fit = NA)
    } else {
      pred_ <- data.frame(fit = fit)
    }
    return(pred_)
  }
  return(ret_fn)
}
