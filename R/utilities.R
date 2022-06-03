cast2vec <- function(obj) {
  # Attempt to cast various objects to a vector.
  # Currently supported inputs:
  # - a matrix with one column
  # - a one-dimensional array

  if (is.vector(obj)) {
    vec <- obj
  } else if (is.matrix(obj) && ncol(obj) == 1) {
    vec <- as.vector(obj)
  } else if (is.array(obj) && length(dim(obj)) == 1) {
    # e.g., mgcv::predict.gam()
    vec <- as.vector(obj)
  } else if (is.data.frame(obj) && ncol(obj) == 1) {
    vec <- obj[, 1]
  } else {
    stop("'obj' cannot be cast to a vector.")
  }
  return(vec)
}


#' @importFrom stats terms
validate_formula <-function(formula, data) {
  if (any(!(labels(terms(formula)) %in% colnames(data)))) {
       stop("X variable not found in data")
  }
  if (length(all.vars(formula)) == length(labels(terms(formula)))) {
    #  indicates a regression (validate presence of y)
    if (!(all.vars(formula)[1] %in% colnames(data))) {
      stop("Y variable not found in data")
    }
  }
}
