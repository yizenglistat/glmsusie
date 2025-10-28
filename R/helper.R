#' Run glmnet with Cross-Validation
#'
#' @description
#' Fits elastic net regularization using glmnet with cross-validation to select
#' the optimal regularization parameter. Formats output to be compatible with
#' gSuSiE result structure.
#'
#' @param X Design matrix of predictors (n × p).
#' @param y Response vector or survival object.
#' @param family Response family. Can be gaussian(), binomial(), poisson(), 
#'   Gamma(), or "cox" for Cox regression.
#' @param alpha Elastic net mixing parameter: 0 = ridge, 1 = lasso, 0 < alpha < 1 = elastic net.
#' @param standardize Logical. Should variables be standardized? Default TRUE.
#' @param nfolds Number of cross-validation folds. Default 10.
#'
#' @return A list containing:
#' \describe{
#'   \item{cs}{List with 'sets' component containing singleton credible sets for selected variables}
#'   \item{theta}{Vector of estimated coefficients (excluding intercept)}
#'   \item{pip}{Binary vector indicating variable selection (1 = selected, 0 = not selected)}
#'   \item{lambda_1se}{Lambda value within 1 standard error of minimum CV error}
#'   \item{lambda_min}{Lambda value that minimizes CV error}
#'   \item{selected}{Indices of selected variables}
#' }
#'
#' @details
#' This function serves as a wrapper around \code{glmnet::cv.glmnet} to provide
#' output compatible with gSuSiE results. Uses lambda.1se for coefficient estimation,
#' which provides more conservative variable selection than lambda.min.
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n <- 100; p <- 50
#' X <- matrix(rnorm(n * p), n, p)
#' y <- X[,1:3] %*% c(1, -1, 0.5) + rnorm(n)
#' 
#' # Run elastic net
#' result <- run_glmnet(X, y, alpha = 0.5)
#' print(result$selected)
#' }
#'
#' @importFrom glmnet cv.glmnet
#' @export
run_glmnet <- function(X, y, family = gaussian(), alpha = 1, 
                       standardize = TRUE, nfolds = 10) {
  
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required but not installed.")
  }
  
  p <- ncol(X)
  
  # Handle family conversion
  if (is.character(family) && family == "cox") {
    glmnet_family <- "cox"
  } else if (!is.character(family)) {
    glmnet_family <- family$family
  } else {
    glmnet_family <- family
  }
  
  # For Cox regression, need to handle differently
  if (glmnet_family == "cox") {
    # Assuming y is a Surv object or matrix with time and status
    cv_fit <- glmnet::cv.glmnet(X, y, family = "cox", alpha = alpha, 
                                standardize = standardize, nfolds = nfolds)
  } else {
    cv_fit <- glmnet::cv.glmnet(X, y, family = glmnet_family, alpha = alpha,
                                standardize = standardize, nfolds = nfolds)
  }
  
  # Get coefficients at lambda.1se (more conservative) or lambda.min
  coef_1se <- as.numeric(coef(cv_fit, s = "lambda.1se"))
  
  # For Cox, coef() doesn't include intercept
  if (glmnet_family == "cox") {
    theta <- coef_1se
  } else {
    theta <- coef_1se[-1]  # Remove intercept
  }
  
  # Find selected variables (non-zero coefficients)
  selected <- which(theta != 0)
  
  # Create credible sets (singleton sets for each selected variable)
  cs_sets <- if (length(selected) > 0) {
    lapply(selected, function(x) x)  # Each selected variable forms its own set
  } else {
    list()  # Empty list if no variables selected
  }
  
  # Create binary PIPs (1 for selected, 0 for unselected)
  # Alternative: could set to NA since glmnet doesn't produce true PIPs
  pip <- rep(0, p)
  pip[selected] <- 1
  
  return(list(
    cs = list(sets = cs_sets),
    theta = theta,
    pip = pip,
    lambda_1se = cv_fit$lambda.1se,
    lambda_min = cv_fit$lambda.min,
    selected = selected
  ))
}

#' Run LASSO Regression with Cross-Validation
#'
#' @description
#' Convenience wrapper for LASSO regression (L1 regularization) using glmnet
#' with cross-validation. Equivalent to \code{run_glmnet} with \code{alpha = 1}.
#'
#' @param X Design matrix of predictors (n × p).
#' @param y Response vector or survival object.
#' @param family Response family. Can be gaussian(), binomial(), poisson(), 
#'   Gamma(), or "cox" for Cox regression.
#' @param standardize Logical. Should variables be standardized? Default TRUE.
#' @param nfolds Number of cross-validation folds. Default 10.
#'
#' @return A list containing the same components as \code{run_glmnet}.
#'
#' @details
#' LASSO (Least Absolute Shrinkage and Selection Operator) performs both
#' regularization and variable selection by setting some coefficients to exactly zero.
#' This function provides a convenient interface specifically for LASSO regression.
#'
#' @examples
#' \dontrun{
#' # Simulate sparse data
#' set.seed(123)
#' n <- 100; p <- 50
#' X <- matrix(rnorm(n * p), n, p)
#' y <- X[,1:3] %*% c(1, -1, 0.5) + rnorm(n)
#' 
#' # Run LASSO
#' lasso_result <- run_lasso(X, y)
#' print(lasso_result$selected)
#' }
#'
#' @seealso \code{\link{run_glmnet}}, \code{\link{run_elastic_net}}
#' @export
run_lasso <- function(X, y, family = gaussian(), standardize = TRUE, nfolds = 10) {
  run_glmnet(X, y, family = family, alpha = 1, standardize = standardize, nfolds = nfolds)
}

#' Run Elastic Net Regression with Cross-Validation
#'
#' @description
#' Convenience wrapper for elastic net regression using glmnet with cross-validation.
#' Combines L1 and L2 penalties, balancing variable selection and grouping effects.
#'
#' @param X Design matrix of predictors (n × p).
#' @param y Response vector or survival object.
#' @param family Response family. Can be gaussian(), binomial(), poisson(), 
#'   Gamma(), or "cox" for Cox regression.
#' @param alpha Elastic net mixing parameter. Default 0.5 (equal L1/L2 weighting).
#' @param standardize Logical. Should variables be standardized? Default TRUE.
#' @param nfolds Number of cross-validation folds. Default 10.
#'
#' @return A list containing the same components as \code{run_glmnet}.
#'
#' @details
#' Elastic net regression combines the variable selection capability of LASSO
#' with the grouping effect of ridge regression. When predictors are correlated,
#' elastic net tends to select groups of correlated variables rather than
#' arbitrarily choosing one from each group.
#'
#' @examples
#' \dontrun{
#' # Simulate correlated predictors
#' set.seed(123)
#' n <- 100; p <- 50
#' X <- matrix(rnorm(n * p), n, p)
#' # Add correlation between first 5 variables
#' X[,2:5] <- X[,2:5] + 0.8 * X[,1]
#' y <- X[,1:3] %*% c(1, -1, 0.5) + rnorm(n)
#' 
#' # Run elastic net
#' enet_result <- run_elastic_net(X, y, alpha = 0.5)
#' print(enet_result$selected)
#' }
#'
#' @seealso \code{\link{run_glmnet}}, \code{\link{run_lasso}}
#' @export
run_elastic_net <- function(X, y, family = gaussian(), alpha = 0.5, 
                           standardize = TRUE, nfolds = 10) {
  run_glmnet(X, y, family = family, alpha = alpha, standardize = standardize, nfolds = nfolds)
}