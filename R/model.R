#' Likelihood-based Additive Single-Effect Regression for Generalized Linear Models
#'
#' @description
#' Fits the LASER (Likelihood-based Additive Single-Effect Regression) model via 
#' block coordinate ascent to perform variable selection under multicollinearity 
#' and build confidence sets at a specified coverage level. Supports both generalized
#' linear models (GLMs) and Cox proportional hazards regression.
#'
#' @details
#' The LASER approach decomposes the coefficient vector into a sum of L sparse
#' effects, each fitted using a Bayesian model averaging approach. This implementation:
#' 
#' 1. Fits the LASER model using block coordinate ascent
#' 2. Identifies statistically significant effects
#' 3. Computes confidence sets with desired posterior coverage
#' 4. Filters confidence sets based on correlation purity (if requested)
#' 5. Calculates posterior marginal probabilities of inclusion
#' 
#' The algorithm is particularly effective for high-dimensional, multicollinear
#' data where traditional methods struggle to identify relevant variables.
#'
#' @param X A numeric matrix with n rows (observations) and p columns (predictors).
#' @param y Response variable:
#'   - For GLMs: a numeric vector of length n
#'   - For Cox: a matrix with n rows and 2 columns (time, status)
#' @param family A GLM family object (e.g., \code{gaussian()}, \code{binomial()},
#'   \code{poisson()}, \code{Gamma()}, \code{inverse.gaussian()}) or the 
#'   string \code{"cox"} for Cox proportional hazards.
#' @param L Integer. Number of single-effect components to include in the model.
#'   Default: 10, or min(10, p) if p < 10.
#' @param coverage Numeric in (0,1). Desired posterior mass for confidence sets.
#'   Default: 0.95.
#' @param standardize Logical. Whether to standardize predictors before fitting.
#'   Default: TRUE.
#' @param ties String. Tie-handling method for Cox regression: "efron" (more accurate) 
#'   or "breslow" (faster). Default: "efron".
#' @param algorithm String. Coordinate ascent update strategy:
#'   - "shuffle": Update effects in random order each iteration (default)
#'   - "greedy": Select update giving maximum improvement
#'   - "cyclic": Update effects in fixed order
#' @param max_iter Integer. Maximum number of coordinate ascent iterations.
#'   Default: 100.
#' @param step_size Numeric > 0. Multiplicative step-size for updates. Values < 1
#'   provide more conservative updates. Default: 1.0.
#' @param min_abs_corr Numeric in [0,1). Minimum absolute correlation threshold
#'   for purity filtering of confidence sets. Default: 0 (no filtering).
#' @param tol Numeric > 0. Convergence tolerance on log-likelihood change.
#'   Default: 1e-6.
#' @param seed Integer or NULL. Random seed for reproducibility of random operations.
#'   Default: NULL (no seed).
#'
#' @return A list of class "glmcs" with components:
#' \describe{
#'   \item{\code{elapsed}}{Numeric. Total fitting time in seconds.}
#'   \item{\code{fit}}{List. Raw output from \code{\link{get_laser_fit}}.}
#'   \item{\code{pmp}}{Numeric vector. Posterior marginal probabilities of inclusion.}
#'   \item{\code{theta}}{Numeric vector. Estimated effect vector (sum over significant components).}
#'   \item{\code{cs}}{List. Confidence sets (list of integer vectors) from \code{\link{get_cs}}.}
#'   \item{\code{coverage}}{Numeric. Requested coverage level.}
#'   \item{\code{keep}}{Logical vector. Indicators for significant effects from \code{\link{get_included}}.}
#'   \item{\code{X}, \code{y}, \code{family}, \code{ties}, \code{algorithm},
#'     \code{step_size}, \code{seed}}{The original input parameters.}
#'   \item{\code{convergence}}{Logical. Whether the algorithm converged within max_iter.}
#' }
#'
#' @seealso
#' \code{\link{get_laser_fit}} for the low-level C++ fit routine,
#' \code{\link{get_included}} to select significant effects,
#' \code{\link{get_cs}} to build confidence sets
#'
#' @examples
#' # Generate synthetic data with 2 true variables
#' set.seed(42)
#' n <- 100; p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' X[,1:2] <- X[,1:2] %*% matrix(c(1, 0.9, 0.9, 1), 2, 2) # Create correlation
#' beta <- c(3, 0, rep(0, p-2)) # Only first variable is truly active
#' y <- X %*% beta + rnorm(n)
#'
#' # Fit with LASER model
#' res <- glmcs(
#'   X           = X,
#'   y           = y,
#'   family      = gaussian(),
#'   L           = 5L,
#'   coverage    = 0.9,
#'   standardize = TRUE,
#'   algorithm   = "greedy"
#' )
#' 
#' # Examine results
#' print(res$theta)       # Estimated coefficients
#' print(res$cs)          # Confidence sets
#' print(res$pmp)         # Posterior marginal probabilities
#' 
#' # Logistic regression example
#' y_bin <- rbinom(n, 1, plogis(X %*% beta))
#' res_bin <- glmcs(X, y_bin, family = binomial())
#'
#' @importFrom stats cor gaussian
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
glmcs <- function(X,
                  y,
                  family        = gaussian("identity"),
                  L             = 10L,
                  coverage      = 0.95,
                  standardize   = TRUE,
                  ties          = c("efron", "breslow"),
                  algorithm     = c("shuffle", "greedy", "cyclic"),
                  max_iter      = 100L,
                  step_size     = 1.0,
                  min_abs_corr  = 0.0,
                  tol           = 1e-6,
                  seed          = NULL) 
{
  # Validate inputs
  if (!is.matrix(X)) {
    stop("X must be a numeric matrix")
  }
  
  if (!is.numeric(L) || L <= 0) {
    warning("L must be a positive integer. Setting to default: min(10, ncol(X))")
    L <- min(10L, ncol(X))
  } else {
    L <- as.integer(min(L, ncol(X)))
  }
  
  if (!is.numeric(coverage) || coverage <= 0 || coverage >= 1) {
    warning("coverage must be between 0 and 1. Setting to default: 0.95")
    coverage <- 0.95
  }
  
  if (!is.numeric(step_size) || step_size <= 0) {
    warning("step_size must be positive. Setting to default: 1.0")
    step_size <- 1.0
  }
  
  if (!is.numeric(min_abs_corr) || min_abs_corr < 0 || min_abs_corr >= 1) {
    warning("min_abs_corr must be in [0,1). Setting to default: 0")
    min_abs_corr <- 0.0
  }
  
  if (!is.numeric(tol) || tol <= 0) {
    warning("tol must be positive. Setting to default: 1.0")
    tol <- 1e-6
  }
  
  # Match tie handling method for Cox models
  ties <- match.arg(ties)
  
  # Match algorithm type
  algorithm <- match.arg(algorithm)
  
  # Ensure seed is either NULL or integer
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    set.seed(seed)
  }
  
  # Record start time for timing
  start_time <- proc.time()[3]
  
  # 1) Fit LASER model
  tryCatch({
    fit <- get_laser_fit(
      X           = X,
      y           = y,
      family      = family,
      L           = L,
      standardize = standardize,
      ties        = ties,
      algorithm   = algorithm,
      max_iter    = max_iter,
      step_size   = step_size,
      tol         = tol,
      seed        = seed
    )
  }, error = function(e) {
    stop("Error in LASER model fitting: ", e$message)
  })

  # 2) Record elapsed time
  elapsed <- proc.time()[3] - start_time

  # 3) Determine which effects are significant
  tryCatch({
    keep <- get_included(
      fit    = fit,
      X      = X,
      y      = y,
      family = family,
      alpha  = 0.05  # Significance level for effect inclusion
    )
  }, error = function(e) {
    warning("Error identifying significant effects: ", e$message)
    keep <- rep(TRUE, L)  # Default to keeping all effects if there's an error
  })
  
  # Early exit with simplified result if no effects are significant
  if (!any(keep)) {
    message("No significant effects found. Returning simplified result.")
    out <- list(
      elapsed     = elapsed,
      fit         = fit,
      pmp         = rep(0, nrow(fit$theta)),
      theta       = rep(0, nrow(fit$theta)),
      cs          = NULL,
      coverage    = coverage,
      keep        = keep,
      X           = X,
      y           = y,
      family      = family,
      ties        = ties,
      algorithm   = algorithm,
      step_size   = step_size,
      seed        = seed,
      convergence = (fit$n_iter < max_iter)
    )
    class(out) <- "glmcs"
    return(out)
  }

  # 4) Compute confidence sets with purity filtering
  tryCatch({
    # Compute correlation matrix for purity filtering if min_abs_corr > 0
    R <- if (min_abs_corr > 0) cor(X) else NULL
    
    cs <- get_cs(
      posterior_mat = fit$posterior,
      keep          = keep,
      coverage      = coverage,
      X             = NULL,
      R             = R,
      min_abs_corr  = min_abs_corr
    )
  }, error = function(e) {
    warning("Error computing confidence sets: ", e$message)
    cs <- NULL
  })

  # 5) Calculate posterior marginal probabilities
  # P(include) = 1 − ∏ℓ (1 − p_{j,ℓ})
  tryCatch({
    pmp <- 1 - apply(
      X      = 1 - fit$posterior[, keep, drop = FALSE],
      MARGIN = 1,
      FUN    = prod
    )
  }, error = function(e) {
    warning("Error computing posterior marginal probabilities: ", e$message)
    pmp <- rep(NA, nrow(fit$theta))
  })

  # 6) Assemble output
  out <- list(
    elapsed     = elapsed,
    fit         = fit,
    pmp         = pmp,
    theta       = rowSums(fit$theta[, keep, drop = FALSE]),
    cs          = cs,
    coverage    = coverage,
    keep        = keep,
    X           = X,
    y           = y,
    family      = family,
    ties        = ties,
    algorithm   = algorithm,
    step_size   = step_size,
    seed        = seed,
    convergence = (fit$n_iter < max_iter)
  )
  class(out) <- "glmcs"
  return(out)
}