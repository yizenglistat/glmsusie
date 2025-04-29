#' LASER: Likelihood-based Additive Single-Effect Regression algorithm for generalized linear regression problems
#'
#' glmcs fits the LASER model via blockwise coordinate ascent, perform variable selection under multicollinearity,
#' and build confidence sets at a significant level. glmcs supports general regression problems such as GLM built-in families and cox regression.
#'
#' @param X A numeric \code{n x p} predictor matrix.
#' @param y A response: numeric vector (for GLM) or \code{n x 2} matrix (time, status) for Cox.
#' @param family A GLM family object (e.g.\ \code{gaussian("identity")}, \code{binomial()})
#'   or the string \code{"cox"}.
#' @param L Integer number of single-effect components (default \code{10L}
#'   or \code{min(10, p)} if \code{p < 10}).
#' @param coverage Desired posterior mass for confidence sets (default \code{0.95}).
#' @param standardize Logical; should predictors be standardized? (default \code{TRUE}).
#' @param ties Character; tie-handling method for Cox: \code{"efron"} or \code{"breslow"}
#'   (default \code{"efron"}).
#' @param method Character; coordinate-ascent update order:
#'   \code{"greedy"} or \code{"shuffle"} (default \code{"greedy"}).
#' @param max_iter Integer; maximum number of coordinate-ascent iterations
#'   (default \code{100L}).
#' @param step_size Numeric; multiplicative step-size for updates (default \code{1.0}).
#' @param min_abs_corr Minimum absolute correlation threshold for purity filtering
#'   (default \code{0}).
#' @param tol Numeric; convergence tolerance on log-likelihood change
#'   (default \code{1e-6}).
#' @param seed Optional integer seed for reproducibility (default \code{NULL}).
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{elapsed}}{Total fitting time (seconds).}
#'   \item{\code{fit}}{Raw output from \code{\link{get_laser_fit}}.}
#'   \item{\code{pmp}}{Posterior marginal probabilities of inclusion.}
#'   \item{\code{theta}}{Estimated effect vector (sum over \eqn{L} components).}
#'   \item{\code{cs}}{Confidence sets (list of integer vectors) from
#'     \code{\link{get_cs}}.}
#'   \item{\code{coverage}}{Requested coverage level.}
#'   \item{\code{keep}}{Logical vector of length \eqn{L}, indicating significant effects
#'     via \code{\link{get_included}}.}
#'   \item{\code{X}, \code{y}, \code{family}, \code{ties}, \code{method},
#'     \code{step_size}, \code{seed}}{Echoed inputs.}
#' }
#'
#' @seealso
#'  * \code{\link{get_laser_fit}} for the low‐level C++ fit routine
#'  * \code{\link{get_included}} to select significant effects
#'  * \code{\link{get_cs}} to build confidence sets
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 100; p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(3, -2, rep(0, p - 2))
#' y <- X %*% beta + rnorm(n)
#'
#' res <- glmcs(
#'   X            = X,
#'   y            = y,
#'   family       = gaussian("identity"),
#'   L            = 5L,
#'   coverage     = 0.9,
#'   standardize  = TRUE,
#'   method       = "greedy",
#'   seed         = 123
#' )
#' print(res$theta)
#' print(res$cs)
#' }
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
                  ties          = "efron",
                  method        = "greedy",
                  max_iter      = 100L,
                  step_size     = 1.0,
                  min_abs_corr  = 0.00,
                  tol           = 1e-6,
                  seed          = NULL) 
{
  #--- 1) Fit LASER model -----------------------------------------------------
  fit <- get_laser_fit(
    X           = X,
    y           = y,
    family      = family,
    L           = L,
    standardize = standardize,
    ties        = ties,
    method      = method,
    max_iter    = max_iter,
    step_size   = step_size,
    tol         = tol,
    seed        = seed
  )

  #--- 2) Time elapsed --------------------------------------------------------
  elapsed <- fit$total_elapsed

  #--- 3) Which effects are significant? -------------------------------------
  keep <- get_included(
    fit    = fit,
    X      = X,
    y      = y,
    family = family,
    alpha  = tol  # you may choose a separate alpha if preferred
  )

  #--- 4) Confidence sets + purity filtering ----------------------------------
  cs <- get_cs(
    posterior_mat = fit$posterior,
    keep          = keep,
    coverage      = coverage,
    X             = NULL,
    R             = cor(X),
    min_abs_corr  = min_abs_corr
  )

  #--- 5) Posterior marginal probabilities -----------------------------------
  #    P(include) = 1 − ∏ℓ (1 − p_{j,ℓ})
  pmp <- 1 - apply(
    X      = 1 - fit$posterior[, keep, drop = FALSE],
    MARGIN = 1,
    FUN    = prod
  )

  #--- 6) Assemble output -----------------------------------------------------
  out <- list(
    elapsed     = elapsed,
    fit         = fit,
    pmp         = pmp,
    theta       = fit$theta_hat,
    cs          = cs,
    coverage    = coverage,
    keep        = keep,
    X           = X,
    y           = y,
    family      = family,
    ties        = ties,
    method      = method,
    step_size   = step_size,
    seed        = seed
  )
  class(out) <- "glmcs"
  return(out)
}