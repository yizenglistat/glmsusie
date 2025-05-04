#' glmcs: Generalized Linear Model with Confidence Sets
#'
#' @description
#' Fits a sparse likelihood-based additive single-effect (LASER) model 
#' using the iteratively blockwise coordinate ascent alogirthm
#' across various generalized linear model families. Supports Gaussian, binomial, 
#' Poisson, Gamma GLMs, and Cox regression models.
#'
#' @param X Numeric matrix or vector of predictors. If a vector, it will be converted to a matrix.
#' @param y Response variable (vector for GLMs, matrix with time and status for Cox regression).
#' @param L Integer specifying the number of components to fit (must be at least 1 and less than ncol(X)).
#' @param family A family object specifying the model type (e.g., \code{gaussian()}, \code{binomial()},
#'        \code{poisson()}, \code{Gamma()}, or \code{cox()}).
#' @param coverage Numeric in \code{[0,1]} specifying the target coverage probability for credible sets (default: 0.95).
#' @param cor_threshold Numeric in \code{[0,1]} specifying the minimum absolute correlation for variables 
#'        within a credible set (default: 0.5).
#' @param standardize Logical indicating whether to standardize predictors (default: TRUE).
#' @param null_threshold Numeric specifying the threshold below which coefficients are set to zero (default: 1e-6).
#' @param tol Numeric specifying convergence tolerance for log-likelihood (default: 5e-2).
#' @param lambda Numeric specifying convergence tolerance for log-likelihood (default: 5e-2).
#' @param tau Numeric specifying convergence tolerance for log-likelihood (default: 5e-2).
#' @param ties String specifying method for handling tied events in Cox regression: "efron" (default) or "breslow".
#' @param max_iter Integer specifying maximum number of fitting iterations (default: 100).
#' @param seed Integer seed for reproducibility (default: NULL).
#'
#' @return A list with components:
#' \describe{
#'   \item{theta}{Matrix of estimated coefficients (p by L)}
#'   \item{intercept}{Estimated intercept (NULL for Cox regression)}
#'   \item{pmp}{Matrix of posterior model probabilities (p by L)}
#'   \item{dispersion}{Estimated dispersion parameter}
#'   \item{loglik}{Matrix of log-likelihoods (p by L)}
#'   \item{cs}{Credible sets obtained from posterior probabilities}
#'   \item{niter}{Number of iterations performed}
#' }
#'
#' @export
glmcs <- function(X, y, L=10L,
                family = gaussian(),
                coverage = 0.95,
                cor_threshold = 0.5,
                standardize = TRUE,
                null_threshold = 1e-5,
                tol = 5e-2,
                ties = c("efron", "breslow"),
                lambda = 0.0,
                tau = 0.5,
                max_iter = 100L,
                seed = NULL) {
 
 # Set random seed if provided
 if (!is.null(seed)) set.seed(seed)
 
 # Input validation and preprocessing
 if (is.vector(X)) X <- as.matrix(X, ncol = 1)
 if (!is.matrix(X)) stop("X must be a matrix or coercible to a matrix")
 
 if (is.character(family)) {
     stop("Family must be a family object (e.g., gaussian(), binomial(), poisson(), Gamma()) or 'cox'")
 }
 
 if (!is.numeric(coverage) || coverage <= 0 || coverage >= 1) {
   stop("coverage must be between 0 and 1")
 }
 
 if (!is.numeric(cor_threshold) || cor_threshold < 0 || cor_threshold > 1) {
   stop("cor_threshold must be between 0 and 1")
 }
 
 if (!is.logical(standardize)) {
   stop("standardize must be TRUE or FALSE")
 }
 
 # Handle Cox regression y format
 if (family$family == "cox") {
   if (is.vector(y)) {
     stop("For Cox regression, y must be a matrix with columns (time, status)")
   } else if (!is.matrix(y) || ncol(y) != 2) {
     stop("For Cox regression, y must be a matrix with columns (time, status)")
   }
 } else {
   if (is.matrix(y)) y <- drop(y)
 }
 
 # Match ties method for Cox regression
 ties <- match.arg(ties)
 
 # Fit the model
 output <- additive_effect_fit(
   X = X, 
   y = y, 
   L = L, 
   family = family, 
   standardize = standardize, 
   null_threshold = null_threshold, 
   tol = tol, 
   ties = ties, 
   lambda = lambda,
   tau = tau,
   max_iter = max_iter
 )
 
 # Calculate credible sets
 if (!is.null(output$pmp) && !is.null(output$kept)) {
   # Calculate correlation matrix if needed for credible sets
   R <- NULL
   if (cor_threshold > 0) {
     R <- stats::cor(X)
   }
   
   cs <- confidence_set(
     pmp = output$pmp,
     kept = output$kept,
     coverage = coverage, 
     Rmat = R,
     cor_threshold = cor_threshold
   )
   
   output$cs <- cs
 }

 output$marginal <- 1 - apply(X = 1 - output$pmp[, output$kept, drop = FALSE], MARGIN = 1, FUN = prod)
 
 return(output)
}