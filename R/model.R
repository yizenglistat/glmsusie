#' Generalized Linear Models with Confidence Sets (glmsusie)
#'
#' @description
#' Fits a sparse Likelihood-based Additive Single-Effect Regression (LASER) model
#' using iteratively blockwise coordinate ascent. The model represents the coefficient
#' vector as a sum of sparse "single effects" and produces confidence sets for variable
#' selection based on posterior model probabilities.
#'
#' @details
#' The LASER model decomposes the coefficient vector into a sum of L sparse components.
#' At each iteration, the algorithm cyclically updates one component while holding
#' the others fixed. For each component, it fits a univariate model for each predictor,
#' computes model probabilities (via BIC), and updates coefficients as probability-weighted
#' averages. The approach extends traditional GLMs by providing Bayesian-inspired confidence
#' sets for variable selection.
#'
#' Supported model families include:
#' \itemize{
#'   \item Gaussian linear regression
#'   \item Binomial logistic regression
#'   \item Poisson regression 
#'   \item Gamma regression
#'   \item Other GLM family regression
#'   \item Cox proportional hazards regression
#' }
#'
#' @param X Numeric matrix of predictors (n × p). If a vector is provided, it will 
#'        be converted to a single-column matrix.
#' @param y Response variable:
#'   \itemize{
#'     \item For GLMs: numeric vector of length n
#'     \item For Cox models: numeric matrix with 2 columns (time, status) and n rows
#'   }
#' @param L Integer specifying the number of components to fit (default: 10).
#'        Will be truncated to min(10, ncol(X)) if necessary.
#' @param family A family object specifying the model type (e.g., \code{gaussian()}, 
#'        \code{binomial()}, \code{poisson()}, \code{Gamma()}) or \code{list(family="cox")}
#'        for Cox regression.
#' @param coverage Numeric in \eqn{[0,1]} specifying the target coverage probability 
#'        for confidence sets (default: 0.95).
#' @param cor_threshold Numeric in \eqn{[0,1]} specifying the minimum absolute correlation 
#'        required for variables to be grouped in the same confidence set (default: 0.5).
#' @param standardize Logical indicating whether to center and scale predictors 
#'        before fitting (default: TRUE).
#' @param decompose Logical indicating whether to decompose theta in fitting (default: TRUE).
#' @param shrinkage Logical indicating whether to shrinkage parameters using pvals (default: TRUE).
#' @param tol Numeric specifying convergence tolerance for the expected log-likelihood
#'        between iterations (default: 5e-2).
#' @param lambda Numeric penalty weight for the truncated-L1 penalty (default: 0.0).
#'        If 0, no penalization is applied.
#' @param tau Numeric truncation parameter for the truncated-L1 penalty (default: 1e-5).
#'        Controls the transition from L1 to L0 regularization.
#' @param ties Character string specifying the method for handling tied event times 
#'        in Cox regression: "efron" (default) or "breslow".
#' @param max_iter Integer specifying maximum number of coordinate ascent iterations
#'        (default: 100).
#' @param seed Integer seed for reproducibility (default: NULL).
#'
#' @return A list with class "glmsusie" containing:
#' \describe{
#'   \item{call}{The matched call}
#'   \item{X}{The model matrix}
#'   \item{y}{The response vector/matrix}
#'   \item{family}{The family object used}
#'   \item{theta}{p × L matrix of estimated coefficients for each single effect}
#'   \item{intercept}{Estimated intercept (NULL for Cox regression)}
#'   \item{pmp}{p × L matrix of posterior model probabilities}
#'   \item{loglik}{p × L matrix of log-likelihoods}
#'   \item{bic}{p × L matrix of BIC values}
#'   \item{bic_diff}{p × L matrix of BIC differences from null model}
#'   \item{evidence}{p × L matrix of evidence}
#'   \item{bf}{p × L matrix of Bayes factors}
#'   \item{marginal}{Vector of marginal inclusion probabilities for each predictor}
#'   \item{kept}{Logical vector indicating which effects were retained}
#'   \item{cs}{List of confidence sets based on posterior probabilities}
#'   \item{niter}{Number of iterations performed}
#'   \item{max_iter}{Number of maximum iterations}
#'   \item{elapsed}{Elapsed computation time in seconds}
#' }
#'
#' @examples
#' \dontrun{
#' # Gaussian linear regression example
#' set.seed(42)
#' n <- 100
#' p <- 50
#' X <- matrix(rnorm(n*p), n, p)
#' colnames(X) <- paste0("X", 1:p)
#' true_beta <- c(rep(1, 5), rep(0, p-5))
#' y <- X %*% true_beta + rnorm(n)
#' 
#' # Fit model with 3 components
#' fit <- glmsusie(X, y, L = 3, family = gaussian())
#' 
#' # Examine results
#' summary(fit)
#' plot(fit)
#' 
#' # Extract coefficients
#' coef(fit)
#' coef(fit, intercept = TRUE)
#' 
#' # Cox regression example
#' X <- matrix(rnorm(100*10), 100, 10)
#' colnames(X) <- paste0("X", 1:10)
#' times <- rexp(100, rate = exp(0.5 * X[,1] + 0.5 * X[,2]))
#' status <- rbinom(100, 1, 0.7)
#' y_cox <- cbind(times, status)
#' 
#' fit_cox <- glmsusie(X, y_cox, L = 2, family = list(family = "cox"))
#' summary(fit_cox)
#' }
#'
#' @seealso
#' \code{\link{summary.glmsusie}} for summarizing model results,
#' \code{\link{coef.glmsusie}} for extracting coefficients,
#' \code{\link{plot.glmsusie}} for plotting results
#'
#' @export
glmsusie <- function(X, y, L = 10L,
                  family = gaussian(),
                  coverage = 0.95,
                  cor_threshold = 0.5,
                  standardize = TRUE,
                  decompose = TRUE,
                  shrinkage = TRUE,
                  tol = 5e-2,
                  lambda = 0.0,
                  tau = 1e-5,
                  ties = c("efron", "breslow"),
                  max_iter = 500L,
                  seed = NULL) {
  
  # Capture the call
  cl <- match.call()

  # Set random seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Input validation and preprocessing
  if (is.vector(X)) X <- as.matrix(X, ncol = 1)
  if (!is.matrix(X)) stop("X must be a matrix or coercible to a matrix")
  
  if (is.character(family)) {
    stop("Family must be a family object (e.g., gaussian(), binomial(), poisson(), Gamma()) or list(family='cox')")
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

  # Store column names if available
  col_names <- colnames(X)
  if (is.null(col_names)) {
    col_names <- paste0("X", 1:ncol(X))
  }
  
  # Handle Cox regression y format
  is_cox <- FALSE
  if (is.list(family) && !is.null(family$family) && family$family == "cox") {
    is_cox <- TRUE
    if (is.vector(y)) {
      stop("For Cox regression, y must be a matrix with columns (time, status)")
    } else if (!is.matrix(y) || ncol(y) != 2) {
      stop("For Cox regression, y must be a matrix with columns (time, status)")
    }
  } else {
    if (is.matrix(y)) y <- drop(y)
  }
  
  if(L > ncol(X)) L <- ncol(X)

  # Match ties method for Cox regression
  ties <- match.arg(ties)
  
  # Fit the model
  out <- additive_effect_fit(
    X = X, 
    y = y, 
    L = L, 
    family = family, 
    standardize = standardize, 
    decompose = decompose, 
    shrinkage = shrinkage,
    tol = tol, 
    ties = ties, 
    lambda = lambda,
    alpha = 1 - coverage,
    tau = tau,
    max_iter = max_iter
  )

  kept <- iskept(out$pmp, out$pval_theta, 1-coverage)

  # Collect everything we'll need later
  result <- list(
    call = cl,
    X = X,
    y = y,
    family = family,
    theta = out$theta,
    intercept = sum(out$intercept),
    pmp = out$pmp,
    loglik = out$loglik,
    bic = out$bic,
    bic_diff = out$bic_diff,
    evidence = out$evidence,
    pval_wald = out$pval_wald,
    std_err = out$std_err,
    bf = out$bf,
    pval_intercept = out$pval_intercept,
    pval_theta = out$pval_theta,
    kept = kept,
    niter = out$niter,
    max_iter = max_iter,
    elapsed = out$elapsed_time
  )
  
  # Add row and column names
  rownames(result$theta) <- col_names
  rownames(result$pmp) <- col_names
  colnames(result$theta) <- paste0("Effect", 1:ncol(result$theta))
  colnames(result$pmp) <- paste0("Effect", 1:ncol(result$pmp))
  colnames(result$pval_intercept) <- paste0("Effect", 1:ncol(result$pval_intercept))
  colnames(result$pval_theta) <- paste0("Effect", 1:ncol(result$pval_theta))

  
  # Calculate confidence sets
  if (!is.null(result$pmp) && !is.null(result$kept)) {
    # Calculate correlation matrix if needed for confidence sets
    R <- NULL
    if (cor_threshold > 0) {
      R <- stats::cor(X)
    }
    
    cs <- confidence_set(
      pmp = result$pmp,
      kept = result$kept,
      coverage = coverage, 
      Rmat = R,
      cor_threshold = cor_threshold
    )
    
    result$cs <- cs
    result$coverage <- coverage
  }

  # Calculate marginal pvalues
  combine_simes <- function(pvals) {
    L <- length(pvals)
    pvals_sorted <- sort(pvals)
    min(L * pvals_sorted / seq_len(L))
  }

  combine_fisher <- function(pvals) {
    stat <- -2 * sum(log(pmax(pvals, 1e-300)))
    pchisq(stat, df = 2 * length(pvals), lower.tail = FALSE)
  }

  if(sum(kept) == 0){
    p <- ncol(X)
    result$pval_simes <- rep(1, p)
    result$pval_fisher <- rep(1, p)
  }else{
    result$pval_simes <- apply(result$pval_theta[, kept, drop=FALSE], 1, combine_simes)
    result$pval_fisher <- apply(result$pval_theta[, kept, drop=FALSE], 1, combine_fisher)
  }
  names(result$pval_simes) <- col_names
  names(result$pval_fisher) <- col_names
  
  # Calculate posterior inclusion probabilities
  if(sum(kept) == 0){
    p <- ncol(X)
    result$pip <- rep(1/p, p)
  }else{ 
    result$pip <- 1 - apply(
      X = 1 - result$pmp[, kept, drop=FALSE], 
      MARGIN = 1, 
      FUN = prod
    )
  }
  names(result$pip) <- col_names
  
  # Assign class and return
  class(result) <- "glmsusie"
  result
}