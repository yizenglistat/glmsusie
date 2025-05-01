#' Compute Univariate Log-Likelihood for GLM or Cox Model
#'
#' @description
#' Calculates the log-likelihood of a single covariate effect under either a
#' generalized linear model (GLM) family or a Cox proportional hazards model.
#'
#' @param x Numeric vector of covariate values (length \(n\)).
#' @param y Response:
#'   - For GLMs: numeric vector of length \(n\).
#'   - For Cox: numeric \eqn{n\times2} matrix with columns (time, status).
#' @param family A GLM family object (e.g. \code{gaussian()}, \code{binomial()},
#'   \code{poisson()}, etc.) or a Cox family list with \code{family = "cox"}.
#' @param theta Numeric scalar. The coefficient for the univariate model
#'   (default: 0).
#' @param offset Numeric scalar or vector of length \(n\). Known offset for the
#'   linear predictor (default: 0).
#' @param ties Character. Ties method for Cox partial likelihood:
#'   \code{"efron"} (default) or \code{"breslow"}.
#'
#' @return
#' A single numeric value: the log-likelihood of the univariate fit.
#'
#' @examples
#' ## Gaussian GLM
#' set.seed(1)
#' x <- rnorm(100)
#' y <- 2 * x + rnorm(100)
#' univariate_loglik(x, y, family = gaussian(), theta = 2)
#'
#' ## Binomial GLM (logit link)
#' set.seed(2)
#' x <- rnorm(200)
#' eta <- -1 + 1.5 * x
#' p  <- 1 / (1 + exp(-eta))
#' y  <- rbinom(200, size = 1, prob = p)
#' univariate_loglik(x, y, family = binomial(link = "logit"), theta = 1.5)
#'
#' @export
univariate_loglik <- function(x, y, family, theta=0, offset=0, ties="efron") 
{
	n <- length(x)
    if (length(offset)==1) offset <- rep(offset, n)
  
	if(family$family == "cox") {
    	loglik <- univariate_loglik_cox(x=x, y=y, theta=theta, offset=offset, ties=ties)
    }else{
        eta <- offset + theta * x
        mu  <- family$linkinv(eta)
        w   <- rep(1, n)
        loglik <- -sum(family$dev.resids(y, mu, w))/2
    }

    return(loglik)
}

#' Fit a Univariate Effect via IRLS for GLM or Cox Models
#'
#' @description
#' Estimates a single regression coefficient θ for one predictor by solving
#' the iteratively reweighted least squares (IRLS) equations.  Works for both
#' generalized linear models (via \code{univariate_irls_glm}) and Cox
#' proportional hazards (via \code{univariate_irls_cox}).
#'
#' @param x Numeric vector of length n.  The single predictor.
#' @param y Response:
#'   \itemize{
#'     \item Numeric vector (length n) for GLMs.
#'     \item Numeric matrix (n × 2) with columns (time, status) for Cox.
#'   }
#' @param family A GLM family object (e.g. \code{gaussian()}, \code{binomial()}, \code{poisson()}, …)
#'   or a list with \code{family = "cox"} for Cox regression.
#' @param offset Numeric scalar or vector of length n.  Known offset in the
#'   linear predictor (default: 0).
#' @param standardize Logical.  If \code{TRUE}, center and scale \code{x}
#'   before fitting (default: \code{TRUE}), so that θ is on the original scale.
#' @param ties Character.  Tie-handling method for Cox IRLS:
#'   \code{"efron"} (default) or \code{"breslow"}.  Ignored for GLMs.
#' @param null_threshold Numeric.  Estimates with absolute value ≤ this
#'   threshold are set to zero (default: 1e-6).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{theta}}{Estimated regression coefficient θ (on original scale).}
#'   \item{\code{loglik}}{Log-likelihood evaluated at θ.}
#'   \item{\code{bic}}{Bayesian Information Criterion,
#'     \eqn{-2*loglik + 2\,\log(n)}.}
#' }
#'
#' @seealso
#' \code{\link{glm}}, \code{\link[glmnet]{glmnet}}, \code{\link{univariate_irls_glm}},
#' \code{\link{univariate_loglik}}
#'
#' @examples
#' ## Gaussian example
#' set.seed(1)
#' x <- rnorm(100)
#' y <- 1 + 2 * x + rnorm(100)
#' univarite_fit(x, y, family = gaussian())
#'
#' ## Logistic regression
#' set.seed(2)
#' x <- rnorm(200)
#' eta <- -0.5 + 1.5 * x
#' p <- 1 / (1 + exp(-eta))
#' y <- rbinom(200, 1, p)
#' univarite_fit(x, y, family = binomial())
#'
#' @export
univarite_fit <- function(x, y, family,
                          offset = 0,
                          standardize = TRUE,
                          ties = "efron",
                          null_threshold = 1e-6)
{
  # validation
  if (is.matrix(y)) n <- nrow(y) else n <- length(y)
  if (length(x) == 1L)     x <- rep(0, n)
  if (length(offset) == 1L) offset <- rep(offset, n)

  # standardize predictor if requested
  if (standardize) {
    x_mean <- mean(x)
    x_norm <- sqrt(sum((x - x_mean)^2))
    if (x_norm == 0) {
      x_mean <- 0; x_norm <- 1; x_std <- x
    } else {
      x_std <- (x - x_mean) / x_norm
    }
  } else {
    x_mean <- 0; x_norm <- 1; x_std <- x
  }

  # fit via IRLS
  if (family$family == "cox") {
    if (all(x == 0)) {
      theta_std <- 0
    } else {
      theta_std <- univariate_irls_cox(
        x      = x_std,
        y      = y,
        offset = offset,
        ties   = ties
      )
    }
    theta <- theta_std / x_norm
  } else {
    if (all(x == 0)) {
      theta_std <- 0
    } else {
      theta_std <- univariate_irls_glm(
        x          = x_std,
        y          = y,
        family     = family,
        offset     = offset,
        dispersion = family$dispersion
      )
    }
    theta <- theta_std / x_norm
  }

  # apply threshold for sparsity
  theta <- ifelse(abs(theta) <= null_threshold, 0, theta)

  # compute loglik and BIC
  loglik <- univariate_loglik(x, y, family, theta, offset)
  bic    <- -2 * loglik + 2 * log(n)

  list(theta = theta, loglik = loglik, bic = bic)
}

#' Fit Single-Effect Regression Across All Predictors
#'
#' @description
#' For each column of the design matrix \code{X}, fits a univariate model to
#' the response \code{y} (using \code{\link{univarite_fit}}), computes BIC
#' differences relative to a null model, converts these into Bayes factors
#' and posterior model probabilities (PMP), and returns both the MAP estimates
#' and the PMP-weighted expectations.
#'
#' @param X Numeric matrix (n × p) of predictors.
#' @param y Response:
#'   \itemize{
#'     \item Numeric vector of length n, for GLMs.
#'     \item Numeric matrix (n × 2) with columns \code{time,status}, for Cox.
#'   }
#' @param family A GLM family object (e.g. \code{gaussian()}, \code{binomial()}, \code{poisson()}, …)
#'   or a list with \code{family = "cox"} for Cox regression.
#' @param offset Numeric scalar or vector of length n.  Known offset term.
#' @param standardize Logical; if \code{TRUE}, center and scale each column of \code{X}
#'   before fitting (default: \code{TRUE}).
#' @param ties Character; tie-handling method for Cox IRLS (\code{"efron"} or \code{"breslow"}),
#'   ignored for GLMs (default: \code{"efron"}).
#' @param null_threshold Numeric; posterior probabilities or coefficient estimates
#'   below this threshold are set to zero (default: \code{1e-6}).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{loglik}}{Numeric vector of length p; log-likelihood for each univariate fit.}
#'   \item{\code{bic}}{Numeric vector of length p; BIC for each fit.}
#'   \item{\code{bic_diff}}{Numeric vector; BIC difference from null model, shifted so min = 0.}
#'   \item{\code{bf}}{Numeric vector; Bayes factors \eqn{\exp(-½\,\Delta{\rm BIC})}.}
#'   \item{\code{pmp}}{Numeric vector; posterior model probabilities (sum to 1).}
#'   \item{\code{theta}}{Numeric vector of MAP coefficient estimates (thresholded).}
#'   \item{\code{expect_theta}}{Numeric vector; PMP-weighted expectation of θ.}
#'   \item{\code{expect_variance}}{Numeric scalar; PMP-weighted average of θ² across predictors.}
#' }
#'
#' @examples
#' \dontrun{
#' # Gaussian single-effect fit
#' X <- matrix(rnorm(200*5), 200, 5)
#' y <- rnorm(200)
#' se <- single_effect_fit(X, y, family = gaussian())
#'
#' # Logistic single-effect fit
#' eta <- X[,1] * 1.2
#' p <- 1/(1+exp(-eta))
#' y_bin <- rbinom(200, 1, p)
#' se_bin <- single_effect_fit(X, y_bin, family = binomial())
#' }
#' 
#' @export
single_effect_fit <- function(X, y, family, offset,
                              standardize = TRUE,
                              ties = "efron",
                              null_threshold = 1e-6) 
{
  n <- nrow(X); p <- ncol(X)
  theta    <- numeric(p)
  loglik   <- numeric(p)
  bic      <- numeric(p)
  bic_diff <- numeric(p)

  # null model BIC
  res_null <- univarite_fit(
    x             = rep(0, n),
    y             = y,
    family        = family,
    offset        = offset,
    standardize   = standardize,
    ties          = ties,
    null_threshold = null_threshold
  )

  for (j in seq_len(p)) {
    res <- univarite_fit(
      x             = X[, j],
      y             = y,
      family        = family,
      offset        = offset,
      standardize   = standardize,
      ties          = ties,
      null_threshold = null_threshold
    )
    theta[j]      <- res$theta
    loglik[j]     <- res$loglik
    bic[j]        <- res$bic
    bic_diff[j]   <- res$bic - res_null$bic
  }

  # Bayes factors & posterior model probabilities
  bic_diff  <- bic_diff - min(bic_diff)
  bf        <- exp(-0.5 * bic_diff)
  pmp       <- bf / sum(bf)

  # threshold small values
  theta            <- ifelse(pmp <= null_threshold, 0, theta)
  pmp              <- ifelse(pmp <= null_threshold, 0, pmp)
  expect_theta     <- pmp * theta
  expect_variance  <- mean(pmp * theta^2)

  expect_theta    <- ifelse(expect_theta    <= null_threshold, 0, expect_theta)
  expect_variance <- ifelse(expect_variance <= null_threshold, 0, expect_variance)

  list(
    loglik          = loglik,
    bic             = bic,
    bic_diff        = bic_diff,
    bf              = bf,
    pmp             = pmp,
    theta           = theta,
    expect_theta    = expect_theta,
    expect_variance = expect_variance
  )
}

#' Fit Likelihood-based Additive Single-Effect Regression (LASER) Model
#'
#' @description
#' Fits the LASER model, representing the coefficient vector as a sum of \code{L}
#' sparse "single effects."  At each iteration, it cyclically applies
#' \code{\link{single_effect_fit}} to update one effect while holding the others
#' fixed, then updates the intercept (and dispersion, if applicable), until
#' convergence or \code{max_iter} is reached.
#'
#' @param X Numeric matrix (n × p) of predictors.
#' @param y Response:
#'   \itemize{
#'     \item Numeric vector of length \code{n} for GLMs.
#'     \item Numeric matrix (n × 2) with \code{(time,status)} columns for Cox.
#'   }
#' @param L Integer; number of single effects to include (will be set to \code{min(10, L)}).
#' @param family GLM family object (e.g. \code{gaussian()}, \code{binomial()}, \code{poisson()})
#'   or list with \code{family="cox"} for Cox models.
#' @param standardize Logical; if \code{TRUE}, center and scale each column of \code{X}
#'   before fitting (default: \code{TRUE}).
#' @param null_threshold Numeric; any posterior probability or coefficient below this
#'   value is set to zero (default: \code{1e-6}).
#' @param tol Numeric; convergence tolerance on the change in expected log-likelihood
#'   (default: \code{5e-2}).
#' @param ties Character; tie-handling method for Cox IRLS (\code{"efron"} or \code{"breslow"}),
#'   ignored for GLMs (default: \code{"efron"}).
#' @param max_iter Integer; maximum number of coordinate-ascent iterations (default: \code{100}).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{niter}}{Number of iterations performed.}
#'   \item{\code{loglik}}{p × L matrix of univariate log-likelihoods.}
#'   \item{\code{expect_loglik}}{Vector of length \code{niter} giving the expected
#'     log-likelihood at each iteration.}
#'   \item{\code{final_loglik}}{Expected log-likelihood at convergence.}
#'   \item{\code{intercept}}{Estimated intercept term.}
#'   \item{\code{dispersion}}{Estimated dispersion parameter (for Gaussian/Gamma).}
#'   \item{\code{theta}}{p × L matrix of fitted single-effect coefficients.}
#'   \item{\code{pmp}}{p × L matrix of posterior model probabilities.}
#'   \item{\code{bic}}{p × L matrix of BIC values.}
#'   \item{\code{bic_diff}}{p × L matrix of BIC differences from null.}
#'   \item{\code{bf}}{p × L matrix of Bayes factors.}
#'   \item{\code{expect_variance}}{Length-L vector of PMP-weighted variances.}
#'   \item{\code{kept}}{Logical vector of length L; \code{TRUE} for effects retained.}
#'   \item{\code{elapsed_time}}{Numeric; total computation time in seconds.}
#' }
#'
#' @examples
#' \dontrun{
#' # Gaussian example with 5 effects
#' X <- matrix(rnorm(100*20), 100, 20)
#' y <- rnorm(100)
#' res <- additive_effect_fit(X, y, L = 5, family = gaussian())
#' 
#' # Cox regression example
#' times  <- rexp(100)
#' status <- rbinom(100, 1, 0.5)
#' y_cox  <- cbind(times, status)
#' res_cox <- additive_effect_fit(X, y_cox, L = 3,
#'                                family = list(family = "cox"),
#'                                ties = "breslow")
#' }
#'
#' @export
additive_effect_fit <- function(X, y, L,
                                family = gaussian(),
                                standardize = TRUE,
                                null_threshold = 1e-6,
                                tol = 5e-2,
                                ties = "efron",
                                max_iter = 100)
{ 
  # Start timing the algorithm
  start_time <- proc.time()
  
  n <- nrow(X); p <- ncol(X)
  L <- min(10, L)

  ## Determine if we should estimate dispersion
  if (is.character(family) || family$family == "cox") {
    family <- list(family = "cox", link = "log", dispersion = 1)
  }
  estimate_dispersion <- family$family %in% 
    c("gaussian", "Gamma", "inverse.gaussian", "quasibinomial", "quasipoisson")
  
  intercept <- 0
  family$dispersion <- 1
  offset <- 0
  theta <- matrix(0, p, L)

  loglik   <- matrix(0, p, L)
  bic      <- matrix(0, p, L)
  bic_diff <- matrix(0, p, L)
  bf       <- matrix(0, p, L)
  pmp      <- matrix(0, p, L)
  expect_variance <- numeric(L)

  expect_loglik <- numeric(max_iter)
  
  for (iter in seq_len(max_iter)) {
    # current linear predictor
    linear_predictor <- intercept + rowSums(X %*% theta)

    # update each single effect
    for (l in seq_len(L)) {
      offset <- linear_predictor - X %*% theta[, l]
      res <- single_effect_fit(X, y, family, offset,
                               standardize = standardize,
                               ties         = ties,
                               null_threshold = null_threshold)

      theta[, l]       <- ifelse(res$expect_theta <= null_threshold, 0, res$expect_theta)
      loglik[, l]      <- res$loglik
      bic[, l]         <- res$bic
      bic_diff[, l]    <- res$bic_diff
      bf[, l]          <- res$bf
      pmp[, l]         <- res$pmp
      expect_variance[l] <- res$expect_variance

      # refresh predictor for next effect
      linear_predictor <- offset + X %*% theta[, l]
    }

    # update intercept (GLM only)
    if (family$family != "cox") {
      intercept <- univariate_irls_glm(
        x      = rep(1, n),
        y      = y,
        family = family,
        offset = rowSums(X %*% theta)
      )
    }

    # update dispersion if needed
    if (estimate_dispersion) {
      family$dispersion <- update_dispersion(
        y        = y,
        family   = family,
        offset   = intercept + rowSums(X %*% theta),
        approach = "pearson"
      )
    }

    # compute expected log-likelihood
    expect_loglik[iter] <- mean(pmp * loglik)
    if (iter > 1 && diff(expect_loglik[iter + c(-1,0)]) < tol) break
  }

  kept <- expect_variance > tol
  
  # Calculate total elapsed time
  elapsed_time <- as.numeric((proc.time() - start_time)[3])

  list(
    niter           = iter,
    loglik          = loglik,
    expect_loglik   = expect_loglik[seq_len(iter)],
    final_loglik    = expect_loglik[iter],
    intercept       = intercept,
    dispersion      = family$dispersion,
    theta           = theta,
    pmp             = pmp,
    bic             = bic,
    bic_diff        = bic_diff,
    bf              = bf,
    expect_variance = expect_variance,
    kept            = kept,
    elapsed_time    = elapsed_time
  )
}