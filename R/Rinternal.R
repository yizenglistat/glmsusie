# Internal functions for glmsusie package

#' @useDynLib glmsusie, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Compute Log-Likelihood for Univariate GLM
#'
#' @description 
#' Computes the log-likelihood for a univariate generalized linear model (GLM)
#' with a single covariate, intercept, and optional offset. Supports any family from
#' \code{stats::family()} (e.g.\ Gaussian, Binomial, Poisson).
#'
#' @details
#' The function forms the linear predictor
#' \code{η = intercept + offset + θ * x}, then uses the family's
#' \code{linkinv} and \code{dev.resids} functions to compute deviance
#' residuals \code{d_i}.  The log-likelihood is
#' \deqn{\ell = -\,\sum_i d_i / 2,} profiling out dispersion for families
#' that estimate it (Gaussian, Gamma, inverse.gaussian) and leaving it fixed
#' for others (Binomial, Poisson).
#'
#' @param x         Numeric vector of length \code{n}: the single covariate.
#' @param y         Numeric vector of length \code{n}: response values (for
#'                  Binomial, can be 0/1 or a two-column matrix of counts).
#' @param family    An R \code{family} object (from \code{stats::family()},
#'                  default \code{gaussian()}).
#' @param theta     Numeric scalar: coefficient at which to evaluate the
#'                  log-likelihood.
#' @param offset    Numeric vector of length \code{n}, or scalar, default 0:
#'                  optional offset in the linear predictor.
#' @param intercept Numeric scalar: intercept term in the linear predictor, default 0.
#'
#' @return
#' Numeric scalar: the (profiled) log-likelihood at \code{theta} and \code{intercept}.
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' y <- 1.5 + 2 * x + rnorm(100)
#' univariate_loglik_glm(x, y, family = gaussian(), theta = 2, intercept = 1.5)
#'
#' x <- runif(200, -2, 2)
#' p <- plogis(0.5 + 1.5 * x)
#' y_bin <- rbinom(200, 1, p)
#' univariate_loglik_glm(x, y_bin, family = binomial(), theta = 1.5, intercept = 0.5)
#'
#' x_p <- rpois(150, lambda = 2)
#' mu <- exp(-1 + 0.3 * x_p)
#' y_p <- rpois(150, mu)
#' univariate_loglik_glm(x_p, y_p, family = poisson(), theta = 0.3, intercept = -1)
#' }
#' @name univariate_loglik_glm
#' @export
NULL

#' Fit a Univariate GLM Without Intercept Using IRLS
#'
#' @description 
#' Fits a generalized linear model (GLM) using a single covariate and no intercept term
#' via the iteratively reweighted least squares (IRLS) algorithm.
#'
#' @details
#' The linear predictor is defined as: 
#' \deqn{\eta = \theta \cdot x + \mathrm{offset}}{
#' η = θ * x + offset
#' }
#' where \eqn{g(\mu) = \eta} is the canonical link function from the specified GLM family.
#'
#' The function uses R's \code{family} object (e.g., \code{binomial()}, \code{poisson()}, \code{gaussian()})
#' to evaluate the inverse link, variance, and derivative functions at each iteration.
#'
#' Numerical safeguards are included:
#' \itemize{
#'   \item The linear predictor \code{eta} is clamped for Poisson models to prevent overflow.
#'   \item Variance and derivative evaluations are floored to avoid division by zero.
#'   \item Extremely large weights are capped.
#' }
#'
#' @param x Numeric vector of length \code{n}: predictor values.
#' @param y Numeric vector of length \code{n}: response values.
#' @param family An R \code{family} object, such as \code{binomial()}, \code{poisson()}, or \code{gaussian()}.
#' @param offset Numeric vector of length 1 or \code{n}, default \code{0}. Optional offset in the linear predictor.
#' @param max_iter Integer. Maximum number of IRLS iterations (default = 25).
#' @param tol Numeric. Convergence tolerance (default = 1e-8).
#'
#' @return 
#' Numeric scalar: estimated slope \eqn{\theta}.
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' y <- rbinom(100, 1, plogis(2 * x))
#' univariate_irls_glm_no_intercept(x, y, family = binomial())
#'
#' x <- rnorm(100)
#' y <- 1 + 3 * x + rnorm(100)
#' univariate_irls_glm_no_intercept(x, y, family = gaussian())
#' }
#' 
#' @seealso [glm()], [stats::family()]
#' @name univariate_irls_glm_no_intercept
#' @export
NULL

#' Compute Log-Likelihood for Univariate Cox Model
#'
#' @description 
#' Computes the partial log-likelihood for a univariate Cox proportional hazards model.
#' Handles ties using either the Breslow or Efron approximation.
#'
#' @details
#' The function computes the linear predictor as: \code{lp = offset + theta * x}.
#' The partial log-likelihood is calculated using either the Breslow or Efron method
#' for handling tied event times.
#'
#' @param x Numeric vector of length \code{n}: covariate values for each individual.
#' @param y Numeric matrix of shape \code{n × 2}, where:
#'        - \code{y[,1]} contains event/censoring times
#'        - \code{y[,2]} contains event status (1 = event, 0 = censored)
#' @param theta Numeric scalar: the coefficient to evaluate the log-likelihood at.
#' @param offset Numeric vector of length \code{n}, or scalar. Optional offset in the linear predictor.
#' @param ties Character string: tie-handling method. Must be either \code{"breslow"} or \code{"efron"}.
#'
#' @return Numeric scalar: the partial log-likelihood value.
#'
#' @examples
#' \dontrun{
#' x <- c(1, 2, 3, 4)
#' y <- matrix(c(4,1, 1,1, 3,0, 2,1), ncol = 2, byrow = TRUE)
#' univariate_loglik_cox(x, y, offset = 0, theta = 0.5, ties = "efron")
#' }
#' @name univariate_loglik_cox
#' @export
NULL

#' Compute Univariate Log-Likelihood for GLM or Cox Model
#'
#' @description
#' Calculates the log-likelihood of a single covariate effect under either a
#' generalized linear model (GLM) family or a Cox proportional hazards model.
#'
#' @param x Numeric vector of length n: covariate values.
#' @param y For GLMs: numeric vector of length n; for Cox: numeric matrix with
#'   two columns (time, status) and n rows.
#' @param family A GLM family object (e.g. \code{gaussian()}, \code{binomial()},
#'   \code{poisson()}) or a Cox family list with element \code{family="cox"}.
#' @param theta Numeric scalar: coefficient at which to evaluate the log-likelihood.
#' @param offset Numeric scalar or vector of length n: offset in the linear predictor.
#' @param intercept Numeric scalar: intercept term in the linear predictor (for GLM only).
#'   Default is 0.
#' @param ties Character string: tie-handling method for Cox partial likelihood;
#'   one of \code{"efron"} or \code{"breslow"}.
#'
#' @return
#' A numeric scalar giving the log-likelihood of the univariate fit.
#'
#' @examples
#' \dontrun{
#' # Gaussian example
#' x <- rnorm(100)
#' y <- 1.5 + 2 * x + rnorm(100)
#' univariate_loglik(x, y, family = gaussian(), theta = 2, intercept = 1.5)
#'
#' # Binomial example
#' x <- rnorm(200)
#' eta <- -1 + 1.5 * x
#' p <- plogis(eta)
#' y <- rbinom(200, 1, p)
#' univariate_loglik(x, y, family = binomial(link = "logit"),
#'                   theta = 1.5, intercept = -1)
#'
#' # Cox example (note: intercept not used in Cox models)
#' x <- rnorm(50)
#' times <- rexp(50)
#' status <- rbinom(50, 1, 0.5)
#' ymat <- cbind(times, status)
#' univariate_loglik(x, ymat, family = list(family = "cox"),
#'                   theta = 0.5, offset = 0, ties = "efron")
#' }
#' @name univariate_loglik
#' @export
NULL

#' Estimate Univariate Cox Model via Iteratively Reweighted Least Squares (IRLS)
#'
#' @description
#' Fits a univariate Cox proportional hazards model by maximizing the partial
#' log-likelihood using an Iteratively Reweighted Least Squares (IRLS) approach.
#' Supports both Breslow and Efron approximations for handling ties.
#'
#' @details
#' The function starts from an initial coefficient value of 0 and updates the
#' slope estimate using Newton-Raphson iterations until convergence or until
#' reaching the maximum number of iterations. The linear predictor is
#' \code{lp = offset + theta * x}. The score function and observed information
#' are used to update the estimate.
#'
#' @param x Numeric vector of length \code{n}: covariate values.
#' @param y Numeric matrix with shape \code{n × 2}, where:
#'   - \code{y[,1]} is the observed time
#'   - \code{y[,2]} is the event indicator (1 = event, 0 = censored)
#' @param offset Numeric scalar or vector of length \code{n}. Optional offset for the linear predictor.
#' @param ties Character string specifying the method to handle ties: \code{"breslow"} (default) or \code{"efron"}.
#' @param lambda Numeric penalty weight; if ≤ 0, defaults to \eqn{\sqrt{2\log(n)/n}} (default: 0.0).
#' @param tau Numeric truncation parameter; if ≤ 0, defaults to 0.5 (default: 0.5).
#' @param max_iter Integer: maximum number of IRLS iterations. Default is 25.
#' @param tol Numeric: convergence tolerance on parameter change. Default is 1e-8.
#'
#' @return A numeric scalar representing the estimated regression coefficient \code{theta}.
#'
#' @examples
#' \dontrun{
#' x <- c(1, 2, 3, 4)
#' y <- matrix(c(4,1, 1,1, 3,0, 2,1), ncol = 2, byrow = TRUE)
#' univariate_irls_cox(x, y, offset = 0, ties = "efron", max_iter = 50, tol = 1e-6)
#' }
#' @name univariate_irls_cox
#' @export
NULL

#' Compute Generalized Linear Model Estimate for Univariate Predictor with Intercept
#'
#' @description
#' Fits a GLM with intercept and single predictor using iteratively reweighted least squares (IRLS).
#' Implements the model:
#' \deqn{g(E[y]) = \beta_0 + \beta_1 x + \text{offset}}
#' where \eqn{g} is the link function determined by the selected family.
#'
#' @param x Numeric vector of covariates (length n).
#' @param y Numeric response vector (length n).
#' @param family A stats::family object (e.g. \code{gaussian()},
#'   \code{binomial()}, \code{poisson()}).
#' @param offset Numeric scalar or vector (length n) giving the linear predictor offset.
#' @param max_iter Integer. Maximum number of IRLS iterations.
#' @param tol Numeric convergence tolerance for parameters.
#'
#' @return A list with components:
#'   \item{intercept}{Numeric scalar: the estimated intercept \eqn{\hat\beta_0}.}
#'   \item{theta}{Numeric scalar: the estimated coefficient \eqn{\hat\beta_1}.}
#'
#' @details
#' This function implements the standard IRLS algorithm for GLMs, optimized for the
#' univariate predictor case with intercept. It provides identical results to R's
#' \code{glm(y ~ x, family=family)} but is more computationally efficient for
#' the single-predictor case.
#'
#' Numerical stability measures are implemented for handling extreme values, including
#' special cases for Poisson regression and fallback methods for matrix solution if
#' the standard solve fails.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 100
#' x <- rnorm(n)
#' 
#' # Gaussian example
#' y_gaussian <- 2 + 0.5*x + rnorm(n)
#' result <- univariate_irls_glm(x, y_gaussian, gaussian(), offset=numeric(0))
#' print(result)
#' 
#' # Poisson example
#' eta <- 1 + 0.5*x
#' y_poisson <- rpois(n, exp(eta))
#' result <- univariate_irls_glm(x, y_poisson, poisson(), offset=numeric(0))
#' print(result)
#' 
#' # Binomial example
#' prob <- 1/(1 + exp(-(0.5 + 0.8*x)))
#' y_binomial <- rbinom(n, 1, prob)
#' result <- univariate_irls_glm(x, y_binomial, binomial(), offset=numeric(0))
#' print(result)
#' }
#' @name univariate_irls_glm
#' @export
NULL

#' Compute Univariate Fit with Optional Truncated‐L1 Penalty
#'
#' @description
#' Fits a single‐covariate model (GLM or Cox) with optional standardization and
#' a truncated‐L1 penalty on the coefficient.  For GLMs it uses IRLS with a
#' capped‐L1 update; for Cox it uses a penalized IRLS on the partial likelihood.
#'
#' @param x Numeric vector of covariate values (length n); a scalar expands to zeros.
#' @param y Response:
#'   - For GLMs: numeric vector of length n.
#'   - For Cox: numeric matrix with 2 columns (time, status) and n rows.
#' @param family A stats::family object (e.g. \code{gaussian()}, \code{binomial()}, \code{poisson()})
#'   or a Cox family list with \code{family = "cox"}.
#' @param offset Numeric scalar or vector (length n) giving the linear predictor offset (default: 0).
#' @param standardize Logical: if TRUE, center and scale \code{x} before fitting (default: TRUE).
#' @param ties Character: ties method for Cox partial likelihood ("efron" or "breslow", default: "efron").
#' @param lambda Numeric penalty weight; if \code{NULL} or ≤ 0, defaults to \eqn{\sqrt{2\log(n)/n}}.
#' @param tau Numeric truncation parameter; if \code{NULL} or ≤ 0, defaults to 0.5.
#'
#' @return A list with elements:
#'   \item{intercept}{Estimated intercept (after undoing standardization).}
#'   \item{theta}{Estimated coefficient (after undoing standardization).}
#'   \item{loglik}{Unpenalized log-likelihood at the estimated \code{theta}.}
#'   \item{bic}{Bayesian Information Criterion: \eqn{-2*loglik + 2\log(n)}.}
#'
#' @examples
#' \dontrun{
#' set.seed(101)
#' n <- 50
#' x <- rnorm(n)
#' # Gaussian GLM
#' y_gauss <- 1.5*x + rnorm(n)
#' res1 <- univariate_fit(x, y_gauss, family = gaussian(), offset = 0)
#'
#' # Cox example
#' times  <- rexp(n, rate = exp(0.7*x))
#' status <- rbinom(n, 1, 0.6)
#' y_cox  <- cbind(time=times, status=status)
#' res2 <- univariate_fit(x, y_cox, family = list(family="cox"))
#' }
#' @name univariate_fit
#' @export
NULL

#' Compute Single‐Effect Fit for GLM or Cox with Optional Truncated‐L1 Penalty
#'
#' @description
#' Applies a single‐effect model to each column of a predictor matrix using either
#' GLM (via IRLS + truncated‐L1) or Cox partial likelihood (via penalized IRLS).
#'
#' @param X Numeric matrix (n × p) of predictors. Each column is fit separately.
#' @param y Response:
#'   - For GLMs: numeric vector of length n.
#'   - For Cox: numeric matrix with 2 columns (time, status) and n rows.
#' @param family A stats::family object (e.g. \code{gaussian()}, \code{binomial()}, 
#'   \code{poisson()}) or a Cox family list with \code{family = "cox"}.
#' @param offset Numeric scalar or vector (length n) giving the linear predictor offset (default: 0).
#' @param standardize Logical: if TRUE, center and scale each predictor column before fitting (default: TRUE).
#' @param shrinkage Logical: if TRUE, shrinkage parameters in fitting (default: TRUE).
#' @param ties Character: ties method for Cox partial likelihood ("efron" or "breslow", default: "efron").
#' @param lambda Numeric penalty weight; if ≤ 0, defaults to \eqn{\sqrt{2\log(n)/n}} (default: 0.0).
#' @param tau Numeric truncation parameter; if ≤ 0, defaults to 1.0 (default: 1.0).
#' @param alpha level of significance
#'
#' @return A list of length p, where each element is itself a list with components:
#'   \item{loglik}{Unpenalized log‐likelihood at the fitted coefficient.}
#'   \item{bic}{Bayesian Information Criterion: \eqn{-2*logLik + 2\log(n)}.}
#'   \item{bic_diff}{BIC difference from null model.}
#'   \item{bf}{Bayes factor.}
#'   \item{pmp}{Posterior model probability.}
#'   \item{intercept}{Estimated intercept for that predictor.}
#'   \item{theta}{Estimated coefficient for that predictor.}
#'   \item{pval_intercept}{P-value for estimated intercept for that predictor.}
#'   \item{pval_theta}{P-value for estimated coefficient for that predictor.}
#'   \item{expect_intercept}{Expected value of intercept under model averaging.}
#'   \item{expect_theta}{Expected value of coefficient under model averaging.}
#'   \item{expect_variance}{Expected variance under model averaging.}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 80; p <- 3
#' X <- matrix(rnorm(n*p), n, p)
#' # Gaussian example
#' y_gauss <- X[,1] * 1.2 + rnorm(n)
#' res_glm <- single_effect_fit(X, y_gauss, family = gaussian())
#'
#' # Cox example
#' times <- rexp(n, rate = exp(0.5 * X[,2]))
#' status <- rbinom(n, 1, 0.7)
#' y_cox <- cbind(time=times, status=status)
#' res_cox <- single_effect_fit(X, y_cox, family = list(family="cox"))
#' }
#' @name single_effect_fit
#' @export
NULL

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
#'   - For GLMs: numeric vector of length n.
#'   - For Cox: numeric matrix with 2 columns (time, status) and n rows.
#' @param L Integer; number of single effects to include (will be set to \code{min(10, L)}).
#' @param family A stats::family object (e.g. \code{gaussian()}, \code{binomial()}, 
#'   \code{poisson()}) or a Cox family list with \code{family = "cox"}.
#' @param standardize Logical: if TRUE, center and scale each predictor column before fitting (default: TRUE).
#' @param ties Character: ties method for Cox partial likelihood ("efron" or "breslow", default: "efron").
#' @param lambda Numeric penalty weight; if ≤ 0, defaults to \eqn{\sqrt{2\log(n)/n}} (default: 0.0).
#' @param tau Numeric truncation parameter; if ≤ 0, defaults to 0.5 (default: 0.5).
#' @param decompose Logical: if TRUE, decompose the theta in fitting (default: TRUE).
#' @param shrinkage Logical: if TRUE, shrinkage parameters in fitting (default: TRUE).
#' @param alpha level of significance
#' @param tol Numeric; convergence tolerance on the change in expected log-likelihood (default: 5e-2).
#' @param max_iter Integer; maximum number of coordinate-ascent iterations (default: 100).
#'
#' @return A list with components:
#'   \item{niter}{Number of iterations performed.}
#'   \item{loglik}{p × L matrix of univariate log-likelihoods.}
#'   \item{expect_loglik}{Vector of length \code{niter} giving the expected log-likelihood at each iteration.}
#'   \item{final_loglik}{Expected log-likelihood at convergence.}
#'   \item{intercept}{Estimated intercept term.}
#'   \item{dispersion}{Estimated dispersion parameter (for Gaussian/Gamma).}
#'   \item{theta}{p × L matrix of fitted single-effect coefficients.}
#'   \item{pval_intercept}{p × L matrix of p-values for fitted single-effect intercepts}
#'   \item{pval_theta}{p × L matrix of p-values for fitted single-effect coefficients.}
#'   \item{pmp}{p × L matrix of posterior model probabilities.}
#'   \item{bic}{p × L matrix of BIC values.}
#'   \item{bic_diff}{p × L matrix of BIC differences from null.}
#'   \item{bf}{p × L matrix of Bayes factors.}
#'   \item{expect_variance}{Length-L vector of PMP-weighted variances.}
#'   \item{elapsed_time}{Numeric; total computation time in seconds.}
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
#' @name additive_effect_fit
#' @export
NULL

#' Decompose Row Sums of a Theta Matrix into Single-Effect Components
#'
#' @description
#' Given a \eqn{p \times L} matrix \code{theta}, this function redistributes
#' the row sums into a new matrix of the same shape where each row's total
#' effect is placed into a single column (using round-robin assignment).
#'
#' This operation is useful for decomposing a combined effect vector into
#' sparse single-effect columns (e.g., in spike-and-slab regression or
#' additive effect models).
#'
#' @param theta A numeric matrix of shape \code{p × L}.
#' @param L Integer. Number of columns in the output matrix.
#'
#' @return A numeric matrix of shape \code{p × L}, where each row sum matches
#' the corresponding row sum of \code{theta}, and each row has only one nonzero entry.
#'
#' @examples
#' \dontrun{
#' theta <- matrix(rnorm(12), nrow = 6, ncol = 2)
#' theta_new <- decompose_theta(theta, L = 2)
#' all.equal(rowSums(theta), rowSums(theta_new))  # Should be TRUE
#' }
#'
#' @name decompose_theta
#' @export
NULL