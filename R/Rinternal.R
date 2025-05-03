#' Estimate Dispersion for GLM Family
#'
#' @description
#' Computes the dispersion parameter for a Generalized Linear Model (GLM) using either the Pearson residual method or deviance-based approach.
#'
#' @details
#' This function supports any \code{family} object from \code{stats::family()}, such as \code{gaussian()}, \code{poisson()}, \code{Gamma()}, etc.
#' The linear predictor is given by \code{offset}, and the inverse link function is applied to compute the mean response \code{mu}.
#'
#' If \code{approach = "pearson"}, it computes:
#' \deqn{ \phi = \frac{1}{n - p} \sum_i \left( \frac{y_i - \mu_i}{\sqrt{\text{var}(\mu_i)}} \right)^2 }
#'
#' If \code{approach = "deviance"}, it computes:
#' \deqn{ \phi = \frac{1}{n - p} \sum_i \text{dev}_i }
#'
#' where \code{dev_i} are the deviance residuals.
#'
#' @param y Numeric response vector of length \code{n}.
#' @param family A GLM family object (e.g., \code{gaussian()}, \code{poisson()}, \code{Gamma()}).
#' @param offset Numeric vector or scalar representing the linear predictor \code{eta}; defaults to 0.
#' @param approach Character string: either \code{"pearson"} or \code{"deviance"}.
#'
#' @return A numeric scalar representing the estimated dispersion.
#'
#' @examples
#' \dontrun{
#' y <- rgamma(100, shape = 2, rate = 2)
#' offset <- rep(log(mean(y)), 100)
#' update_dispersion(y, Gamma(), offset = offset, approach = "pearson")
#' }
#'
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @name update_dispersion
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
#' @param offset Numeric vector of length \code{n}, or scalar. Optional offset in the linear predictor.
#' @param theta Numeric scalar: the coefficient to evaluate the log-likelihood at.
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
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
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
#' y <- 2 * x + rnorm(100)
#' univariate_loglik(x, y, family = gaussian(), theta = 2, offset = 0)
#'
#' # Binomial example
#' x <- rnorm(200)
#' eta <- -1 + 1.5 * x
#' p <- plogis(eta)
#' y <- rbinom(200, 1, p)
#' univariate_loglik(x, y, family = binomial(link = "logit"),
#'                   theta = 1.5, offset = 0)
#'
#' # Cox example
#' x <- rnorm(50)
#' times <- rexp(50)
#' status <- rbinom(50, 1, 0.5)
#' ymat <- cbind(times, status)
#' univariate_loglik(x, ymat, family = list(family = "cox"),
#'                   theta = 0.5, offset = 0, ties = "efron")
#' }
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
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
#' @param lambda Numeric penalty weight; if ≤0 defaults to \eqn{\sqrt{2\log(n)/n}}.
#' @param tau Numeric truncation parameter; if ≤0 uses grid \{1/n,…,5/n\} for n≤500 and 1/n otherwise.
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
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @name univariate_irls_cox
#' @export
NULL


#' Compute Truncated‐L1 IRLS Estimate for Univariate GLM
#'
#' @description
#' Fits a single‐covariate GLM using IRLS with a truncated‐L1 penalty.
#' Solves:
#' \deqn{\min_\theta\;-\ell(\theta) \;+\;\lambda\,\min\bigl(1,|\theta|/\tau\bigr)}
#' via iteratively reweighted least squares, with a closed‐form capped‐L1 update each step.
#'
#' @param x Numeric vector of covariates (length n).
#' @param y Numeric response vector (length n).
#' @param family A stats::family object (e.g. \code{gaussian()},
#'   \code{binomial()}, \code{poisson()}).
#' @param offset Numeric scalar or vector (length n) giving the linear predictor offset.
#' @param lambda Numeric penalty weight; if ≤0 defaults to \eqn{\sqrt{2\log(n)/n}}.
#' @param tau Numeric truncation parameter; if ≤0 uses grid \{1/n,…,5/n\} for n≤500 and 1/n otherwise.
#' @param max_iter Integer. Maximum number of IRLS iterations.
#' @param tol Numeric convergence tolerance for θ.
#'
#' @return Numeric scalar: the estimated coefficient \eqn{\hat\theta}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 100
#' x <- rnorm(n)
#' # Gaussian example
#' y <- 2*x + rnorm(n)
#' univariate_irls_glm(x, y, gaussian(), offset=0,
#'                     lambda=0, tau=0, max_iter=25, tol=1e-8)
#' }
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
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
#' @param null_threshold Numeric threshold below which the final \code{theta} is set to zero (default: 1e-6).
#' @param max_iter Integer: maximum number of IRLS iterations (default: 25).
#' @param tol Numeric convergence tolerance on \code{theta} updates (default: 1e-8).
#'
#' @return A list with elements:
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
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
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
#' @param ties Character: ties method for Cox partial likelihood ("efron" or "breslow", default: "efron").
#' @param lambda Numeric penalty weight; if ≤ 0, defaults to \eqn{\sqrt{2\log(n)/n}} (default: 0.0).
#' @param tau Numeric truncation parameter; if ≤ 0, defaults to 1.0 (default: 1.0).
#' @param null_threshold Numeric threshold below which an estimated coefficient is set to zero (default: 1e-6).
#'
#' @return A list of length p, where each element is itself a list with components:
#'   \item{theta}{Estimated coefficient for that predictor.}
#'   \item{loglik}{Unpenalized log‐likelihood at the fitted coefficient.}
#'   \item{bic}{Bayesian Information Criterion: \eqn{-2*logLik + 2\log(n)}.}
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
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
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
#' @param null_threshold Numeric threshold below which an estimated coefficient is set to zero (default: 1e-6).
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
#'   \item{pmp}{p × L matrix of posterior model probabilities.}
#'   \item{bic}{p × L matrix of BIC values.}
#'   \item{bic_diff}{p × L matrix of BIC differences from null.}
#'   \item{bf}{p × L matrix of Bayes factors.}
#'   \item{expect_variance}{Length-L vector of PMP-weighted variances.}
#'   \item{kept}{Logical vector of length L; \code{TRUE} for effects retained.}
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
#' @useDynLib glmcs, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @name additive_effect_fit
#' @export
NULL
