#’ @useDynLib glmcs, .registration=TRUE
#’ @importFrom Rcpp evalCpp
NULL


#' Cox Proportional Hazards "family" Object
#'
#' @description
#' Creates a \code{family} object for univariate Cox regression, compatible
#' with GLM‐style interfaces.  Only the log link is supported, yielding the
#' usual Cox partial‐likelihood formulation.
#'
#' @param link Character string; link function for the Cox model.  Only
#'   \code{"log"} is supported (default).
#'
#' @details
#' This "family" object provides the minimal set of functions required
#' for use with IRLS‐based routines and log-likelihood calculations:
#' \itemize{
#'   \item \code{linkfun} and \code{linkinv} implement the log link.
#'   \item \code{mu.eta} provides the derivative of the inverse link.
#'   \item \code{variance}, \code{dev.resids}, and \code{aic} are placeholders
#'         (not used in Cox but needed for compatibility).
#' }
#'
#' @return
#' A \code{family} object (a list with class \code{"family"}) containing:
#' \describe{
#'   \item{\code{family}}{Always \code{"cox"}.}
#'   \item{\code{link}}{Link name, \code{"log"}.}
#'   \item{\code{linkfun}}{Function transforming \code{mu} to \code{eta}.}
#'   \item{\code{linkinv}}{Inverse link, mapping \code{eta} to \code{mu}.}
#'   \item{\code{mu.eta}}{Derivative of \code{linkinv}.}
#'   \item{\code{variance}}{Variance function (returns 1).}
#'   \item{\code{dev.resids}}{Deviance residuals (zeros).}
#'   \item{\code{aic}}{AIC placeholder (returns \code{-2}).}
#'   \item{\code{validmu}, \code{valideta}}{Validation functions (always TRUE).}
#'   \item{\code{dispersion}}{Always 1.}
#' }
#'
#' @examples
#' # Create the Cox family object
#' fam <- cox(link = "log")
#' stopifnot(fam$family == "cox")
#' # Check that linkinv(exp(eta)) == exp(eta)
#' eta <- c(0, 1, -1)
#' all.equal(fam$linkinv(eta), exp(eta))
#'
#' @export
cox <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  }

  if (linktemp == "log") {
    linkfun  <- function(mu)       log(mu)
    linkinv  <- function(eta)      exp(eta)
    mu.eta   <- function(eta)      exp(eta)
  } else {
    stop("Link '", linktemp, "' not recognized for Cox regression family")
  }

  valideta   <- function(eta) TRUE
  validmu    <- function(mu) all(mu > 0)
  variance   <- function(mu) rep(1, length(mu))
  dev.resids <- function(y, mu, wt) rep(0, length(mu))
  aic        <- function(y, n, mu, wt, dev) -2

  structure(
    list(
      family    = "cox",
      link      = linktemp,
      linkfun   = linkfun,
      linkinv   = linkinv,
      variance  = variance,
      dev.resids= dev.resids,
      aic       = aic,
      mu.eta    = mu.eta,
      validmu   = validmu,
      valideta  = valideta,
      dispersion= 1
    ),
    class = "family"
  )
}

#' Fit a Univariate GLM via IRLS
#'
#' @description
#' Efficiently fits a univariate generalized linear model (GLM) for a single
#' predictor using iteratively reweighted least squares (IRLS), with a fixed
#' offset and known dispersion parameter. Returns only the slope estimate.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values.
#' @param family A GLM family object (e.g., \code{gaussian()}, \code{binomial()},
#'   \code{poisson()}, etc.). Must inherit from class \code{"family"}.
#' @param offset Numeric scalar or vector of known offsets (default: \code{0}).
#' @param dispersion Numeric scalar dispersion parameter (default: \code{1}).
#'   Used for two-parameter families (Gaussian, Gamma, inverse Gaussian).
#' @param max_iter Integer. Maximum number of IRLS iterations (default: \code{25}).
#' @param tol Numeric. Convergence tolerance on the change in slope (default: \code{1e-8}).
#'
#' @details
#' This function solves for the slope coefficient \eqn{\beta} in the model
#' \deqn{\eta = \text{offset} + \beta\,x,\quad \mu = g^{-1}(\eta),}
#' where \eqn{g} is the link function of the specified family.  At each IRLS
#' iteration, it constructs the working response
#' \deqn{z = \eta + (y - \mu) / g'(\eta)}, computes weights
#' \deqn{W = \{g'(\eta)\}^2 / \bigl[\mathrm{Var}(\mu)\times\phi\bigr]}, and
#' then updates \eqn{\beta} by weighted least squares:
#' \deqn{\beta \leftarrow \frac{\sum_i W_i\,x_i\,(z_i - \text{offset}_i)}
#'                       {\sum_i W_i\,x_i^2}.}
#' Iteration stops when the absolute change in \eqn{\beta} is below \code{tol},
#' or when \code{max_iter} is reached.  If a non‐positive or non‐finite
#' denominator arises, a warning is issued and the current estimate is returned.
#'
#' @return
#' A numeric scalar: the estimated slope coefficient \eqn{\beta}.
#'
#' @examples
#' \dontrun{
#' # Gaussian example
#' set.seed(0)
#' x <- rnorm(100)
#' y <- 3*x + rnorm(100)
#' univariate_irls_glm(x, y, family = gaussian(), offset = 0)
#'
#' # Logistic example
#' x <- rnorm(200)
#' p <- 1 / (1 + exp(-x))
#' y <- rbinom(200, 1, p)
#' univariate_irls_glm(x, y, family = binomial(link = "logit"), offset = 0)
#'
#' # Poisson example
#' x <- rnorm(150)
#' lambda <- exp(1 + 2*x)
#' y <- rpois(150, lambda)
#' univariate_irls_glm(x, y, family = poisson(link = "log"), offset = 0)
#' }
#' @export
univariate_irls_glm <- function(x, y,
                                family = gaussian(),
                                offset = 0,
                                dispersion = 1,
                                max_iter = 25,
                                tol = 1e-8) {
  # — Basic sanity checks —
  stopifnot(length(x) == length(y))
  nobs <- length(y)
  if (length(offset) == 1L) offset <- rep(offset, nobs)
  if (!inherits(family, "family")) stop("`family` must be a GLM family object")
  
  weights <- rep(1, nobs)
  etastart <- NULL
  mustart  <- NULL 
  
  eval(family$initialize)
  
  if (!is.null(etastart)) {
    eta <- etastart
  } else {
    eta <- family$linkfun(mustart)
  }
  mu <- family$linkinv(eta)
  
  # — IRLS loop to solve for slope only —
  theta     <- 0
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    var_mu <- family$variance(mu) * rep(dispersion, length(mu))
    gprime <- family$mu.eta(eta)
    
    # working response & weights
    z  <- eta + (y - mu) / gprime
    z0 <- z - offset
    W  <- as.numeric((gprime^2) / var_mu)
    
    # closed‐form update for a single slope β
    num   <- sum(W * x * z0)
    denom <- sum(W * x^2)
    if (denom <= 0 || !is.finite(num/denom)) {
      warning("Non-positive or non-finite denom in IRLS; returning current slope.")
      break
    }
    theta_new <- num / denom
    
    # convergence check
    if (abs(theta_new - theta) < tol) {
      theta     <- theta_new
      converged <- TRUE
      break
    }
    theta <- theta_new
    
    eta <- offset + theta * x
    mu  <- family$linkinv(eta)
  }
  
  return(theta)
}

#' Check Whether True Actives Are Covered by Confidence Sets
#'
#' @description
#' Determines if all true active variables are included in at least one of the
#' provided confidence sets. Returns \code{TRUE} if every element of \code{true_active}
#' appears in the union of the sets, and \code{FALSE} otherwise.
#'
#' @param confidence_sets A list of integer vectors, each representing a confidence set
#'   of selected variable indices. If a single vector is provided, it will be coerced
#'   to a list of length one.
#' @param true_active Integer vector of true active variable indices.
#'
#' @return Logical scalar. \code{TRUE} if all elements of \code{true_active} are
#'   contained in the union of \code{confidence_sets}, \code{FALSE} otherwise.
#'
#' @examples
#' # Single set, true actives 1 and 3 are both covered
#' is_covered(confidence_sets = c(1, 3, 5), true_active = c(1, 3))
#'
#' # Multiple sets, true active 2 appears in one of them
#' sets <- list(cs1 = c(1,4), cs2 = c(2,5), cs3 = c(3))
#' is_covered(confidence_sets = sets, true_active = c(2,3))
#'
#' # Not covered if any true active is missing
#' is_covered(confidence_sets = list(c(1,4)), true_active = c(1,2))
#'
#' @export
is_covered <- function(confidence_sets, true_active) {
  if (!is.list(confidence_sets)) {
    confidence_sets <- list(confidence_sets)
  }
  all_indices <- unique(unlist(confidence_sets, use.names = FALSE))
  all(true_active %in% all_indices)
}

#' Update the Dispersion (Scale) Parameter for GLMs
#'
#' @description
#' Estimate the dispersion (scale) parameter for two-parameter GLM families
#' (Gaussian, Gamma, inverse Gaussian) using either Pearson or deviance residuals.
#' For one-parameter families (binomial, Poisson) and Cox, the dispersion remains 1.
#'
#' @param y Numeric vector of responses.  For Cox models, supply the event times
#'   (first column) if \code{family$family == "cox"}.
#' @param family A \code{stats::family} object (e.g. \code{gaussian()}, \code{Gamma()},
#'   \code{inverse.gaussian()}).  Must include \code{linkinv} and \code{variance} methods,
#'   and—if \code{approach = "deviance"}—a \code{dev.resids} method.
#' @param offset Numeric scalar or vector of linear predictors.  If scalar, recycled
#'   to length \code{length(y)} before applying the inverse link to compute \eqn{\mu}.
#' @param approach Character string, either \code{"pearson"} (default) to use Pearson
#'   residuals, or \code{"deviance"} to use deviance residuals.
#'
#' @return
#' Numeric scalar giving the estimated dispersion:
#' \describe{
#'   \item{Gaussian, Gamma, inverse Gaussian}{Moment‐based estimate}
#'   \item{Binomial, Poisson, Cox}{Always \code{1}}
#' }
#'
#' @details
#' The Pearson estimate is
#' \deqn{\hat{\phi} = \frac{\sum_i (y_i - \mu_i)^2 / V(\mu_i)}{n - p},}
#' and the deviance estimate is
#' \deqn{\hat{\phi} = \frac{\sum_i d_i}{n - p},}
#' where \eqn{d_i} are the deviance residuals returned by
#' \code{family\$}\code{dev.resids}\eqn{(y, \mu, wt)}.  Here, \eqn{p} is the number of estimated
#' parameters (intercept + slopes); for a univariate slope-only model, one may set
#' \eqn{p = 1}, but by default we use \eqn{p = 0} when only updating dispersion.
#'
#' @examples
#' \donttest{
#' # Gaussian with known offset
#' y <- rnorm(100, mean = 2)
#' off <- rep(1, 100)
#' update_dispersion(y, gaussian(), offset = off, approach = "pearson")
#'
#' # Gamma model
#' y <- rgamma(100, shape = 2, scale = 3)
#' off <- rep(0, 100)
#' update_dispersion(y, Gamma(link = "log"), offset = off, approach = "deviance")
#'
#' # Poisson (dispersion fixed at 1)
#' y <- rpois(100, lambda = 5)
#' update_dispersion(y, poisson(), offset = 0)
#' }
#' @export
update_dispersion <- function(y,
                              family = gaussian(),
                              offset = 0,
                              approach = "pearson") {
  # Compute linear predictor's inverse link
  mu <- family$linkinv(offset)
  n  <- length(y)
  p  <- 0  # number of estimated parameters (default 0)
  
  if (approach == "pearson") {
    var_mu    <- family$variance(mu)
    residuals <- (y - mu) / sqrt(var_mu)
    dispersion <- sum(residuals^2) / (n - p)
  } else if (approach == "deviance") {
    dev_resids <- family$dev.resids(y, mu, rep(1, n))
    dispersion <- sum(dev_resids) / (n - p)
  } else {
    stop("`approach` must be either 'pearson' or 'deviance'")
  }
  
  dispersion
}

#' Summarize Confidence Sets Across Simulations
#'
#' @description
#' Given a list of confidence sets (each a vector of selected variable indices)
#' from multiple simulation replicates and the true active set, tabulate
#' the unique sets, their frequencies, proportions, and whether they cover the true set.
#'
#' @param cs_list List of length \code{n_sims}, where each element is an integer
#'   vector (possibly of length zero) of selected variable indices in that simulation.
#' @param true_active Integer vector of the true active variable indices.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{set}{Character representation of each unique confidence set, e.g. "\{1,3\}".}
#'   \item{count}{Number of simulations in which this set was returned.}
#'   \item{percent}{Proportion of simulations with this set (\code{count} / \code{n_sims}).}
#'   \item{cover}{Logical: \code{TRUE} if the set contains all \code{true_active} indices.}
#' }
#' Rows are sorted by decreasing \code{count}.
#'
#' @examples
#' # Suppose over 100 sims we obtained these sets:
#' cs_list <- list(c(1), c(1,2), c(1), integer(0), c(2), c(1))
#' true_active <- 1
#' summarize_cs(cs_list, true_active)
#'
#' @export
summarize_cs <- function(cs_list, true_active) {
  n_sims <- length(cs_list)

  # Helper: normalize any input to a sorted integer vector
  normalize <- function(s) {
    if (is.null(s) || length(s) == 0L) {
      return(integer(0))
    }
    if (is.numeric(s) || is.integer(s)) {
      return(sort(unique(as.integer(s))))
    }
    if (is.list(s)) {
      nums <- unlist(s, recursive = TRUE, use.names = FALSE)
      if (is.numeric(nums)) {
        return(sort(unique(as.integer(nums))))
      }
    }
    warning("Unable to normalize element to integer vector; treating as empty set.")
    integer(0)
  }

  # 1) Normalize each simulation’s set
  norm_list <- lapply(cs_list, normalize)

  # 2) Create string keys for uniqueness
  keys <- vapply(norm_list, function(v) paste(v, collapse = ","), character(1))

  # 3) Identify unique keys and tabulate
  uniq_keys   <- unique(keys)
  factor_keys <- factor(keys, levels = uniq_keys)
  counts      <- as.integer(table(factor_keys))
  percents    <- counts / n_sims

  # 4) Decode unique keys back to integer vectors
  uniq_sets <- lapply(uniq_keys, function(k) {
    if (k == "") integer(0) else as.integer(strsplit(k, ",")[[1]])
  })

  # 5) Compute coverage for each unique set
  covers <- vapply(uniq_sets, is_covered, logical(1), true_active = true_active)

  # 6) Format display strings
  display <- vapply(uniq_sets, function(v) {
    if (length(v) == 0L) "{}" else paste0("{", paste(v, collapse = ","), "}")
  }, character(1))

  # 7) Assemble and sort the data.frame
  df <- data.frame(
    set     = display,
    count   = counts,
    percent = percents,
    cover   = covers,
    stringsAsFactors = FALSE
  )
  df <- df[order(-df$count), ]
  rownames(df) <- NULL
  df
}

#' Summarize Coefficient Estimates Across Simulations
#'
#' @description
#' Given a matrix of simulated coefficient estimates and the true coefficient vector,
#' compute per‐variable summary statistics: the empirical mean and standard deviation
#' of the estimates, alongside the true value.
#'
#' @param sims_coef Numeric matrix of dimension \eqn{p \times n_{\text{sim}}}, where each
#'   row corresponds to one predictor and each column to a simulation replicate.
#' @param true_theta Numeric vector of length \eqn{p}, containing the true coefficient
#'   values for each predictor.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{\code{true}}{The true coefficient values, from \code{true_theta}.}
#'   \item{\code{mean}}{Row‐wise mean of \code{sims_coef}, the average estimated coefficient.}
#'   \item{\code{ssd}}{Row‐wise standard deviation of \code{sims_coef}, the empirical sampling variability.}
#' }
#'
#' @examples
#' # Suppose we ran 100 simulations for 3 predictors
#' set.seed(42)
#' true_theta <- c(1.5, 0, -2)
#' sims_coef  <- matrix(rnorm(3 * 100, mean = rep(true_theta, each = 100), sd = 0.3),
#'                      nrow = 3, byrow = TRUE)
#' summarize_coef(sims_coef, true_theta)
#'
#' @export
summarize_coef <- function(sims_coef, true_theta) {
  if (!is.matrix(sims_coef)) {
    stop("`sims_coef` must be a numeric matrix of size p x n_sims")
  }
  p <- nrow(sims_coef)
  if (length(true_theta) != p) {
    stop("Length of `true_theta` must match the number of rows of `sims_coef`")
  }
  # Compute per-covariate summaries
  mean_est <- rowMeans(sims_coef, na.rm = TRUE)
  ssd_est  <- apply(sims_coef, 1, sd, na.rm = TRUE)
  # Return a data.frame
  data.frame(
    true = true_theta,
    mean = mean_est,
    ssd  = ssd_est,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Construct Confidence Sets from Posterior Model Probabilities
#'
#' @description
#' For each latent effect (column) in a posterior model probability matrix \code{pmp},
#' select the smallest set of variables whose cumulative probability reaches at least
#' \code{coverage}.  Optionally filter out sets whose minimum absolute pairwise
#' correlation (from \code{Rmat}) is below \code{cor_threshold}, and remove duplicates.
#'
#' @param pmp Numeric matrix of dimension \eqn{p \times L}, where each column sums to 1
#'   and contains posterior inclusion probabilities for each variable and effect.
#' @param kept Logical vector of length \eqn{L}, indicating which effects (columns) to process.
#' @param coverage Numeric scalar in \[0,1\], the target cumulative probability for a confidence set.
#'   Defaults to 0.95.
#' @param Rmat Optional \eqn{p \times p} numeric correlation matrix for the predictors.
#'   If supplied, any set whose minimum off-diagonal absolute correlation is below
#'   \code{cor_threshold} is discarded.  Defaults to \code{NULL}.
#' @param cor_threshold Numeric scalar in \[0,1\], the minimum absolute correlation allowed
#'   within a set to pass the correlation filter.  Defaults to 0.5.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{sets}}{A named list of integer vectors.  Each element \code{cs1}, \code{cs2}, …
#'     is a confidence set of variable indices that achieves the target coverage.  If no sets
#'     survive the filters, \code{sets} is \code{NULL}.}
#'   \item{\code{claimed}}{Numeric vector of the actual cumulative probabilities
#'     ("claimed coverage") for each returned set.}
#' }
#'
#' @examples
#' # Simple two-effect example
#' pmp <- matrix(c(0.6, 0.3, 0.1,
#'                 0.2, 0.5, 0.3),
#'               nrow = 3, byrow = FALSE)
#' kept <- c(TRUE, TRUE)
#' cs <- confidence_set(pmp, kept)
#' 
#' # With a correlation filter
#' Rmat <- cor(matrix(rnorm(9), nrow = 3))
#' cs2 <- confidence_set(pmp, kept, coverage = 0.8, Rmat = Rmat, cor_threshold = 0.2)
#'
#' @export
confidence_set <- function(pmp,
                           kept,
                           coverage     = 0.95,
                           Rmat         = NULL,
                           cor_threshold = 0.5) {
  if (!is.matrix(pmp)) stop("`pmp` must be a matrix")
  p <- nrow(pmp); L <- ncol(pmp)
  if (length(kept) != L) stop("`kept` must be length ncol(pmp)")

  sets    <- list()
  claimed <- numeric()

  for (ell in seq_len(L)[kept]) {
    probs <- pmp[, ell]
    ord   <- order(probs, decreasing = TRUE)
    cum   <- cumsum(probs[ord])
    m     <- which(cum >= coverage)[1]
    sel   <- sort(ord[seq_len(m)])
    covel <- sum(probs[sel])

    if (!is.null(Rmat)) {
      if (!all(dim(Rmat) == c(p, p))) {
        stop("`Rmat` must be a p x p matrix where p = nrow(pmp)")
      }
      if (length(sel) > 1) {
        subcorr <- abs(Rmat[sel, sel])
        diag(subcorr) <- NA
        mincorr <- min(subcorr, na.rm = TRUE)
      } else {
        mincorr <- 1
      }
      if (mincorr < cor_threshold) next
    }

    # avoid duplicates
    duplicate <- FALSE
    for (i in seq_along(sets)) {
      if (length(sets[[i]]) == length(sel) && all(sets[[i]] == sel)) {
        duplicate <- TRUE
        break
      }
    }
    if (duplicate) next

    name <- paste0("cs", length(sets) + 1L)
    sets[[name]] <- sel
    claimed[length(sets)] <- covel
  }

  if (length(sets) == 0L) {
    return(list(sets = NULL, claimed = NULL))
  }

  list(sets = sets, claimed = claimed)
}

#' Calculate Cox Proportional Hazards Log-Likelihood
#'
#' @param x Numeric vector of covariate values
#' @param y Matrix with 2 columns (time, status)
#' @param offset Numeric vector or scalar offset
#' @param theta Coefficient value (default: 0)
#' @param ties Method for handling ties ("breslow" or "efron")
#'
#' @return Log-likelihood value
#' @export
univariate_loglik_cox <- function(x, y, offset = numeric(0), theta = 0, ties = "efron") {
  # Type checking
  x <- as.numeric(x)
  y <- as.matrix(y)
  offset <- as.numeric(offset)
  
  # Call C++ function
  .Call("_glmcs_univariate_loglik_cox", PACKAGE = "glmcs", 
        x, y, offset, theta, ties)
}

#' Fit Univariate Cox Regression via IRLS
#'
#' @param x Numeric vector of covariate values
#' @param y Matrix with 2 columns (time, status)
#' @param offset Numeric vector or scalar offset
#' @param ties Method for handling ties ("breslow" or "efron")
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#'
#' @return Estimated coefficient
#' @export
univariate_irls_cox <- function(x, y, offset = numeric(0), ties = "efron", 
                               max_iter = 25, tol = 1e-8) {
  # Type checking
  x <- as.numeric(x)
  y <- as.matrix(y)
  offset <- as.numeric(offset)
  
  # Call C++ function
  .Call("_glmcs_univariate_irls_cox", PACKAGE = "glmcs", 
        x, y, offset, ties, max_iter, tol)
}