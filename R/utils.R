#’ Simulate Data for GLM and Cox Models
#’
#’ @description 
#’ Generate synthetic datasets with controlled correlation structure and known truth
#’ for benchmarking variable‐selection methods (e.g. \code{glmcs}, \code{susie}, \code{glmnet}, …).
#’
#’ @param n Integer. Number of observations (default: 500).
#’ @param p Integer. Number of predictors (default: 10).
#’ @param family Family object (e.g. \code{gaussian()}, \code{binomial()}, or \code{"cox"}).
#’ @param settings Character. Which correlation scenario to use. Supported:
#’   \describe{
#’     \item{\code{"example-1"}}{2 highly correlated variables, only the first active.}
#’     \item{\code{"example-2"}}{A 5×5 block‐correlation among the first 5 predictors, two of which are active.}
#’   }
#’ @param control Named list of simulation settings (all optional):
#’   \describe{
#’     \item{\code{intercept}}{Numeric. True intercept (default: 0). Ignored for Cox.}
#’     \item{\code{dispersion}}{Numeric. Dispersion for Gaussian/Gamma (default: 1).}
#’     \item{\code{rho}}{Numeric. Pairwise correlation for “example-1” (default: 0.9).}
#’     \item{\code{censoring_rate}}{Numeric in [0,1). Censoring fraction for Cox (default: 0.3).}
#’     \item{\code{seed}}{Integer. RNG seed (default: \code{NULL}).}
#’   }
#’
#’ @return A list with components:
#’   \item{X}{\eqn{n\times p} design matrix.}
#’   \item{y}{Response: numeric vector for GLMs, or \eqn{n\times2} \code{(time,status)} for Cox.}
#’   \item{theta}{True coefficient vector (length \code{p}).}
#’   \item{active}{Indices of the nonzero entries in \code{theta}.}
#’   \item{intercept, dispersion, family, settings, control}{Echoed inputs.}
#’
#’ @examples
#’ sim1 <- simulate( n=200, p=10, family=gaussian(), settings="example-1",
#’                        control=list(intercept=1.5, dispersion=2, rho=0.8, seed=123) )
#’ sim2 <- simulate( n=300, p=8,  family=binomial(),  settings="example-2",
#’                        control=list(seed=42) )
#’ sim3 <- simulate( n=150, p=5,  family="cox",        settings="example-1",
#’                        control=list(censoring_rate=0.5, seed=99) )
#’
#’ @importFrom MASS mvrnorm
#’ @importFrom Matrix nearPD
#' @importFrom utils modifyList
#’ @importFrom stats rnorm rpois rbinom rexp runif
#’ @export
simulate <- function(n = 500,
                  p = 10,
                  family = gaussian(),
                  settings = "example-1",
                  control = list()) {

  ## 1) merge defaults
  defaults <- list(
    intercept      = 0,      # true intercept
    dispersion     = 1,      # for Gaussian/Gamma
    rho            = 0.9,    # correlation for example-1
    censoring_rate = 0.3,    # for Cox
    seed           = NULL    # RNG
  )
  ctrl <- utils::modifyList(defaults, control)

  ## 2) reproducibility
  if (!is.null(ctrl$seed)) set.seed(ctrl$seed)

  ## 3) build two scenarios
  X <- matrix(0, nrow = n, ncol = p)
  theta <- numeric(p)

  if (settings == "example-1") {
    if (p < 2) stop("example-1 requires p >= 2")
    Sigma2 <- matrix(c(1, ctrl$rho, ctrl$rho, 1), 2, 2)
    X12 <- MASS::mvrnorm(n, mu = c(0,0), Sigma = Sigma2)
    X[,1:2] <- X12
    if (p>2) X[,3:p] <- matrix(rnorm(n*(p-2)), n, p-2)
    theta[1] <- 1
  } else if (settings == "example-2") {
    if (p < 5) stop("example-2 requires p >= 5")
    S5 <- matrix(c(
      1.0, 0.92, 0.7,  0.7,  0.9,
      0.92,1.0,  0.7,  0.7,  0.7,
      0.7, 0.7,  1.0,  0.92, 0.8,
      0.7, 0.7,  0.92, 1.0,  0.8,
      0.9, 0.7,  0.8,  0.8,  1.0), 5,5, byrow=TRUE)
    S5 <- as.matrix(Matrix::nearPD(S5)$mat)
    X5 <- MASS::mvrnorm(n, mu=rep(0,5), Sigma=S5)
    X[,1:5] <- X5
    if (p>5) X[,6:p] <- matrix(rnorm(n*(p-5)), n, p-5)
    theta[c(2,3)] <- 1
  } else {
    stop("Unknown settings: ", settings)
  }

  ## 4) linear predictor & response
  eta <- ctrl$intercept + X %*% theta
  fam_name <- if (is.character(family)) family else family$family

  if (fam_name=="cox") {
    # survival times
    rate <- exp(eta)
    T_event <- rexp(n, rate=rate)
    Censor  <- runif(n,0, quantile(T_event,1-ctrl$censoring_rate))
    time    <- pmin(T_event, Censor)
    status  <- as.integer(T_event <= Censor)
    y <- cbind(time=time, status=status)
  } else {
    inv_link <- if (is.character(family)) {
      stop("non-cox character families not supported")
    } else family$linkinv
    mu <- inv_link(eta)
    if (fam_name=="gaussian") {
      y <- as.numeric(mu + rnorm(n, sd=sqrt(ctrl$dispersion)))
    } else if (fam_name=="binomial") {
      y <- rbinom(n,1,mu)
    } else if (fam_name=="poisson") {
      y <- rpois(n, mu)
    } else if (fam_name=="Gamma") {
      shape <- 1/ctrl$dispersion
      scale <- mu * ctrl$dispersion
      y <- rgamma(n, shape=shape, scale=scale)
    } else {
      stop("Unsupported family: ", fam_name)
    }
  }

  list(
    X           = X,
    y           = y,
    theta       = theta,
    active      = which(theta!=0),
    intercept   = ctrl$intercept,
    dispersion  = ctrl$dispersion,
    family      = family,
    settings    = settings,
    control     = ctrl
  )
}

#' Modify GLM Family Objects or Create Cox Family with Fixed Dispersion (if any)
#'
#' @title Modify GLM Family Objects or Create Cox Family with Fixed Dispersion
#' 
#' @description
#' Creates or modifies family objects for use in generalized linear models (GLMs) 
#' or Cox proportional hazards models. For two-parameter GLMs (gaussian, Gamma, 
#' inverse.gaussian), this function sets a fixed dispersion parameter, modifying
#' the variance, deviance, and AIC functions accordingly.
#'
#' @details
#' For GLM families, this function:
#' 
#' 1. Validates if the family supports a dispersion parameter
#' 2. For two-parameter families (gaussian, Gamma, inverse.gaussian), modifies:
#'    - `variance()`: Returns dispersion-adjusted variance function
#'    - `dev.resids()`: Returns dispersion-adjusted deviance residuals
#'    - `aic()`: Returns dispersion-adjusted AIC
#' 3. For one-parameter families (binomial, poisson), ignores the dispersion
#'    parameter and sets it to 1
#' 4. Optionally validates the family by testing link functions and other components
#' 
#' For "cox", returns the string "cox" directly, which is recognized by the LASER
#' modeling functions as a special case.
#'
#' @param family A GLM family object (e.g., gaussian(), binomial(), poisson(), 
#'               Gamma(), inverse.gaussian()) or the string "cox".
#' @param dispersion Positive numeric scalar. The fixed dispersion parameter value
#'                   for two-parameter GLM families. Ignored for one-parameter
#'                   families (binomial, poisson) and "cox". Default is 1.
#' @param validate Logical. If TRUE, performs basic consistency checks on the 
#'                 modified family object (link function reversibility, etc.)
#'                 Default is TRUE.
#'
#' @return 
#' For GLM families, returns a modified family object (list) with components:
#'   * `family`: Character string naming the family
#'   * `link`: Character string naming the link function
#'   * `linkfun()`: Link function
#'   * `linkinv()`: Inverse link function
#'   * `variance()`: Variance function (modified for two-parameter families)
#'   * `dev.resids()`: Deviance residuals function (modified for two-parameter families)
#'   * `aic()`: AIC function (modified for two-parameter families)
#'   * `dispersion`: Numeric scalar with the fixed dispersion value
#'   * ... and other components from the original family object
#'
#' For "cox", returns the string "cox".
#'
#' @note
#' The dispersion parameter in R is the inverse of the "shape" parameter 
#' often used in statistical literature, particularly for the Gamma distribution.
#'
#' @examples
#' # Standard Gaussian family with default dispersion
#' gaussian_fam <- get_family(gaussian())
#' 
#' # Gaussian family with custom dispersion
#' gaussian_fam2 <- get_family(gaussian(), dispersion = 2.5)
#' 
#' # Cox proportional hazards
#' cox_fam <- get_family("cox")
#' 
#' # Gamma distribution with dispersion = 0.5 (shape = 2)
#' gamma_fam <- get_family(Gamma(), dispersion = 0.5)
#' 
#' # One-parameter family (ignores dispersion)
#' poisson_fam <- get_family(poisson(), dispersion = 2) # dispersion will be set to 1
#'
#' @seealso \code{\link[stats]{family}}
#' 
#' @export
get_family <- function(family,
                       dispersion = 1,
                       validate = TRUE) {
  # --- Handle string case: must be "cox" ---
  if (is.character(family)) {
    if (!identical(family, "cox")) {
      stop("If 'family' is a string, it must be exactly 'cox'.")
    }
    # Return minimal cox-family object
    return(family)
  }
  
  # --- Validate family is a GLM family object ---
  if (!inherits(family, "family")) {
    stop("`family` must be a GLM family object or the string 'cox'.")
  }
  
  # --- Extract family name and determine if it supports dispersion ---
  fam_name <- family$family
  two_param <- fam_name %in% c("gaussian", "Gamma", "inverse.gaussian")
  
  # --- Validate dispersion parameter ---
  if (!is.numeric(dispersion) || length(dispersion) != 1 || dispersion <= 0) {
    warning("`dispersion` must be a positive scalar; resetting to 1.")
    dispersion <- 1
  }
  
  # --- Modify two-parameter families to use fixed dispersion ---
  if (two_param) {
    if (fam_name == "gaussian") {
      # Gaussian (normal) distribution
      family$variance <- function(mu) rep(dispersion, length(mu))
      family$dev.resids <- function(y, mu, wt) wt * (y - mu)^2 / dispersion
      family$aic <- function(y, n, mu, wt, dev)
        sum(wt) * (log(2 * pi * dispersion) + 1) + dev / dispersion
    } else if (fam_name == "Gamma") {
      # Gamma distribution
      shape <- 1 / dispersion
      family$variance <- function(mu) mu^2 / shape
      family$dev.resids <- function(y, mu, wt)
        -2 * wt * shape * (log(y / mu) - (y - mu) / mu)
      family$aic <- function(y, n, mu, wt, dev)
        -2 * sum(wt * shape * log(shape)) + 2 * sum(wt * shape) +
          2 * sum(wt * lgamma(shape)) + dev / dispersion
    } else if (fam_name == "inverse.gaussian") {
      # Inverse Gaussian distribution
      family$variance <- function(mu) mu^3 / dispersion
      family$dev.resids <- function(y, mu, wt)
        wt * (y - mu)^2 / (mu^2 * y * dispersion)
      family$aic <- function(y, n, mu, wt, dev)
        sum(wt) * (log(2 * pi * dispersion) - log(y)) + dev / dispersion
    }
    # Store dispersion parameter in family object
    family$dispersion <- dispersion
  } else {
    # One-parameter GLMs (binomial, poisson) ignore user dispersion
    family$dispersion <- 1
  }
  
  # --- Optional validation checks ---
  if (validate && fam_name != "cox") {
    try({
      # Create test values appropriate for each family
      test_mu <- switch(fam_name,
                         gaussian = c(-1, 0, 1),
                         binomial = c(0.2, 0.5, 0.8),
                         poisson = c(1, 2, 3),
                         Gamma = c(1, 2, 3),
                         inverse.gaussian = c(1, 2, 3))
      test_y <- test_mu  # Use means as simulated responses for testing
      test_wt <- rep(1, length(test_mu))
      
      # Test variance and deviance functions
      family$variance(test_mu)
      family$dev.resids(test_y, test_mu, test_wt)
      
      # Test link function consistency for continuous families
      if (!fam_name %in% c("binomial", "poisson")) {
        eta <- family$linkfun(test_mu)
        mu_check <- family$linkinv(eta)
        if (max(abs(test_mu - mu_check)) > 1e-8) {
          warning("Link function inconsistency detected.")
        }
      }
    }, silent = TRUE)
  }
  
  # Return the modified family object
  return(family)
}

#’ Check Confidence Set Coverage
#’
#’ Return `TRUE` if all truly active variables are contained in the confidence sets,
#’ otherwise `FALSE`.
#’
#’ @param confidence_sets A list of integer vectors (each a confidence set),  
#’   or a single integer vector.
#’ @param true_active Integer vector of truly active covariate indices.
#’
#’ @return Logical `TRUE` if every element of \code{true_active} appears in at least one set,  
#’   otherwise `FALSE`.
#’
#’ @examples
#’ # multiple confidence sets
#’ sets   <- list(c(3,5,7), c(1,4,7), c(2,6))
#’ true   <- c(1,2,7)
#’ is_covered(confidence_sets = sets, true_active = true)  # TRUE
#’
#’ # single set (e.g. lasso support)
#’ lasso_support <- c(2,7,9)
#’ is_covered(confidence_sets = lasso_support, true_active = true)  # FALSE
#’
#’ @export
is_covered <- function(confidence_sets, true_active) {
  if (!is.list(confidence_sets)) {
    confidence_sets <- list(confidence_sets)
  }
  all_indices <- unique(unlist(confidence_sets, use.names = FALSE))
  all(true_active %in% all_indices)
}

#' Summarize Confidence Sets Across Simulations
#'
#' @title Summarize Confidence Sets Across Simulations
#'
#' @description
#' Aggregate and tabulate the distinct confidence sets produced over multiple
#' simulation replicates, computing their frequencies, relative proportions,
#' and whether they successfully cover all true active predictors.
#'
#' @details
#' Given a list of confidence-set vectors (one per simulation) and the known
#' true active indices, this function:
#' \enumerate{
#'   \item Normalizes each set (sorts and deduplicates, treating \code{NULL}
#'         or zero-length as the empty set \{\});
#'   \item Identifies unique configurations in order of first appearance;
#'   \item Counts how often each configuration occurs and computes its percentage
#'         of the total simulations;
#'   \item Determines whether each unique set covers all true actives by calling
#'         \code{is_covered()};
#'   \item Returns a data.frame sorted by descending frequency, with no row names.
#' }
#'
#' @param cs_list List of length \code{n_sims}.  Each element must be either:
#'   \itemize{
#'     \item An integer vector of selected covariate indices, or
#'     \item \code{NULL} / zero-length, interpreted as the empty set \{\}.
#'   }
#' @param true_active Integer vector. The indices of truly active covariates.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{\code{set}}{Character string representation of the set, e.g.\ \code{"\{1,2\}"} or \code{"\{\}"}.}
#'   \item{\code{count}}{Integer. Number of simulations in which this exact set occurred.}
#'   \item{\code{percent}}{Numeric. \code{count} divided by \code{length(cs_list)}.}
#'   \item{\code{cover}}{Logical. \code{TRUE} if \code{is_covered(set, true_active)}; \code{FALSE} otherwise.}
#' }
#'
#' @seealso \code{is_covered} for checking whether a single confidence set
#' contains all true active indices.
#'
#' @examples
#' \dontrun{
#' # Simulated list of sets (NULL or integer(0) indicate empty set)
#' cs_list    <- list(c(1,2), NULL, c(3), c(2,1), integer(0))
#' true_active <- c(1,2)
#'
#' # Summarize their frequencies and coverage
#' summarize_cs(cs_list, true_active)
#' }
#'
#' @export
summarize_cs <- function(cs_list, true_active) {
  n_sims <- length(cs_list)

  # 1) Normalize each to sorted integer vector (NULL → integer(0))
  norm_list <- lapply(cs_list, function(s) {
    if (is.null(s) || length(s) == 0L) integer(0)
    else sort(unique(as.integer(s)))
  })

  # 2) Create keys for uniqueness
  keys      <- vapply(norm_list, paste, "", collapse = ",")
  uniq_keys <- unique(keys)
  factor_keys <- factor(keys, levels = uniq_keys)

  # 3) Tabulate counts and percentages
  counts   <- as.integer(table(factor_keys))
  percents <- counts / n_sims

  # 4) Decode unique key strings back to integer vectors
  uniq_sets <- lapply(uniq_keys, function(k) {
    if (k == "") integer(0)
    else as.integer(strsplit(k, ",")[[1]])
  })

  # 5) Evaluate coverage for each unique set
  covers <- vapply(uniq_sets, is_covered, logical(1), true_active = true_active)

  # 6) Format display strings for each set
  disp_sets <- vapply(uniq_sets, function(v) {
    if (length(v) == 0L) return("{}")
    paste0("{", paste(v, collapse = ","), "}")
  }, character(1))

  # 7) Assemble and order result
  df <- data.frame(
    set     = disp_sets,
    count   = counts,
    percent = percents,
    cover   = covers,
    stringsAsFactors = FALSE
  )
  df <- df[order(-df$count), , drop = FALSE]
  rownames(df) <- NULL
  df
}








