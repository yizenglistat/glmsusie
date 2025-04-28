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