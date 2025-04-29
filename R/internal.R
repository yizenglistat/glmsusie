#' Kullback-Leibler Divergence from Uniform Distribution
#'
#' @description
#' Calculates the Kullback-Leibler (KL) divergence between a given probability
#' distribution `p` and a uniform distribution of equal length. KL divergence
#' measures how much one probability distribution diverges from another
#' expected distribution.
#'
#' In the LASER modeling framework, this function helps quantify how much a
#' posterior distribution deviates from uniformity, which is useful for
#' identifying significant effects.
#'
#' @param p Numeric vector. Probability distribution vector that should sum to 1.
#' @param log_base Numeric. Base of the logarithm used in calculation (default: e = 2.718...).
#'                Common alternatives include 2 (bits) or 10.
#'
#' @return Numeric. KL divergence value (always non-negative). A value of 0 indicates
#'        that `p` is exactly uniform, while larger values indicate greater
#'        divergence from uniformity.
#'
#' @details
#' The KL divergence is calculated as:
#' \deqn{D_{KL}(P||U) = \sum_{i=1}^{n} P(i) \log\left(\frac{P(i)}{U(i)}\right)}
#' where \eqn{U(i) = 1/n} is the uniform distribution.
#'
#' For numerical stability, this implementation:
#' 1. Only processes positive probability elements (zeros contribute nothing)
#' 2. Allows for different logarithm bases
#'
#' The KL divergence is always non-negative and equals zero if and only if
#' the distributions are identical (i.e., `p` is perfectly uniform).
#'
#' @examples
#' # Uniform distribution (KL = 0)
#' p1 <- rep(0.25, 4)
#' kl_divergence(p1)  # Should return 0
#'
#' # Somewhat skewed distribution
#' p2 <- c(0.1, 0.2, 0.3, 0.4)
#' kl_divergence(p2)  # Positive value
#'
#' # Very concentrated distribution
#' p3 <- c(0.01, 0.01, 0.97, 0.01)
#' kl_divergence(p3)  # Large positive value
#'
#' # Using binary logarithm (base 2)
#' kl_divergence(p3, log_base = 2)
#'
#' @references
#' Kullback, S., & Leibler, R. A. (1951). On information and sufficiency.
#' The Annals of Mathematical Statistics, 22(1), 79-86.
#'
#' @seealso \code{\link{get_included}} which uses KL divergence to identify significant effects
#' @name kl_divergence
#' @export
NULL

#' Compute Log-Likelihood for LASER Model
#'
#' @description 
#' Computes the model log-likelihood for various statistical models supported by LASER.
#' Handles both Generalized Linear Models (GLMs) and Cox Proportional Hazards models.
#' For GLMs, supports gaussian, binomial, poisson, gamma, and inverse gaussian families.
#' 
#' @details
#' The linear predictor is computed as: intercept + X * sum(theta, 1), where
#' sum(theta, 1) represents the sum across the latent dimensions (column-wise).
#' For each family, the appropriate link function and log-likelihood formula are applied.
#' 
#' @param X NumericMatrix (n × p) predictor matrix, where:
#'        n = number of observations
#'        p = number of predictors
#' @param y Response variable:
#'        - For GLMs: numeric vector of length n
#'        - For Cox: matrix with 2 columns (time, status), where status = 1 indicates event
#' @param family Model family specification:
#'        - String: "cox", "gaussian", "binomial", "poisson", "Gamma", or "inverse.gaussian"
#'        - List with element "family" or "name" containing one of the above strings
#' @param theta Matrix (p × L) estimated coefficient matrix, where:
#'        p = number of predictors
#'        L = number of latent dimensions
#' @param intercept Numeric scalar. Estimated intercept term in the model.
#' @param dispersion Numeric scalar. Estimated dispersion parameter (only used for certain GLMs).
#'                   Defaults to 1.0.
#' 
#' @return Numeric scalar: total log-likelihood of the model.
#' 
#' @examples
#' \dontrun{
#' # Example for gaussian model with 100 observations, 10 predictors, 2 latent dimensions
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' y <- rnorm(100)
#' theta <- matrix(rnorm(10 * 2), 10, 2)
#' get_loglike(X, y, "gaussian", theta, 0.5, 1.2)
#' }
#' @name get_loglike
#' @export
NULL

#' Calculate Purity Statistics for LASER Model Variables
#'
#' @description 
#' Computes purity statistics (minimum, mean, and median of absolute correlation values)
#' for a set of selected variables. Purity measures help assess variable independence
#' within selected groups.
#' 
#' @details
#' Purity is calculated based on the absolute correlation values between pairs of
#' selected variables. The function can either use a provided correlation matrix or
#' compute correlations from the raw data matrix. For large datasets, a random subset
#' of observations can be used to calculate correlations more efficiently.
#' 
#' @param pos IntegerVector containing indices of selected variables (1-based indexing).
#' @param X Matrix. Original predictor matrix (n x p).
#'           If provided, correlations are calculated from this matrix.
#'           Not required if R is provided.
#' @param R Matrix. Pre-calculated correlation matrix (p x p).
#'               Not required if X is provided.
#' @param squared Logical. If TRUE, squared correlation values are used. Default is FALSE.
#' @param n_purity Integer. Maximum number of observations to use when calculating 
#'                 correlations from X. Default is 100.
#' 
#' @return NumericVector with three elements:
#'         \itemize{
#'           \item Minimum absolute correlation among selected variables
#'           \item Mean absolute correlation among selected variables
#'           \item Median absolute correlation among selected variables
#'         }
#'         If only one variable is selected, returns (1.0, 1.0, 1.0).
#' 
#' @examples
#' \dontrun{
#' # With raw data matrix
#' X <- matrix(rnorm(1000 * 20), 1000, 20)
#' selected_vars <- c(1, 5, 10, 15)
#' get_purity(selected_vars, X = X)
#' 
#' # With correlation matrix
#' Xcorr <- cor(X)
#' get_purity(selected_vars, R = Xcorr, squared = TRUE)
#' }
#' @name get_purity
#' @export
NULL

#' Determine Significant Components in a LASER Model
#'
#' @description
#' Tests which latent components in a LASER model are statistically significant
#' using likelihood ratio tests.
#' 
#' @details
#' For each latent dimension ℓ, the function performs a likelihood ratio test by:
#' 1. Setting the coefficients for dimension ℓ to zero (null model)
#' 2. Computing the log-likelihood of this restricted model
#' 3. Comparing with the full model's log-likelihood using a chi-squared test
#' 4. Determining significance based on the provided alpha level
#' 
#' @param fit List containing LASER model fit objects, including:
#'        - theta: coefficient matrix
#'        - intercept: model intercept
#'        - dispersion: dispersion parameter
#'        - final_loglike: log-likelihood of the full model
#' @param X NumericMatrix (n × p) predictor matrix.
#' @param y Response variable: vector for GLM, matrix (time, status) for Cox.
#' @param family List (for GLM) or string "cox".
#' @param alpha Numeric scalar. Significance level for the tests. Default is 0.05.
#' 
#' @return LogicalVector indicating which latent dimensions are significant.
#' 
#' @examples
#' \dontrun{
#' # Assuming 'laser_fit' is a fitted LASER model with 3 latent dimensions
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' y <- rnorm(100)
#' significant_dims <- get_included(laser_fit, X, y, "gaussian", alpha = 0.01)
#' }
#' 
#' @seealso \code{\link{get_loglike}}
#' @name get_included
#' @export
NULL

#' Compute Confidence Sets for LASER Model Parameters
#'
#' @description
#' Constructs confidence sets for parameters in a LASER model based on posterior probabilities
#' and optionally filters them based on correlation purity metrics.
#' 
#' @details
#' This function constructs confidence sets for each selected latent dimension based on
#' posterior inclusion probabilities. It then:
#' 1. Identifies unique sets to avoid duplication
#' 2. Computes the claimed coverage for each set
#' 3. If correlation data is provided, computes purity metrics for each set
#' 4. Filters sets that don't meet minimum correlation thresholds
#' 5. Returns the filtered sets ordered by purity
#' 
#' Purity measures assess the independence of variables within each confidence set.
#' Higher purity indicates stronger independent signals rather than correlated variables.
#' 
#' @param posterior_mat NumericMatrix (p × L), posterior inclusion probabilities for each variable
#'                      and latent dimension, where:
#'                      p = number of variables
#'                      L = number of latent dimensions
#' @param keep LogicalVector of length L indicating which latent dimensions to include.
#' @param coverage Numeric scalar. Target coverage probability for confidence sets (default: 0.95).
#' @param X NumericMatrix. Original predictor matrix for computing correlations.
#'           Not required if R is provided.
#' @param R NumericMatrix. Pre-computed correlation matrix.
#'               Not required if X is provided.
#' @param check_symmetric Logical. If TRUE, forces the correlation matrix to be symmetric 
#'                        (default: TRUE).
#' @param min_abs_corr Numeric scalar. Minimum absolute correlation threshold for filtering
#'                     confidence sets (default: 0.5).
#' @param n_purity Integer. Maximum number of observations to use for correlation computation
#'                 (default: 100).
#' @param squared Logical. If TRUE, uses squared correlations (default: FALSE).
#' 
#' @return A List containing:
#'         \itemize{
#'           \item sets: List of confidence sets (each as an IntegerVector of variable indices)
#'           \item coverage: NumericVector of claimed coverage for each set
#'         }
#'         If no sets meet the filtering criteria, returns NULL for both elements.
#' 
#' @examples
#' \dontrun{
#' # With posterior matrix and raw data
#' posterior <- matrix(runif(1000 * 3), 1000, 3)
#' keep <- c(TRUE, TRUE, FALSE)
#' X <- matrix(rnorm(100 * 1000), 100, 1000)
#' result <- get_cs(posterior, keep, coverage = 0.9, X = X, min_abs_corr = 0.3)
#' 
#' # With posterior matrix and correlation matrix
#' Xcorr <- cor(X)
#' result <- get_cs(posterior, keep, R = Xcorr, squared = TRUE)
#' }
#' 
#' @seealso \code{\link{get_purity}}
#' @name get_cs
#' @export
NULL

#' Fit a Univariate Generalized Linear Model (GLM) Fast
#'
#' @description
#' Efficiently fit a univariate GLM for predictor \code{x} and response \code{y},
#' supporting Gaussian (closed-form) and Binomial, Poisson, Gamma, Inverse Gaussian (IRLS).
#'
#' @param x Numeric vector. Predictor.
#' @param y Numeric vector. Response.
#' @param family List. Family object with elements \code{name} and optional \code{dispersion}.
#' @param offset Numeric vector. Optional offset.
#' @param standardize Logical. If \code{TRUE}, standardize \code{x} (default \code{TRUE}).
#' @param max_iter Integer. Max IRLS iterations (default 25).
#' @param tol Numeric. Convergence tolerance (default 1e-8).
#'
#' @return A named list: \code{theta}, \code{se}, \code{p_value}, \code{logLik}, \code{BIC}.
#' @name get_glm_fit
#' @export
NULL

#' Fit a Univariate Cox Proportional Hazards Model
#'
#' @description
#' Efficiently fit a univariate Cox model for predictor \code{x} and survival time \code{y}.
#'
#' @param x Numeric vector. Predictor.
#' @param y Numeric matrix (n x 2). First column = time, second = status (1 = event).
#' @param offset Numeric vector. Known offset.
#' @param standardize Logical. If \code{TRUE}, standardize \code{x} (default \code{TRUE}).
#' @param ties Character. Tie method ("efron" or "breslow").
#' @param max_iter Integer. Max Newton-Raphson iterations (default 25).
#' @param tol Numeric. Convergence tolerance (default 1e-8).
#'
#' @return A named list: \code{theta}, \code{se}, \code{p_value}, \code{logLik}, \code{BIC}.
#' @name get_cox_fit
#' @export
NULL

#' Fit a Null (Offset-Only) GLM Model (Fast)
#'
#' @description
#' Fit a generalized linear model with no covariates (only offset), and return
#' the fitted log-likelihood and BIC. Supports Gaussian (closed-form) and
#' Binomial, Poisson, Gamma, Inverse Gaussian via IRLS.
#'
#' @param y Numeric vector. Response values.
#' @param family List. Family object, containing \code{name} (string),
#' and optionally \code{dispersion} (for Gaussian etc.).
#' @param offset Numeric vector. Known offset.
#' @param max_iter Integer. Maximum IRLS iterations (default 25).
#' @param tol Numeric. Convergence tolerance for IRLS (default 1e-8).
#'
#' @return A named list with elements:
#'         \itemize{
#'           \item logLik: Model log-likelihood
#'           \item BIC: Bayesian Information Criterion
#'         }
#' @name null_glm_fit
#' @export
NULL

#' Fit a Null (Offset-Only) Cox Proportional Hazards Model
#'
#' @description
#' Fit a Cox model with only a known offset (no covariates), returning the
#' partial log-likelihood and Bayesian Information Criterion (BIC).
#'
#' @param y NumericMatrix (n × 2). First column: follow-up time; second column: event indicator (1=event, 0=censoring).
#' @param offset Numeric vector. Known offset vector (length 1 or n).
#' @param ties String. Ties handling method: "efron" (default) or "breslow".
#'
#' @return A named list with elements:
#'         \itemize{
#'           \item logLik: Partial log-likelihood
#'           \item BIC: Bayesian Information Criterion
#'         }
#' @name null_cox_fit
#' @export
NULL

#' Update Intercept Given Offset
#'
#' @description
#' Estimate the intercept term in a GLM given a known offset.
#' For Cox models, returns 0.
#' 
#' @param y Response vector or matrix (for Cox).
#' @param family List containing "family" or "name" field.
#' @param offset Known offset vector.
#' 
#' @return Numeric scalar: estimated intercept.
#' @name update_intercept
#' @export
NULL

#' Update Dispersion Parameter
#'
#' @description
#' Given a known offset (linear predictor), update the dispersion parameter.
#' For Gaussian, Gamma, Inverse Gaussian GLM families, computes dispersion from residuals.
#' For binomial, Poisson, or Cox families, returns dispersion = 1.
#'
#' @param y Response vector or matrix (time, status) for Cox.
#' @param family Family description as a List: must contain "family" or "name".
#' @param offset Known offset (linear predictor).
#' 
#' @return Updated dispersion (numeric scalar).
#' @name update_dispersion
#' @export
NULL

#' Fit Single-Effect Regression (SER) Across Predictors
#'
#' @description
#' Performs single-effect regression for each predictor variable against the response,
#' supporting both Generalized Linear Models (GLMs) and Cox Proportional Hazards models.
#'
#' @details
#' For each predictor, this function:
#' 1. Fits a model with just that predictor (and optional offset)
#' 2. Computes coefficient estimate, standard error, and p-value
#' 3. Calculates model fit statistics (log-likelihood, BIC)
#' 4. Derives Bayes factors and posterior probabilities
#'
#' All models are compared against a null model with only the offset.
#'
#' @param X NumericMatrix (n × p). Predictor matrix with:
#'        n = number of observations
#'        p = number of predictors
#' @param y Response variable:
#'        - For GLMs: numeric vector of length n
#'        - For Cox: matrix with 2 columns (time, status)
#' @param family Model family specification:
#'        - String: "cox", "gaussian", "binomial", "poisson", etc.
#'        - List with element "family" or "name" containing family string
#' @param offset NumericVector. Known offset term for the linear predictor:
#'        - Length 1: Same offset applied to all observations
#'        - Length n: Observation-specific offsets
#' @param standardize Logical. Whether to standardize predictors before fitting (default: TRUE).
#' @param ties String. Method for handling ties in Cox regression:
#'        - "efron": Efron approximation (default, more accurate)
#'        - "breslow": Breslow approximation (faster)
#'
#' @return DataFrame with columns:
#'   - theta: Estimated coefficient for each predictor
#'   - se: Standard error of the coefficient estimate
#'   - p_value: Two-sided p-value for coefficient significance
#'   - logLik: Log-likelihood of the model
#'   - BIC: Bayesian Information Criterion
#'   - delta_BIC: Difference in BIC compared to minimum BIC (stabilized)
#'   - bf: Bayes factor derived from delta_BIC
#'   - posterior: Posterior probability for variable importance
#'
#' @examples
#' \dontrun{
#' # Example for Gaussian regression
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' y <- rnorm(100)
#' results <- get_ser_fit(X, y, "gaussian", rep(0, 100))
#' 
#' # Example for Cox regression
#' time <- rexp(100)
#' status <- rbinom(100, 1, 0.5)
#' y_cox <- cbind(time, status)
#' results_cox <- get_ser_fit(X, y_cox, "cox", rep(0, 100), ties = "efron")
#' }
#' @name get_ser_fit
#' @export
NULL

#' Fit LASER Model 
#' 
#' @description
#' Fits the LASER model, which represents regression coefficients as a sum of L single-effects,
#' using blockwise coordinate ascent. The model supports both GLM families and Cox 
#' regression and includes early stopping for efficient computation.
#'
#' @details
#' The LASER model represents the linear predictor as:
#' 
#'    η = intercept + X * Σ(θᵢ)
#' 
#' where θᵢ are L sparse vectors of coefficients. The model is fitted using a block 
#' coordinate ascent algorithm that iteratively updates each effect while holding 
#' others constant. Three update strategies are available:
#' 
#' 1. "cyclic": Update each effect in order 
#' 2. "shuffle": Update effects in random order each iteration (default)
#' 3. "greedy": Select the single effect update giving maximum improvement
#' 
#' For each update, the function uses single-effect regression (SER) to compute posterior
#' probabilities and coefficient estimates. The dispersion parameter and intercept are
#' updated after each complete iteration.
#'
#' @param X NumericMatrix (n × p) predictor matrix, where:
#'        n = number of observations
#'        p = number of predictors
#' @param y Response variable:
#'        - For GLMs: numeric vector of length n
#'        - For Cox: matrix with 2 columns (time, status)
#' @param family Model family specification:
#'        - String: "cox", "gaussian", "binomial", "poisson", etc.
#'        - List with element "family" or "name" containing family string
#' @param L Integer. Number of latent effects to include (default = 10).
#'        If ≤ 0 or > p, automatically set to min(10, p).
#' @param standardize Logical. Whether to standardize predictors (default: TRUE).
#' @param ties String. Method for handling ties in Cox regression (default: "efron").
#' @param algorithm String. Update strategy: "cyclic", "shuffle", or "greedy" (default: "shuffle").
#' @param max_iter Integer. Maximum number of coordinate ascent iterations (default: 100).
#' @param step_size Double. Step size multiplier for updates (default: 1.0).
#'        Values < 1 provide more conservative updates.
#' @param tol Double. Convergence tolerance on log-likelihood change (default: 1).
#' @param seed Integer or NULL. Random seed for reproducibility (default: NULL).
#'        If NULL or NA, no seed is set.
#'
#' @return A list containing:
#'         \itemize{
#'           \item total_elapsed: Total execution time (seconds)
#'           \item elapsed_per_ser: Average time per SER computation
#'           \item n_iter: Number of iterations performed
#'           \item final_loglike: Final model log-likelihood
#'           \item intercept: Estimated intercept term
#'           \item dispersion: Estimated dispersion parameter (for applicable GLMs)
#'           \item theta: Matrix (p×L) of estimated coefficients for each effect
#'           \item posterior: Matrix (p×L) of posterior probabilities
#'           \item p_values: Matrix (p×L) of p-values
#'           \item theta_hat: Vector (length p) of summed coefficients across effects
#'         }
#'
#' @examples
#' \dontrun{
#' # Gaussian example with 10 latent effects
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' laser_fit <- get_laser_fit(X, y, "gaussian", L = 5, algorithm = "cyclic")
#' 
#' # Binomial example with fewer effects
#' y_bin <- rbinom(100, 1, 0.5)
#' laser_bin <- get_laser_fit(X, y_bin, "binomial", L = 3, seed = 123)
#' 
#' # Cox regression example
#' time <- rexp(100)
#' status <- rbinom(100, 1, 0.5)
#' y_cox <- cbind(time, status)
#' laser_cox <- get_laser_fit(X, y_cox, "cox", L = 5, algorithm = "greedy")
#' }
#' 
#' @note
#' The "greedy" algorithm may be slower but often produces higher-quality models,
#' especially for smaller datasets. For large datasets, "cyclic" offers a good
#' balance of speed and quality.
#' 
#' @seealso \code{\link{get_ser_fit}} \code{\link{get_included}} \code{\link{get_cs}}
#' @name get_laser_fit
#' @export
NULL

#' @importFrom stats rexp quantile rbinom rpois rgamma rnorm runif
NULL