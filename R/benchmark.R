#' Benchmark Multiple Variable Selection Methods
#'
#' @description
#' Performs comprehensive benchmarking of variable selection methods including
#' gSuSiE, SuSiE, LASSO, and Elastic Net across multiple simulated datasets.
#' Compares coefficient estimation, variable selection performance, and 
#' credible set coverage.
#'
#' @param settings Data generation settings passed to the \code{generate} function.
#' @param n_sims Number of simulation replicates. Default 1000.
#' @param n Sample size for each simulation. Default 500.
#' @param true_theta Vector of true regression coefficients.
#' @param intercept True intercept value. Default 0.
#' @param family Response family. Can be gaussian(), binomial(), poisson(), 
#'   Gamma(), or "cox" for Cox regression.
#' @param coverage Nominal coverage level for credible sets. Default 0.95.
#' @param L Maximum number of single effects for gSuSiE and SuSiE. Default 10.
#' @param standardize Logical. Should variables be standardized? Default TRUE.
#' @param decompose Logical. Use matrix decomposition in gSuSiE? Default TRUE.
#' @param shrinkage Logical. Apply shrinkage in gSuSiE? Default TRUE.
#' @param tol Convergence tolerance for gSuSiE. Default 5e-2.
#' @param max_iter Maximum iterations for gSuSiE. Default 100.
#' @param ties Method for handling ties in Cox regression. Default "efron".
#' @param lambda Ridge penalty parameter for gSuSiE. Default 0.0.
#' @param tau Threshold parameter for gSuSiE. Default 0.5.
#' @param cor_threshold Correlation threshold for credible sets. Default 0.5.
#' @param include_lasso Logical. Include LASSO in benchmark? Default TRUE.
#' @param include_enet Logical. Include Elastic Net in benchmark? Default TRUE.
#' @param enet_alpha Elastic net mixing parameter. Default 0.5.
#' @param nfolds Number of cross-validation folds for glmnet methods. Default 10.
#' @param seed Random seed for reproducibility. Default 42.
#' @param verbose Logical. Show progress messages? Default TRUE.
#' @param parallel Logical. Use parallel processing? Default FALSE.
#' @param ... Additional arguments passed to the \code{generate} function.
#'
#' @return A list containing benchmark results for each method:
#' \describe{
#'   \item{glmsusie}{List with coef_summary, pmp_summary, and cs_summary}
#'   \item{susie}{List with coef_summary, pip_summary, and cs_summary (Gaussian only)}
#'   \item{lasso}{List with coef_summary, pip_summary, and cs_summary (if included)}
#'   \item{enet}{List with coef_summary, pip_summary, and cs_summary (if included)}
#' }
#' 
#' Each summary contains:
#' \describe{
#'   \item{coef_summary}{Coefficient estimation performance metrics}
#'   \item{pmp_summary/pip_summary}{Inclusion probability performance metrics}
#'   \item{cs_summary}{Credible set coverage and size statistics}
#' }
#'
#' @details
#' This function performs a comprehensive comparison of variable selection methods
#' by running multiple simulations with the same data generation process. For each
#' simulation, it:
#' \enumerate{
#'   \item Generates data using the specified settings
#'   \item Fits all requested methods (gSuSiE, SuSiE, LASSO, Elastic Net)
#'   \item Collects coefficient estimates, inclusion probabilities, and credible sets
#'   \item Summarizes performance across all simulations
#' }
#' 
#' The function automatically handles different response families and excludes
#' SuSiE for non-Gaussian families since it only supports linear regression.
#' 
#' When \code{parallel = TRUE}, simulations are run in parallel using the
#' \code{future} framework. Progress is tracked with a progress bar when
#' \code{verbose = TRUE} and parallel processing is disabled.
#'
#' @examples
#' \dontrun{
#' # Simple benchmark with correlated predictors
#' settings <- list(
#'   rho = 0.8,
#'   block_size = 5
#' )
#' 
#' # Quick benchmark with few simulations
#' results <- benchmark(
#'   settings = settings,
#'   n_sims = 100,
#'   n = 200,
#'   true_theta = c(1, -0.5, rep(0, 18)),
#'   family = gaussian(),
#'   verbose = TRUE
#' )
#' 
#' # Print coefficient estimation results
#' print(results$glmsusie$coef_summary)
#' print(results$lasso$coef_summary)
#' 
#' # Logistic regression benchmark
#' logistic_results <- benchmark(
#'   settings = settings,
#'   n_sims = 50,
#'   n = 500,
#'   true_theta = c(0.5, -0.3, rep(0, 8)),
#'   family = binomial(),
#'   include_lasso = TRUE,
#'   include_enet = FALSE
#' )
#' }
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom susieR susie
#' @export
benchmark <- function(
  settings,
  n_sims       = 1000,
  n            = 500,
  true_theta   = c(1, 0),
  intercept    = 0,
  family       = gaussian(),
  coverage     = 0.95,
  L            = 10,
  standardize  = TRUE,
  decompose    = TRUE,
  shrinkage    = TRUE,
  tol          = 5e-2,
  max_iter     = 100,
  ties         = "efron",
  lambda       = 0.0,
  tau          = 0.5,
  cor_threshold = 0.5,
  include_lasso = TRUE,
  include_enet  = TRUE,
  enet_alpha   = 0.5,
  nfolds       = 10,
  seed         = 42,
  verbose      = FALSE,
  parallel     = FALSE,
  ...
) {
  start_time <- Sys.time()
  set.seed(seed)
  p <- length(true_theta)
  active <- which(true_theta != 0)

  # Prepare parallel backend
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
     stop("Package 'future.apply' is required for parallel processing. Please install it with: install.packages('future.apply')")
    }
    future::plan(future::multisession)
    lapply_fun <- function(X, FUN) future.apply::future_lapply(X, FUN, future.seed = TRUE)
  } else {
    future::plan(future::sequential)
    lapply_fun <- lapply
  }

  # Initialize containers for all methods
  glmsusie_cs    <- vector("list", n_sims)
  glmsusie_theta <- matrix(0, p, n_sims)
  glmsusie_pmp   <- matrix(0, p, n_sims)
  
  susie_cs    <- vector("list", n_sims)
  susie_theta <- matrix(0, p, n_sims)
  susie_pip   <- matrix(0, p, n_sims)
  
  if (include_lasso) {
    lasso_cs    <- vector("list", n_sims)
    lasso_theta <- matrix(0, p, n_sims)
    lasso_pip   <- matrix(0, p, n_sims)
  }
  
  if (include_enet) {
    enet_cs    <- vector("list", n_sims)
    enet_theta <- matrix(0, p, n_sims)
    enet_pip   <- matrix(0, p, n_sims)
  }

  if (!requireNamespace("progress", quietly = TRUE)) {
    message("Install 'progress' for progress bars: install.packages('progress')")
    # Gracefully disable the feature
    verbose <- FALSE
  } else {
    verbose <- TRUE
  }

  # Progress bar
  if (verbose && !parallel) {
    family_name <- if (is.character(family) && family == "cox") "Cox" else family$family
    methods_count <- 2 + include_lasso + include_enet
    pb <- progress::progress_bar$new(
      format = paste0(
        "  Benchmarking ", family_name, " (", methods_count, " methods) [:bar] :percent | Sim: :current/:total | ETA: :eta"
      ),
      total = n_sims, clear = FALSE, width = 90
    )
    pb$tick(0)
  }

  # Simulation
  results <- lapply_fun(seq_len(n_sims), function(i) {
    sim_data <- generate(
      n = n, theta = true_theta, intercept = intercept,
      settings = settings, family = family, ...
    )
    X <- sim_data$X
    y <- sim_data$y

    # Run glmsusie
    res <- glmsusie(
      X = X, y = y, L = L, family = family,
      standardize = standardize, decompose = decompose, shrinkage = shrinkage,
      tol = tol, coverage = coverage, cor_threshold = cor_threshold,
      max_iter = max_iter, ties = ties, lambda = lambda, tau = tau
    )

    # Run susie (only for non-Cox)
    sr <- NULL
    if ((is.character(family) && family != "cox") || (!is.character(family) && family$family != "cox")) {
      sr <- susieR::susie(X, y, L = L, max_iter = 1000)
    }
    
    # Run lasso
    lasso_res <- NULL
    if (include_lasso) {
      lasso_res <- run_lasso(X, y, family = family, standardize = standardize, nfolds = nfolds)
    }
    
    # Run elastic net
    enet_res <- NULL
    if (include_enet) {
      enet_res <- run_elastic_net(X, y, family = family, alpha = enet_alpha, 
                                  standardize = standardize, nfolds = nfolds)
    }

    if (verbose && !parallel) pb$tick()

    list(
      glmsusie_cs    = res$cs$sets,
      glmsusie_theta = rowSums(res$theta),
      glmsusie_pmp   = res$pip,
      susie_cs    = if (!is.null(sr)) sr$sets$cs else NULL,
      susie_theta = if (!is.null(sr)) coef(sr)[-1] else rep(NA, p),
      susie_pip   = if (!is.null(sr)) sr$pip else rep(NA, p),
      lasso_cs    = if (!is.null(lasso_res)) lasso_res$cs$sets else NULL,
      lasso_theta = if (!is.null(lasso_res)) lasso_res$theta else rep(NA, p),
      lasso_pip   = if (!is.null(lasso_res)) lasso_res$pip else rep(NA, p),
      enet_cs     = if (!is.null(enet_res)) enet_res$cs$sets else NULL,
      enet_theta  = if (!is.null(enet_res)) enet_res$theta else rep(NA, p),
      enet_pip    = if (!is.null(enet_res)) enet_res$pip else rep(NA, p)
    )
  })

  # Aggregate results
  for (i in seq_len(n_sims)) {
    glmsusie_cs[[i]]    <- results[[i]]$glmsusie_cs
    glmsusie_theta[, i] <- results[[i]]$glmsusie_theta
    glmsusie_pmp[, i]   <- results[[i]]$glmsusie_pmp
    
    susie_cs[[i]]    <- results[[i]]$susie_cs
    susie_theta[, i] <- results[[i]]$susie_theta
    susie_pip[, i]   <- results[[i]]$susie_pip
    
    if (include_lasso) {
      lasso_cs[[i]]    <- results[[i]]$lasso_cs
      lasso_theta[, i] <- results[[i]]$lasso_theta
      lasso_pip[, i]   <- results[[i]]$lasso_pip
    }
    
    if (include_enet) {
      enet_cs[[i]]    <- results[[i]]$enet_cs
      enet_theta[, i] <- results[[i]]$enet_theta
      enet_pip[, i]   <- results[[i]]$enet_pip
    }
  }

  if (verbose) {
    end_time <- Sys.time()
    elapsed  <- difftime(end_time, start_time, units = "secs")
    message(sprintf("\n Simulation complete in %.2f seconds.", as.numeric(elapsed)))
  }

  # Prepare results
  results_list <- list(
    glmsusie = list(
      coef_summary = summarize_coef(glmsusie_theta, true_theta, 3),
      pmp_summary  = summarize_coef(glmsusie_pmp, true_theta, 3),
      cs_summary   = summarize_cs(glmsusie_cs, active)
    ),
    susie = if (!is.character(family) && family$family == "gaussian") list(
      coef_summary = summarize_coef(susie_theta, true_theta, 3),
      pip_summary  = summarize_coef(susie_pip, true_theta, 3),
      cs_summary   = summarize_cs(susie_cs, active)
    ) else NULL
  )
  
  if (include_lasso) {
    results_list$lasso <- list(
      coef_summary = summarize_coef(lasso_theta, true_theta, 3),
      pip_summary  = summarize_coef(lasso_pip, true_theta, 3),  # Binary PIPs
      cs_summary   = summarize_cs(lasso_cs, active)
    )
  }
  
  if (include_enet) {
    results_list$enet <- list(
      coef_summary = summarize_coef(enet_theta, true_theta, 3),
      pip_summary  = summarize_coef(enet_pip, true_theta, 3),   # Binary PIPs
      cs_summary   = summarize_cs(enet_cs, active)
    )
  }
  
  return(results_list)
}