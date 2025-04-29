#' Run Example 2 Benchmark Simulation
#'
#' @description
#' Runs a comprehensive benchmark simulation for Example 2 settings, comparing 
#' the performance of multiple variable selection methods. This function is a 
#' convenience wrapper around the \code{\link{benchmark}} function with pre-configured
#' settings for the "ex2" simulation scenario.
#'
#' @details
#' This function configures and executes a simulation benchmark with the following characteristics:
#' 
#' 1. Generates synthetic data with controlled correlation structure (rho), 
#'    intercept, dispersion, and optional censoring for survival data
#' 2. Compares method performance across multiple simulation replicates:
#'    - **glmcs**: Uses the LASER algorithm for variable selection and confidence sets
#'    - **susie**: Uses the Sum of Single Effects model (Gaussian family only)
#' 3. Evaluates methods on:
#'    - Coefficient estimation accuracy
#'    - Confidence set coverage
#'    - Runtime performance
#'
#' The "ex2" setting represents a specific simulation scenario with predefined 
#' sparsity and coefficient structure designed to test confidence set construction
#' in challenging, correlated settings.
#'
#' @param n_sims Integer. Number of simulation replicates to run. Default: 50.
#' @param n Integer. Number of observations in each simulation. Default: 600.
#' @param p Integer. Number of predictors. Default: 5.
#' @param L Integer. Number of components to include in the model. Default: 10.
#' @param family An R family object specifying the response distribution. Default: gaussian().
#' @param rho Numeric in (0,1). Correlation between predictors. Default: 0.9.
#' @param intercept Numeric. True intercept value. Default: -1.
#' @param dispersion Numeric. Dispersion parameter for the model. Default: 9.
#' @param censoring_rate Numeric in [0,1). Censoring rate for survival models. Default: 0.3.
#' @param coverage Numeric in (0,1). Target coverage level for confidence sets. Default: 0.95.
#' @param methods Character vector. Methods to evaluate: "glmcs", "susie", or both. Default: c("glmcs", "susie").
#' @param standardize Logical. Whether to standardize predictors. Default: TRUE.
#' @param ties Character. Method for handling ties in Cox models: "breslow" or "efron". Default: "efron".
#' @param algorithm Character. Algorithm to use: "greedy", "shuffle", or "cyclic" . Default: "greedy".
#' @param max_iter Integer. Maximum number of iterations for optimization. Default: 100.
#' @param step_size Numeric. Step size for the gradient-based optimization. Default: 1.0.
#' @param tol Numeric. Convergence tolerance. Default: 1.0.
#' @param min_abs_corr Numeric in [0,1). Minimum absolute correlation threshold. Default: 0.0.
#' @param parallel Logical. Whether to run simulations in parallel. Default: TRUE.
#' @param cores Integer. Number of cores to use if parallel=TRUE. Default: detectCores() - 1.
#' @param save_path Character or NULL. Path to save intermediate results. Default: NULL.
#' @param save_freq Integer. Frequency of saving intermediate results. Default: 50.
#' @param seed Integer or NULL. Random seed for reproducibility. Default: 42.
#' @param progress Character or logical. Progress tracking method: "bar", "detailed", "text", FALSE. Default: "bar".
#' @param verbose Logical. Whether to display additional information. Default: TRUE.
#'
#' @return A list with benchmark results containing:
#' \describe{
#'   \item{\code{coef_summary}}{List of data frames with coefficient recovery statistics.}
#'   \item{\code{cs_summary}}{List of data frames with confidence set statistics.}
#'   \item{\code{cover_summary}}{List of coverage rates for each method.}
#'   \item{\code{timing_summary}}{Data frame with timing information for each method.}
#'   \item{\code{elapsed_time}}{Total benchmark execution time.}
#'   \item{\code{methods_used}}{Character vector of methods actually used.}
#'   \item{\code{settings}}{List containing all benchmark parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' # Run with default settings (Gaussian response)
#' results <- run_ex2()
#' 
#' # Run with binomial response
#' bin_results <- run_ex2(
#'   n = 500,
#'   p = 10,
#'   family = binomial(),
#'   methods = "glmcs",  # SuSiE not compatible with non-Gaussian
#'   n_sims = 20
#' )
#' 
#' # Run with survival data
#' surv_results <- run_ex2(
#'   family = survival::coxph,
#'   censoring_rate = 0.4,
#'   methods = "glmcs",
#'   n_sims = 30
#' )
#' 
#' # Extract results
#' print(results$cover_summary)
#' print(results$timing_summary)
#' }
#'
#' @seealso \code{\link{benchmark}}, \code{\link{simulate}}
#' @importFrom parallel detectCores
#' @export
run_ex2 <- function(
  # Simulation parameters
  n_sims         = 50L,
  n              = 600L,
  p              = 5L,
  L              = 10L,
  family         = gaussian(),
  rho            = 0.9,
  intercept      = -1,
  dispersion     = 9,
  censoring_rate = 0.3,
  coverage       = 0.95,
  
  # Method selection
  methods        = c("glmcs", "susie"),
  
  # Model fitting parameters
  standardize    = TRUE,
  ties           = "efron",
  algorithm      = "greedy",
  max_iter       = 100L,
  step_size      = 1.0,
  tol            = 1.0,
  min_abs_corr   = 0.0,
  
  # Execution parameters
  parallel       = TRUE,
  cores          = parallel::detectCores() - 1L,
  save_path      = NULL,
  save_freq      = 50L,
  seed           = 42L,
  progress       = "bar",
  verbose        = TRUE
) {  
  # Fixed setting for this function
  settings <- "ex2"
  
  # Start timing
  start_time <- Sys.time()
  
  # Input validation with informative error messages
  validate_inputs <- function() {
    # Validate numeric parameters 
    if (!is.numeric(n) || n <= 0) {
      stop("Parameter 'n' must be a positive integer")
    }
    
    if (!is.numeric(p) || p <= 0) {
      stop("Parameter 'p' must be a positive integer")
    }
    
    if (!is.numeric(L) || L <= 0) {
      stop("Parameter 'L' must be a positive integer")
    }
    
    if (!is.numeric(rho) || rho <= 0 || rho >= 1) {
      warning("Parameter 'rho' must be in (0,1). Using default value 0.9.")
      rho <- 0.9
    }
    
    if (!is.numeric(coverage) || coverage <= 0 || coverage >= 1) {
      warning("Parameter 'coverage' must be in (0,1). Using default value 0.95.")
      coverage <- 0.95
    }
    
    if (!is.numeric(censoring_rate) || censoring_rate < 0 || censoring_rate >= 1) {
      warning("Parameter 'censoring_rate' must be in [0,1). Using default value 0.3.")
      censoring_rate <- 0.3
    }
    
    # Validate categorical parameters
    if (!algorithm %in% c("greedy", "shuffle", "cyclic")) {
      warning("Parameter 'algorithm' must be either 'greedy', 'shuffle', 'cyclic'. Using default 'greedy'.")
      algorithm <- "greedy"
    }
    
    if (!ties %in% c("efron", "breslow")) {
      warning("Parameter 'ties' must be either 'efron' or 'breslow'. Using default 'efron'.")
      ties <- "efron"
    }
    
    # Method compatibility check
    if (is.character(family)) {
      fam_name <- family
    } else if (inherits(family, "family")) {
      fam_name <- family$family
    } else {
      fam_name <- deparse(substitute(family))
    }
    
    # Check SuSiE compatibility
    if ("susie" %in% methods && !identical(fam_name, "gaussian")) {
      warning("SuSiE method only works with Gaussian family. Removing 'susie' from methods.")
      methods <- setdiff(methods, "susie")
    }
    
    # Ensure at least one method remains
    if (length(methods) == 0) {
      stop("No compatible methods specified for the given family.")
    }
    
    # Return any modified parameters
    list(
      rho = rho,
      coverage = coverage,
      censoring_rate = censoring_rate,
      algorithm = algorithm,
      ties = ties,
      methods = methods
    )
  }
  
  # Run validation
  validated <- validate_inputs()
  
  # Update parameters with validated values
  rho <- validated$rho
  coverage <- validated$coverage
  censoring_rate <- validated$censoring_rate
  algorithm <- validated$algorithm
  ties <- validated$ties
  methods <- validated$methods
  
  # Set global seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
    if (verbose) message("Random seed set to: ", seed)
  }
  
  # Assemble simulate() controls
  simulate_control <- list(
    n        = n,
    p        = p,
    family   = family,
    settings = settings,
    control  = list(
      intercept      = intercept,
      dispersion     = dispersion,
      rho            = rho,
      censoring_rate = censoring_rate
    )
  )
  
  # LASER fit controls
  glmcs_control <- list(
    L            = L,
    coverage     = coverage,
    standardize  = standardize,
    ties         = ties,
    algorithm    = algorithm,
    max_iter     = max_iter,
    step_size    = step_size,
    min_abs_corr = min_abs_corr,
    tol          = tol
  )
  
  # SuSiE fit controls (include only if needed)
  susie_control <- if ("susie" %in% methods) {
    list(
      L            = L,
      coverage     = coverage,
      standardize  = standardize,
      min_abs_corr = min_abs_corr,
      max_iter     = max_iter
    )
  } else {
    list()  # Empty list if SuSiE not used
  }
  
  # Configure parallel processing
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' not available. Using sequential processing.")
      parallel <- FALSE
    } else {
      cores <- min(cores, parallel::detectCores())
      if (cores < 1) cores <- 1
      if (verbose) message("Running in parallel with ", cores, " cores.")
    }
  }
  
  # Configure save_path handling
  if (!is.null(save_path)) {
    # Create directory if it doesn't exist
    save_dir <- dirname(save_path)
    if (!dir.exists(save_dir) && save_dir != ".") {
      dir.create(save_dir, recursive = TRUE)
      if (verbose) message("Created directory: ", save_dir)
    }
    
    # Add timestamp to filename if not already specified
    if (!grepl("\\d{4}-\\d{2}-\\d{2}", basename(save_path))) {
      file_ext <- tools::file_ext(save_path)
      file_base <- if (file_ext != "") {
        sub(paste0("\\.", file_ext, "$"), "", save_path)
      } else {
        save_path
      }
      timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
      save_path <- if (file_ext != "") {
        paste0(file_base, "_", timestamp, ".", file_ext)
      } else {
        paste0(file_base, "_", timestamp, ".rds")
      }
      if (verbose) message("Save path with timestamp: ", save_path)
    }
  }
  
  # Call the benchmark function with error handling
  results <- tryCatch({
    benchmark(
      n_sims           = n_sims,
      simulate_control = simulate_control,
      methods          = methods,
      parallel         = parallel,
      cores            = cores,
      glmcs_control    = glmcs_control,
      susie_control    = susie_control,
      save_path        = save_path,
      save_freq        = save_freq,
      progress         = progress
    )
  }, error = function(e) {
    message("Error in benchmark execution: ", e$message)
    
    # Try to provide helpful error diagnostics
    if (grepl("could not find function \"benchmark\"", e$message)) {
      stop("Function 'benchmark' not found. Make sure the package is properly installed and loaded.")
    } else if (grepl("could not find function \"simulate\"", e$message)) {
      stop("Function 'simulate' not found. Make sure the package is properly installed and loaded.")
    } else if (grepl("missing value", e$message)) {
      stop("Error possibly due to NA values in data generation. Check dispersion and correlation parameters.")
    } else {
      stop("Benchmark execution failed: ", e$message)
    }
  })
  
  # Calculate total execution time
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  
  # Add run_ex2 specific information to results
  results$run_ex2_settings <- list(
    n              = n,
    p              = p,
    L              = L,
    family         = family,
    rho            = rho,
    intercept      = intercept,
    dispersion     = dispersion,
    censoring_rate = censoring_rate,
    seed           = seed
  )
  
  if (verbose) {
    message(sprintf("run_ex2 completed in %.2f minutes.", as.numeric(total_time)))
    
    # Print brief summary of results
    methods_used <- results$methods_used
    
    if (length(methods_used) > 0) {
      message("Results summary:")
      for (m in methods_used) {
        cover <- results$cover_summary[m]
        message(sprintf("  - Method '%s': Coverage = %.3f", m, cover))
      }
    }
  }
  
  return(results)
}