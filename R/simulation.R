#’ Run the “Example 1” Benchmark Scenario
#’
#’ @description
#’ A convenience wrapper for the predefined “Example 1” benchmark shipped with **glmcs**.
#’ It constructs sensible defaults for simulation and fitting, sources the packaged
#’ `inst/scripts/benchmark.R` script in an isolated environment, and invokes its
#’ `benchmark()` function.
#’
#’ @param n_sims Integer. Number of simulation replicates. Default: 50.
#’ @param n Integer. Number of observations per simulation. Default: 500.
#’ @param p Integer. Number of predictors. Default: 2.
#’ @param L Integer. Number of single–effect components (LASER). Default: 10.
#’ @param family A GLM family object (e.g. \code{gaussian()}, \code{binomial()}) or
#’   the string \code{"cox"}. Default: \code{gaussian()}.
#’ @param rho Numeric in (0,1). Pairwise correlation for the two–variable setting.
#’   Default: 0.95.
#’ @param intercept Numeric. True intercept for data generation (ignored for Cox). Default: –1.
#’ @param dispersion Numeric. Dispersion for Gaussian/Gamma simulations. Default: 9.
#’ @param censoring_rate Numeric in [0,1). Censoring fraction for Cox simulations. Default: 0.3.
#’ @param coverage Numeric in (0,1). Desired credible–set coverage level. Default: 0.95.
#’ @param methods Character vector. Which methods to run. Supported: \code{"glmcs"}, \code{"susie"}. Default: both.
#’ @param standardize Logical. Should predictors be standardized in LASER fits? Default: \code{TRUE}.
#’ @param ties Character; Cox tie–handling for LASER: \code{"efron"} or \code{"breslow"}. Default: \code{"efron"}.
#’ @param method Character; LASER update order: \code{"greedy"} or \code{"shuffle"}. Default: \code{"greedy"}.
#’ @param max_iter Integer. Maximum coordinate–ascent iterations (LASER). Default: 100.
#’ @param step_size Numeric. Multiplicative step‐size for LASER updates. Default: 1.
#’ @param tol Numeric. Convergence tolerance on log‐likelihood change. Default: 1e-6.
#’ @param parallel Logical. Run simulations in parallel via \code{mclapply}? Default: \code{FALSE}.
#’ @param cores Integer. Number of parallel workers (if \code{parallel=TRUE}). Default: \code{detectCores() - 1}.
#’ @param cors_path Character or \code{NULL}. If a filepath, intermediate results are
#’   saved every \code{save_freq} runs via \code{saveRDS()}. Default: \code{NULL}.
#’ @param save_freq Integer. How often to save intermediate results (in replicates). Default: 50.
#’ @param seed Integer. Seed for reproducibility (affects both simulation and LASER). Default: 42.
#’
#’ @return
#’ The exact output of \code{\link{benchmark}()}: a list with components
#’ \code{coef_summary}, \code{cs_summary}, and \code{raw}.
#’
#’ @details
#’ This helper hides all the boilerplate of setting up control lists,
#’ finding and sourcing the packaged benchmarking script, and scoping
#’ into a temporary environment so that no intermediate objects leak
#’ into your workspace.
#’
#’ @seealso
#’ \itemize{
#’   \item \code{\link{benchmark}}: core benchmarking engine
#’   \item \code{\link{run_ex2}}: equivalent wrapper for “Example 2” (if provided)
#’ }
#’
#’ @examples
#’ \dontrun{
#’ result <- run_ex1(
#’   n_sims         = 100,
#’   n              = 200,
#’   p              = 20,
#’   L              = 5,
#’   family         = gaussian(),
#’   rho            = 0.8,
#’   intercept      = 1,
#’   dispersion     = 2,
#’   censoring_rate = 0.4,
#’   coverage       = 0.9,
#’   methods        = c("glmcs","susie"),
#’   parallel       = TRUE,
#’   cores          = 4,
#’   save_path      = "benchmark_ex1.rds",
#’   save_freq      = 20,
#’   seed           = 42
#’ )
#’ print(result$coef_summary$glmcs)
#’ print(result$cs_summary$susie)
#’ }
#’
#’ @export
run_ex1 <- function(
  n_sims         = 50L,
  n              = 600L,
  p              = 2L,
  L              = 10L,
  family         = gaussian(),
  rho            = 0.95,
  intercept      = -1,
  dispersion     = 9,
  censoring_rate = 0.3,
  coverage       = 0.95,
  methods        = c("glmcs", "susie"),
  standardize    = TRUE,
  ties           = c("efron", "breslow"),
  method         = c("greedy", "shuffle"),
  max_iter       = 100L,
  step_size      = 1.0,
  tol            = 1e-3,
  parallel       = FALSE,
  cores          = parallel::detectCores() - 1L,
  save_path      = NULL,
  save_freq      = 50L,
  seed           = 42L
) {
  # enforce choices
  ties   <- match.arg(ties)
  method <- match.arg(method)
  
  # set global seed
  set.seed(seed)
  
  # assemble simulate() controls
  simulate_control <- list(
    n        = n,
    p        = p,
    family   = family,
    settings = "example-1",
    control  = list(
      intercept      = intercept,
      dispersion     = dispersion,
      rho            = rho,
      censoring_rate = censoring_rate,
      seed           = seed
    )
  )
  
  # LASER fit controls
  fit_control <- list(
    L           = L,
    coverage    = coverage,
    standardize = standardize,
    ties        = ties,
    method      = method,
    max_iter    = max_iter,
    step_size   = step_size,
    tol         = tol,
    seed        = seed
  )
  
  # susie fit controls passed through susie_control arg
  susie_control <- list(max_iter = max_iter)  # or pass in from user
  
  # locate the packaged benchmark script
  script_path <- system.file("scripts", "benchmark.R", package = "glmcs")
  if (script_path == "") {
    stop("Cannot find 'inst/scripts/benchmark.R' in the glmcs package.")
  }
  
  # source in isolated environment
  bench_env <- new.env(parent = baseenv())
  sys.source(script_path, envir = bench_env)
  
  # invoke benchmark()
  bench_env$benchmark(
    n_sims           = n_sims,
    simulate_control = simulate_control,
    methods          = methods,
    parallel         = parallel,
    cores            = cores,
    fit_control      = fit_control,
    susie_control    = susie_control,
    save_path        = save_path,
    save_freq        = save_freq
  )
}