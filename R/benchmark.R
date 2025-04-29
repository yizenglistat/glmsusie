#' Benchmark Variable Selection Methods
#'
#' @description
#' Run comparative simulations to evaluate multiple variable selection methods
#' on their ability to recover true coefficients and construct valid confidence sets.
#' The framework automatically handles method compatibility with different response types.
#'
#' @details
#' This function performs the following steps for each simulation replicate:
#' 
#' 1. Generate synthetic data via \code{\link{simulate}()} with controlled correlation structure
#' 2. Fit each compatible method to the same data:
#'    - **glmcs**: LASER algorithm for all GLM families and Cox models
#'    - **susie**: Sum of Single Effects for Gaussian responses only
#' 3. Collect and evaluate results on:
#'    - Coefficient recovery (bias, variance)
#'    - Confidence set construction (frequency, coverage)
#'    
#' The benchmark supports progress tracking and parallel execution for faster evaluation
#' of large simulation studies. Results are saved periodically if requested.
#'
#' @param n_sims Integer. Number of simulation replicates. Default: 50.
#' @param simulate_control List. Arguments forwarded to \code{\link{simulate}()},
#'   must include at least \code{n}, \code{p}, \code{family}, \code{settings}.
#' @param methods Character vector. Methods to evaluate: \code{"glmcs"}, \code{"susie"}.
#'   Default: c("glmcs", "susie"). Note: SuSiE works with Gaussian family only.
#' @param glmcs_control List. Arguments forwarded to \code{\link{glmcs}()}.
#'   Default parameters include L=10, coverage=0.95, standardize=TRUE, etc.
#' @param susie_control List. Arguments forwarded to \code{susieR::susie()}.
#'   Default parameters include L=10, max_iter=1000, etc.
#' @param parallel Logical. Whether to run simulations in parallel. Default: FALSE.
#' @param cores Integer. Number of parallel workers if parallel=TRUE.
#'   Default: parallel::detectCores() - 1.
#' @param save_path Character or NULL. Path to save intermediate results.
#'   If provided, results are saved every save_freq iterations. Default: NULL.
#' @param save_freq Integer. How often to save intermediate results. Default: 50.
#' @param progress Logical or character. Progress tracking method:
#'   - TRUE or "bar": Progress bar (requires 'progress' package)
#'   - "text": Text-based updates
#'   - FALSE: No progress reporting
#'   Default: TRUE.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{coef_summary}}{List of data frames (one per method) summarizing coefficient recovery.}
#'   \item{\code{cs_summary}}{List of data frames (one per method) summarizing confidence sets.}
#'   \item{\code{cover_summary}}{List of coverage rates (one value per method).}
#'   \item{\code{elapsed_time}}{Total benchmark execution time.}
#'   \item{\code{methods_used}}{Character vector of methods actually used.}
#' }
#'
#' @examples
#' \dontrun{
#' # Gaussian example (both glmcs and susie compatible)
#' sim_ctrl <- list(
#'   n = 200,
#'   p = 20,
#'   family = gaussian(),
#'   settings = "ex1",
#'   control = list(intercept = 1, dispersion = 2, rho = 0.8)
#' )
#' 
#' # Run with default settings
#' results <- benchmark(
#'   n_sims = 50,
#'   simulate_control = sim_ctrl,
#'   methods = c("glmcs", "susie")
#' )
#' 
#' # Extract performance metrics
#' print(results$coef_summary$glmcs)
#' print(results$cover_summary)
#' 
#' # Binomial example (only glmcs compatible)
#' sim_logit <- list(
#'   n = 200,
#'   p = 10,
#'   family = binomial(),
#'   settings = "ex1"
#' )
#' 
#' # Run with parallel processing
#' logit_results <- benchmark(
#'   n_sims = 100,
#'   simulate_control = sim_logit,
#'   methods = "glmcs",
#'   parallel = TRUE,
#'   cores = 4
#' )
#' }
#'
#' @importFrom parallel detectCores mclapply
#' @importFrom utils setTxtProgressBar txtProgressBar packageVersion
#' @export
benchmark <- function(n_sims           = 50,
                      simulate_control = list(),
                      methods          = c("glmcs", "susie"),
                      glmcs_control    = list(),
                      susie_control    = list(),
                      parallel         = FALSE,
                      cores            = parallel::detectCores() - 1,
                      save_path        = NULL,
                      save_freq        = 50,
                      progress         = TRUE) 
{
  # Validate core inputs
  if (!is.list(simulate_control)) {
    stop("simulate_control must be a list")
  }
  
  if (!all(c("n", "p", "family", "settings") %in% names(simulate_control))) {
    stop("simulate_control must include at least n, p, family, and settings")
  }
  
  if (!is.numeric(n_sims) || n_sims < 1) {
    stop("n_sims must be a positive integer")
  }
  n_sims <- as.integer(n_sims)
  
  # Check if family is compatible with requested methods
  family <- simulate_control$family
  fam_name <- if (is.character(family)) family else family$family
  
  # Initialize notes collection
  notes <- character()
  
  # Check method compatibility
  methods_compatible <- methods
  if ("susie" %in% methods && !identical(fam_name, "gaussian")) {
    notes <- c(notes, "Note: SuSiE method only works with Gaussian family and was not used.")
    methods_compatible <- setdiff(methods_compatible, "susie")
  }
  
  # If no compatible methods, stop
  if (length(methods_compatible) == 0) {
    stop("No compatible methods for the specified family.")
  }
  
  # Check required packages
  if ("susie" %in% methods_compatible && !requireNamespace("susieR", quietly = TRUE)) {
    notes <- c(notes, "Note: Package 'susieR' not available. SuSiE method was not used.")
    methods_compatible <- setdiff(methods_compatible, "susie")
  }
  
  # Setup progress tracking
  use_progress_bar <- FALSE
  use_text_progress <- FALSE
  
  if (is.character(progress)) {
    if (progress == "bar") {
      use_progress_bar <- TRUE
    } else if (progress == "text") {
      use_text_progress <- TRUE
    }
  } else if (isTRUE(progress)) {
    # Try to use progress bar, fall back to text
    use_progress_bar <- TRUE
  }
  
  # Check for progress package if needed
  if (use_progress_bar) {
    if (!requireNamespace("progress", quietly = TRUE)) {
      notes <- c(notes, "Note: Package 'progress' not available. Using text progress updates.")
      use_progress_bar <- FALSE
      use_text_progress <- TRUE
    }
  }
  
  # Set up default control parameters if not provided
  if ("glmcs" %in% methods_compatible) {
    default_glmcs <- list(
      L = 10,
      coverage = 0.95,
      standardize = TRUE,
      ties = "efron",
      algorithm = "greedy",
      max_iter = 100,
      step_size = 1.0,
      min_abs_corr = 0.0,
      tol = 1.0
    )
    glmcs_control <- utils::modifyList(default_glmcs, glmcs_control)
  }
  
  if ("susie" %in% methods_compatible) {
    default_susie <- list(
      L = 10,
      max_iter = 1000
    )
    susie_control <- utils::modifyList(default_susie, susie_control)
  }
  
  # Initialize progress tracking
  pb <- NULL
  if (use_progress_bar) {
    pb <- progress::progress_bar$new(
      format = "  Benchmarking [:bar] :percent | Sim :current/:total | Elapsed: :elapsed | ETA: :eta",
      total = n_sims, 
      clear = FALSE,
      width = 80
    )
    pb$tick(0)
  } else if (use_text_progress) {
    message("Running benchmark with methods: ", paste(methods_compatible, collapse=", "))
    if (length(notes) > 0) {
      message(paste(notes, collapse = "\n"))
    }
    message("Starting ", n_sims, " simulations.")
  }

  # Configure parallel processing
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' not available. Using sequential processing.")
      parallel <- FALSE
    } else {
      cores <- min(cores, parallel::detectCores())
      if (cores < 1) cores <- 1
    }
  }
  
  # Choose apply function
  apply_fun <- if (parallel) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    lapply
  }

  # Create container for results
  raw <- vector("list", n_sims)
  
  # Record start time
  start_time <- Sys.time()

  # Main simulation loop
  for (i in seq_len(n_sims)) {
    # Update progress
    if (use_progress_bar) {
      pb$tick()
    } else if (use_text_progress && (i %% 10 == 0 || i == 1 || i == n_sims)) {
      elapsed <- format(difftime(Sys.time(), start_time), digits = 2)
      pct_done <- sprintf("%.1f%%", i / n_sims * 100)
      message(sprintf("Simulation %d of %d (%s) | Elapsed: %s", i, n_sims, pct_done, elapsed))
    }
    
    # 1) Simulate data
    tryCatch({
      sim <- do.call(simulate, simulate_control)
    }, error = function(e) {
      stop("Error generating simulation data (iteration ", i, "): ", e$message)
    })

    # 2) Fit each method
    out_i <- list(simulation = sim)
    
    # Fit glmcs if requested
    if ("glmcs" %in% methods_compatible) {
      tryCatch({
        glmcs_fit <- do.call(glmcs, c(
          list(X = sim$X, y = sim$y, family = sim$family),
          glmcs_control
        ))
        out_i$glmcs <- list(
          theta = glmcs_fit$theta,
          cs    = glmcs_fit$cs$sets,
          cover = is_covered(glmcs_fit$cs$sets, sim$active)
        )
      }, error = function(e) {
        message("Error fitting glmcs in simulation ", i, ": ", e$message)
        out_i$glmcs <- list(
          theta = rep(NA, length(sim$theta)),
          cs    = NULL,
          cover = FALSE
        )
      })
    }
    
    # Fit susie if requested and compatible
    if ("susie" %in% methods_compatible) {
      tryCatch({
        susie_fit <- do.call(susieR::susie, c(
          list(X = sim$X, y = sim$y),
          susie_control
        ))
        sets <- susie_fit$sets$cs
        out_i$susie <- list(
          theta = tryCatch(coef(susie_fit)[-1], error = function(e) rep(NA, length(sim$theta))),
          cs    = sets,
          cover = is_covered(sets, sim$active)
        )
      }, error = function(e) {
        message("Error fitting susie in simulation ", i, ": ", e$message)
        out_i$susie <- list(
          theta = rep(NA, length(sim$theta)),
          cs    = NULL,
          cover = FALSE
        )
      })
    }

    raw[[i]] <- out_i

    # 3) Save intermediate results
    if (!is.null(save_path) && (i %% save_freq == 0L)) {
      saveRDS(raw, save_path)
      if (use_text_progress) message("Intermediate results saved to ", save_path)
    }
  }

  # Calculate total execution time
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  
  if (use_text_progress || use_progress_bar) {
    message(sprintf("Benchmark completed in %.2f minutes.", as.numeric(total_time)))
  }

  # Extract reference values for summarization
  tryCatch({
    true_theta  <- raw[[1]]$simulation$theta
    true_active <- raw[[1]]$simulation$active
  }, error = function(e) {
    stop("Error extracting reference values from simulations: ", e$message)
  })

  # Initialize summary containers
  coef_summary <- setNames(vector("list", length(methods_compatible)), methods_compatible)
  cs_summary   <- setNames(vector("list", length(methods_compatible)), methods_compatible)
  cover_summary <- setNames(vector("numeric", length(methods_compatible)), methods_compatible)

  if (use_text_progress) message("Summarizing results...")
  
  # Summarize results for each method
  for (m in methods_compatible) {
    # Check if method results exist
    method_results <- !all(sapply(raw, function(r) is.null(r[[m]])))
    
    if (method_results) {
      # 1) Coefficient recovery summary
      tryCatch({
        theta_list <- lapply(raw, function(r) r[[m]]$theta)
        theta_mat <- do.call(cbind, theta_list)
        coef_summary[[m]] <- summarize_coef(theta_mat, true_theta)
      }, error = function(e) {
        warning("Error summarizing coefficients for method ", m, ": ", e$message)
        coef_summary[[m]] <- NULL
      })

      # 2) Confidence set summary
      tryCatch({
        cs_list <- lapply(raw, function(r) r[[m]]$cs)
        cs_summary[[m]] <- summarize_cs(cs_list, true_active)
      }, error = function(e) {
        warning("Error summarizing confidence sets for method ", m, ": ", e$message)
        cs_summary[[m]] <- NULL
      })

      # 3) Overall coverage rate
      tryCatch({
        cover_vals <- sapply(raw, function(r) r[[m]]$cover)
        cover_summary[m] <- mean(cover_vals, na.rm = TRUE)
      }, error = function(e) {
        warning("Error calculating coverage rate for method ", m, ": ", e$message)
        cover_summary[m] <- NA
      })
    } else {
      message("No results available for method: ", m)
      coef_summary[[m]] <- NULL
      cs_summary[[m]] <- NULL
      cover_summary[m] <- NA
    }
  }

  if (use_text_progress) message("Benchmark completed successfully!")
  
  # Return final results
  list(
    coef_summary  = coef_summary,
    cs_summary    = cs_summary,
    cover_summary = cover_summary,
    elapsed_time  = total_time,
    methods_used  = methods_compatible,
    notes         = notes
  )
}
#run_ex1(n_sims=1000,n=100000, family="cox",intercept=0, algorithm='greedy',tol=1, max_iter=1000, seed=NULL)