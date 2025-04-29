#’ Benchmark Variable‐Selection Methods
#’
#’ @description
#’ Run repeated simulations to compare one or more variable‐selection methods
#’ (currently “glmcs” and “susie”) on their ability to recover true coefficients
#’ and confidence sets.
#’
#’ @details
#’ For each replicate, \code{\link{simulate}()} is called with \code{simulate_control},
#’ then each specified method is fit:
#’ - **glmcs** via \code{\link{glmcs}()}, controlled by \code{glmcs_control}.
#’ - **susie** via \code{susieR::susie}(), controlled by \code{susie_control}.
#’ 
#’ After all replicates, the function summarizes:
#’ 1. Coefficient recovery: mean ± sd over all simulated coefficients.
#’ 2. Confidence‐set behavior: unique sets, their frequencies (%), and coverage.
#’
#’ @param n_sims Integer. Number of simulation replicates. Default: 50.
#’ @param simulate_control List. Arguments forwarded to \code{\link{simulate}()};
#’   must include at least \code{n}, \code{p}, \code{family}, \code{settings}, etc.
#’ @param methods Character vector. Which methods to evaluate. Supported:
#’   \code{"glmcs"}, \code{"susie"}. Default: \code{c("glmcs", "susie")}.
#’ @param glmcs_control List. Arguments forwarded to \code{\link{glmcs}()}; e.g.
#’   \code{list(L=10, coverage=0.95, standardize=TRUE, ties="efron", method="greedy",
#’   max_iter=100, step_size=1, tol=1e-6, seed=NULL)}.
#’ @param susie_control List. Arguments forwarded to \code{susieR::susie()}; e.g.
#’   \code{list(L=10, max_iter=1000)}.
#’ @param parallel Logical. Should simulations be run in parallel? Default: \code{FALSE}.
#’ @param cores Integer. Number of parallel workers (only if \code{parallel = TRUE}).
#’   Default: \code{parallel::detectCores() - 1}.
#’ @param save_path Character or \code{NULL}. If non‐NULL, intermediate \code{raw}
#’   results are saved via \code{saveRDS()} every \code{save_freq} replicates.
#’ @param save_freq Integer. How often (in replicates) to save the intermediate results.
#’   Default: 50.
#’
#’ @return A named list with components:
#’ \describe{
#’   \item{\code{coef_summary}}{Named list of data.frames (one per method) with
#’     mean and standard deviation of each coefficient over all simulations.}
#’   \item{\code{cs_summary}}{Named list of data.frames (one per method) listing
#’     unique confidence sets, their counts, percentages, and coverage status.}
#’   \item{\code{raw}}{List (length = \code{n_sims}) of per‐replicate results,
#’     each containing the simulated data and each method’s output
#’     (\code{theta}, \code{cs}, \code{cover}).}
#’ }
#’
#’ @examples
#’ \dontrun{
#’ # Define simulation settings
#’ sim_ctrl <- list(
#’   n        = 200,
#’   p        = 20,
#’   family   = gaussian(),
#’   settings = "ex1",
#’   control  = list(intercept=1, dispersion=2, rho=0.8, seed=123)
#’ )
#’
#’ # Fit controls for each method
#’ glmcs_ctrl   <- list(L=5, coverage=0.9, standardize=TRUE,
#’                      ties="efron", method="greedy",
#’                      max_iter=50, step_size=1, tol=1e-6, seed=NULL)
#’ susie_ctrl <- list(L=5, max_iter=1000)
#’
#’ # Run benchmark
#’ res <- benchmark(
#’   n_sims           = 100,
#’   simulate_control = sim_ctrl,
#’   methods          = c("glmcs", "susie"),
#’   glmcs_control    = glmcs_ctrl,
#’   susie_control    = susie_ctrl,
#’   parallel         = FALSE
#’ )
#’
#’ # Inspect summaries
#’ print(res$coef_summary$glmcs)
#’ print(res$cs_summary$susie)
#’ }
#’
#’ @importFrom parallel detectCores mclapply
#’ @export
benchmark <- function(n_sims           = 50,
                      simulate_control = list(),
                      methods          = c("glmcs", "susie"),
                      glmcs_control    = list(),
                      susie_control    = list(),
                      parallel         = FALSE,
                      cores            = parallel::detectCores() - 1,
                      save_path        = NULL,
                      save_freq        = 50) 
{
  ## choose apply function
  apply_fun <- if (parallel) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    lapply
  }

  raw <- vector("list", n_sims)

  for (i in seq_len(n_sims)) {
    ## 1) simulate data
    sim <- do.call(simulate, simulate_control)

    ## 2) fit each method
    out_i <- list(simulation = sim)
    if ("glmcs" %in% methods) {
      glmcs_fit <- do.call(glmcs, c(
        list(X = sim$X, y = sim$y, family = sim$family),
        glmcs_control
      ))
      out_i$glmcs <- list(
        theta = glmcs_fit$theta,
        cs    = glmcs_fit$cs$sets,
        cover = is_covered(glmcs_fit$cs$sets, sim$active)
      )
    }
    if ("susie" %in% methods) {
      susie_fit <- do.call(susieR::susie, c(
        list(X = sim$X, y = sim$y),
        susie_control
      ))
      sets <- susie_fit$sets$cs
      out_i$susie <- list(
        theta = coef(susie_fit)[-1],  # drop intercept
        cs    = sets,
        cover = is_covered(sets, sim$active)
      )
    }

    raw[[i]] <- out_i

    ## 3) save intermittently
    if (!is.null(save_path) && (i %% save_freq == 0L)) {
      saveRDS(raw, save_path)
    }
  }

  ## 4) summarize across simulations
  true_theta  <- raw[[1]]$simulation$theta
  true_active <- raw[[1]]$simulation$active

  coef_summary <- setNames(vector("list", length(methods)), methods)
  cs_summary   <- setNames(vector("list", length(methods)), methods)

  for (m in methods) {
    # build p × n_sims matrix of theta
    theta_list       <- lapply(raw, `[[`, m) %>% 
                          lapply(`[[`, "theta")
    theta_mat        <- do.call(cbind, theta_list)
    coef_summary[[m]] <- summarize_coef(theta_mat, true_theta)

    # build list of confidence‐sets
    cs_list          <- lapply(raw, `[[`, m) %>% 
                          lapply(`[[`, "cs")
    cs_summary[[m]]   <- summarize_cs(cs_list, true_active)
  }

  list(
    coef_summary = coef_summary,
    cs_summary   = cs_summary,
    raw          = raw
  )
}