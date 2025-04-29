#' Benchmark Variable‐Selection Methods
#'
#' @description
#' Run repeated simulations to compare one or more variable‐selection methods
#' (e.g. “glmcs”, “susie”) on their ability to recover true coefficients and
#' confidence sets.
#'
#' @param n_sims Integer. Number of simulation replicates. Default: 50.
#' @param simulate_control List. Arguments forwarded to \code{\link{simulate}()}; 
#'   must include at least \code{n}, \code{p}, \code{family}, \code{settings}, etc.
#' @param methods Character. Which methods to evaluate. Supported: \code{"glmcs"}, 
#'   \code{"susie"}. Default: \code{c("glmcs", "susie")}.
#' @param parallel Logical. Should simulations be run in parallel? Default: \code{FALSE}.
#' @param cores Integer. Number of parallel workers. Default: \code{detectCores() - 1}.
#' @param fit_control List. Arguments forwarded to \code{glmcs()}. Should include  
#'    \code{L, coverage, standardize, ties, method, max_iter, step_size, tol, seed}.
#' @param susie_control List. Arguments forwarded to \code{susieR::susie()}. E.g. \code{list(max_iter=1000)}.
#' @param save_path Character or \code{NULL}. If non-NULL, intermediate results are
#'   saved (via \code{saveRDS()}) every \code{save_freq} iterations.
#' @param save_freq Integer. Save frequency (in replicates). Default: 50.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{coef_summary}}{List of data.frames (one per method) summarizing
#'     the mean and sd of each coefficient over all simulations.}
#'   \item{\code{cs_summary}}{List of data.frames (one per method) listing unique
#'     confidence sets, their frequencies, percentages, and overall coverage.}
#'   \item{\code{raw}}{List of length \code{n_sims}; each element contains the
#'     simulated data plus each method’s output (\code{theta}, \code{cs}, \code{cover}).}
#' }
#'
#' @examples
#' \dontrun{
#' ctrl <- list(
#'   n = 200, p = 20,
#'   family   = gaussian(),
#'   settings = "example-1",
#'   control  = list(intercept=1, dispersion=2, rho=0.8, seed=123)
#' )
#' fit_ctrl   <- list(L=5, coverage=0.9, standardize=TRUE,
#'                    ties="efron", method="greedy",
#'                    max_iter=50, step_size=1, tol=1e-6, seed=NULL)
#' susie_ctrl <- list(max_iter=1000)
#'
#' res <- benchmark(
#'   n_sims           = 100,
#'   simulate_control = ctrl,
#'   methods          = c("glmcs","susie"),
#'   parallel         = FALSE,
#'   fit_control      = fit_ctrl,
#'   susie_control    = susie_ctrl
#' )
#'
#' # Coefficient summaries
#' res$coef_summary$glmcs
#' res$coef_summary$susie
#'
#' # Confidence‐set summaries
#' res$cs_summary$glmcs
#' res$cs_summary$susie
#' }
#'
#' @importFrom parallel detectCores mclapply
#' @export
benchmark <- function(n_sims           = 50,
                      simulate_control = list(),
                      methods          = c("glmcs", "susie"),
                      parallel         = FALSE,
                      cores            = parallel::detectCores() - 1,
                      fit_control      = list(),
                      susie_control    = list(),
                      save_path        = NULL,
                      save_freq        = 50) 
{
  # choose apply function
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
    for (m in methods) {
      if (m == "glmcs") {
        glmcs_fit <- do.call(glmcs, c(
          list(X = sim$X, y = sim$y, family = sim$family),
          fit_control
        ))
        out_i$glmcs <- list(
          theta = glmcs_fit$theta,
          cs    = glmcs_fit$cs$sets,
          cover = is_covered(glmcs_fit$cs$sets, sim$active)
        )

      } else if (m == "susie") {
        susie_fit <- do.call(susieR::susie, c(
          list(X = sim$X, y = sim$y, L = fit_control$L),
          susie_control
        ))
        sets <- susie_fit$sets$cs
        out_i$susie <- list(
          theta = coef(susie_fit)[-1],  # drop intercept
          cs    = sets,
          cover = is_covered(sets, sim$active)
        )

      } else {
        warning("Unknown method: ", m)
      }
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
    # p × n_sims matrix of theta
    theta_list      <- lapply(raw, function(res) res[[m]]$theta)
    theta_mat       <- do.call(cbind, theta_list)
    coef_summary[[m]] <- summarize_coef(theta_mat, true_theta)

    # list of confidence‐sets
    cs_list           <- lapply(raw, function(res) res[[m]]$cs)
    cs_summary[[m]]   <- summarize_cs(cs_list, true_active)
  }

  list(
    coef_summary = coef_summary,
    cs_summary   = cs_summary,
    raw          = raw
  )
}









# library(glmcs)
# library(susieR)

# ctrl <- list(
#   n = 500, p = 2,
#   family   = gaussian(),
#   settings = "example-1",
#   control  = list(intercept=1, dispersion=2, rho=0.8, seed=123)
# )
# fit_ctrl    <- list(L=5, coverage=0.95, standardize=TRUE,
#                     ties="efron", method="greedy",
#                     max_iter=100, step_size=1, tol=1e-6, seed=NULL)
# susie_ctrl  <- list(max_iter=100)

# res <- benchmark(
#   n_sims           = 2,
#   simulate_control = ctrl,
#   methods          = c("glmcs","susie"),
#   parallel         = FALSE,
#   fit_control      = fit_ctrl,
#   susie_control    = susie_ctrl
# )

# # View summaries
# res$coef_summary$glmcs
# res$coef_summary$susie

# res$cs_summary$glmcs
# res$cs_summary$susie