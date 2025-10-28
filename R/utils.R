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

#' Summarize Multi‐Set Confidence Regions Across Simulations
#'
#' @description
#' Given a list of simulation results where each run may return multiple
#' confidence sets (e.g. \code{cs1}, \code{cs2}, …), this function enumerates
#' the distinct patterns of sets, counts how often each occurs, computes the
#' fraction of simulations, and tests whether the **union** of all sets in
#' a pattern covers the true active variables.  Within each pattern, the
#' individual sets are ordered first by their size (smallest first), then
#' lexicographically by their indices.
#'
#' @param cs_list List of length \code{n_sims}. Each element is either:
#'   \itemize{
#'     \item an R list (e.g. \code{list(cs1=..., cs2=..., …)}) of integer vectors,
#'     \item \code{NULL} or an empty list, meaning no confidence sets returned.
#'   }
#'   Each inner integer vector holds the indices in a single confidence set.
#' @param true_active Integer vector of the true active variable indices.
#' @param decimal Nonnegative integer. 
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{set}{Character: each unique simulation pattern rendered as
#'     comma‐separated set literals, e.g.\  "\{1,4\}, \{2,3\}".  An empty
#'     pattern is shown as "\{\}".}
#'   \item{count}{Integer: number of simulations that produced exactly that
#'     pattern of sets.}
#'   \item{percent}{Numeric: \code{count} divided by \code{n_sims}.}
#'   \item{cover}{Logical: \code{TRUE} if the union of all sets in the pattern
#'     contains *every* element of \code{true_active}.}
#' }
#' Rows are sorted in descending order of \code{count}.
#'
#' @export
summarize_cs <- function(cs_list, true_active, decimal=2) {
  n_sims <- length(cs_list)

  # normalize a single simulation's list-of-sets into a canonical, sorted string
  keyify <- function(L) {
    if (is.null(L) || length(L) == 0L) {
      return("{}")
    }
    # extract and normalize each set
    sets <- lapply(L, function(s) {
      s <- unique(as.integer(s))
      if (length(s) == 0L) return(integer(0))
      sort(s)
    })
    # order by size, then lexicographically
    lengths <- vapply(sets, length, integer(1))
    keys    <- vapply(sets, function(v) paste(v, collapse=","), character(1))
    ord     <- order(lengths, keys)
    sets    <- sets[ord]
    # render back to "{...}" and join
    rendered <- vapply(sets, function(v) {
      if (length(v) == 0L) return("{}")
      paste0("{", paste(v, collapse=","), "}")
    }, character(1))
    paste(rendered, collapse = ", ")
  }

  # build keys for all sims
  keys <- vapply(cs_list, keyify, character(1))

  # tabulate unique patterns
  uniq_keys <- unique(keys)
  counts    <- as.integer(table(factor(keys, levels = uniq_keys)))
  percents  <- counts / n_sims

  # compute coverage: union of all sets in a pattern covers true_active?
  is_covered <- function(L) {
    if (is.null(L) || length(L) == 0L) return(FALSE)
    u <- unique(unlist(L, recursive = FALSE, use.names = FALSE))
    all(true_active %in% u)
  }
  # find a representative for each unique key
  first_occ <- match(uniq_keys, keys)
  covers    <- vapply(first_occ, function(i) is_covered(cs_list[[i]]), logical(1))

  # assemble data.frame
  df <- data.frame(
    set     = uniq_keys,
    count   = counts,
    percent = round(percents, decimal),
    cover   = covers,
    stringsAsFactors = FALSE
  )
  df[order(-df$count), , drop = FALSE]
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
#' @param decimal Nonnegative integer. 
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
summarize_coef <- function(sims_coef, true_theta, decimal=2) {
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
    mean = round(mean_est,decimal),
    ssd  = round(ssd_est,decimal),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#’ Summarize Confidence Sets from a glmsusie Fit
#’
#’ @param object A fitted model object of class “glmsusie”
#’ @return A named list with:
#’   * \code{n_sets}: total number of confidence sets  
#’   * \code{avg_size}: average size of those sets  
#’   * \code{pct_singleton}: percent of sets of size 1  
#’   * \code{max_size}: size of the largest set  
#’ @export
cs_statistics <- function(object) {
  # pull out the list of sets
  sets <- object$cs$sets
  
  # if no sets, return all zeros
  if (is.null(sets) || length(sets) == 0L) {
    return(list(
      n_sets        = 0L,
      avg_size      = 0,
      pct_singleton = 0,
      max_size      = 0L
    ))
  }
  
  # compute sizes
  sizes <- vapply(sets, length, integer(1))
  n_sets        <- length(sizes)
  avg_size      <- mean(sizes)
  pct_singleton <- mean(sizes == 1) * 100
  max_size      <- max(sizes)
  
  list(
    n_sets        = n_sets,
    avg_size      = avg_size,
    pct_singleton = pct_singleton,
    max_size      = max_size
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