#' Combine p-values using Fisher's Method
#'
#' @description
#' Aggregates a vector of p-values using Fisher's method, which tests the union null
#' hypothesis that all component tests are null. This method is powerful when multiple
#' weak signals are present across independent or weakly correlated tests.
#'
#' @param pvals A numeric vector of p-values (typically for the same predictor across models or effects).
#'
#' @return A scalar: the Fisher-combined p-value.
#'
#' @details
#' The test statistic is \eqn{T = -2 \sum \log(p_i)}, which follows a
#' chi-squared distribution with \eqn{2k} degrees of freedom under the null,
#' where \eqn{k} is the number of p-values.
#'
#' @importFrom stats pchisq
#' 
#' @examples
#' \dontrun{
#' p <- c(0.01, 0.2, 0.04)
#' combine_fisher(p)
#' }
#' @name combine_fisher
#' @export
combine_fisher <- function(pvals) {
  stat <- -2 * sum(log(pmax(pvals, 1e-300)))
  stats::pchisq(stat, df = 2 * length(pvals), lower.tail = FALSE)
}

#' Combine p-values using Simes' Method
#'
#' @description
#' Aggregates a vector of p-values using Simes' method, which controls the family-wise error rate
#' under the union null hypothesis. This is particularly useful when testing whether at least one
#' of several related tests shows a significant effect.
#'
#' @param pvals A numeric vector of p-values (e.g., from L single-effect models for one predictor).
#'
#' @return A scalar: the Simes-combined p-value.
#'
#' @details
#' The Simes p-value is computed as:
#' \deqn{\min_{i=1}^L \frac{L \cdot p_{(i)}}{i}},
#' where \eqn{p_{(i)}} are the ordered p-values.
#'
#' @examples
#' \dontrun{
#' p <- c(0.01, 0.2, 0.04)
#' combine_simes(p)
#' }
#' @name combine_simes
#' @export
combine_simes <- function(pvals) {
  L <- length(pvals)
  pvals_sorted <- sort(pvals)
  min(L * pvals_sorted / seq_len(L))
}

#' Select Columns with Significant and Non-Diffuse Credible Sets
#'
#' @description
#' This function identifies columns in a probability matrix that are:
#' (1) statistically significant based on a column-wise p-value aggregation method (e.g., Simes),
#' and (2) concentrated, in the sense that the credible set formed from their PIPs does not
#' include all variables. If no column is significant, a fallback selects the column whose
#' p-values are most non-uniform.
#'
#' @details
#' Columns are first assessed using \code{combine_simes} to combine p-values across rows.
#' Among the significant columns (p-value < \code{alpha}), each is further evaluated by forming
#' a credible set from its posterior inclusion probabilities (PIPs). The credible set consists
#' of variables ranked by PIP and included until the cumulative sum reaches at least \code{1 - alpha}.
#'
#' If the credible set includes **all variables**, the column is excluded as too diffuse.
#' This prevents keeping effects that do not concentrate signal on a subset of variables.
#'
#' If no columns remain after filtering, the function selects the column with the most
#' non-uniform p-value distribution (based on variance from mean).
#'
#' @param pmp A numeric matrix of posterior inclusion probabilities (PIPs), dimension \code{p × L}.
#' @param pval A numeric matrix of p-values, dimension \code{p × L}, where each column corresponds to a single effect.
#' @param alpha Significance threshold for both p-value aggregation and credible set construction (default = 0.05).
#'
#' @return A logical vector of length \code{L} indicating selected columns:
#' \itemize{
#'   \item \code{TRUE} if a column is retained (significant and non-diffuse),
#'   \item \code{FALSE} otherwise.
#' }
#' Column names are preserved if present. At least one column is always selected.
#'
#' @note
#' Assumes each column of \code{pval} and \code{pmp} corresponds to the same single-effect model.
#'
#' @name iskept
#' @export
iskept <- function(pmp, pval, alpha = 0.05) {
  # Step 1: Aggregate p-values column-wise using Simes' method
  pvals <- apply(pval, 2, combine_simes)
  keep <- pvals < alpha

  # # Step 2: For kept columns, form credible sets and remove overly diffuse ones
  # if (any(keep)) {
  #   p <- nrow(pmp)
  #   for (j in which(keep)) {
  #     pip_j <- pmp[, j]
  #     ranked <- order(pip_j, decreasing = TRUE)
  #     cum_pip <- cumsum(pip_j[ranked])
  #     threshold_idx <- which(cum_pip >= (1 - alpha))[1]
  #     selected_idx <- ranked[seq_len(threshold_idx)]

  #     # If credible set contains all variables, mark column as not kept
  #     if (length(selected_idx) == p) {
  #       keep[j] <- FALSE
  #     }
  #   }
  # }

  # Step 3: Fallback – keep the most non-uniform column if none are kept
  if (!any(keep)) {
    dists <- apply(pval, 2, function(p) sum((p - mean(p))^2, na.rm = TRUE))
    keep[which.max(dists)] <- TRUE
  }

  # Preserve column names
  if (!is.null(colnames(pval))) {
    names(keep) <- colnames(pval)
  }

  return(keep)
}



#' Extract Model Coefficients from a glmsusie Object
#'
#' @description
#' Extracts coefficients from a glmsusie model fit by summing the coefficients
#' across all single effects for each predictor.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param intercept Logical; if TRUE, include the intercept term as the first element
#'        of the returned vector. Default is FALSE.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of coefficients, with names corresponding to the column
#'         names of the design matrix. If \code{intercept = TRUE}, the first element
#'         will be named "(Intercept)".
#'
#' @export
coef.glmsusie <- function(object, intercept = FALSE, ...) {
  # Sum coefficients across all kept effects
  if (is.null(object$kept) || !any(object$kept)) {
    # If no effects kept, return all zeros
    theta <- rep(0, nrow(object$theta))
    names(theta) <- rownames(object$theta)
  } else {
    theta <- rowSums(object$theta)
  }
  
  # Include intercept if requested and available
  if (intercept && !is.null(object$intercept)) {
    result <- c(object$intercept, theta)
    names(result) <- c("(Intercept)", names(theta))
    return(result)
  } else {
    return(theta)
  }
}

#' Model Predictions from a glmsusie Object
#'
#' @description
#' Produces predictions from a glmsusie model for new data or the original training data.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param newdata An optional data frame or matrix in which to look for variables with
#'        which to predict. If omitted, the fitted values are returned.
#' @param type Character string specifying the type of prediction:
#'        \itemize{
#'          \item "link": predictions on the scale of the linear predictors
#'          \item "response": predictions on the scale of the response
#'          \item "terms": the contributions of individual terms
#'        }
#' @param ... Additional arguments passed to the family's link function.
#'
#' @return Vector or matrix of predictions, depending on the type.
#'
#' @export
predict.glmsusie <- function(object, newdata = NULL, type = c("link", "response", "terms"), ...) {
  type <- match.arg(type)
  
  # Get coefficients
  theta <- coef(object)
  
  # If no newdata, use original data
  if (is.null(newdata)) {
    X <- object$X
  } else {
    # Otherwise, check that newdata has the right variables
    if (is.data.frame(newdata)) {
      X <- as.matrix(newdata[, rownames(object$theta), drop = FALSE])
    } else if (is.matrix(newdata)) {
      X <- newdata
    } else {
      stop("newdata must be a data frame or matrix")
    }
  }
  
  # Calculate linear predictor
  lp <- X %*% theta
  if (!is.null(object$intercept)) {
    lp <- lp + object$intercept
  }
  
  # Return predictions based on type
  if (type == "link") {
    return(drop(lp))
  } else if (type == "response") {
    # For Cox regression, return exp(lp) (relative hazard)
    if (is.list(object$family) && object$family$family == "cox") {
      return(drop(exp(lp)))
    } else {
      # For GLMs, use the inverse link function
      return(drop(object$family$linkinv(lp)))
    }
  } else if (type == "terms") {
    # Return contribution of each term
    terms <- X * rep(theta, each = nrow(X))
    colnames(terms) <- names(theta)
    if (!is.null(object$intercept)) {
      terms <- cbind("(Intercept)" = rep(object$intercept, nrow(X)), terms)
    }
    return(terms)
  }
}


#' Extract Residuals from a glmsusie Object
#'
#' @description
#' Extracts residuals from a glmsusie model fit.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param type Character string indicating the type of residuals to be returned:
#'        \itemize{
#'          \item "deviance": deviance residuals
#'          \item "pearson": Pearson residuals
#'          \item "response": response residuals (observed - fitted)
#'          \item "working": working residuals
#'        }
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of residuals.
#'
#' @export
residuals.glmsusie <- function(object, type = c("deviance", "pearson", "response", "working"), ...) {
  type <- match.arg(type)
  
  # Get fitted values
  fitted_values <- predict(object, type = "response")
  
  # Extract response
  y <- object$y
  
  # Handle different model families
  if (is.list(object$family) && object$family$family == "cox") {
    stop("Residuals are not yet implemented for Cox models")
  } else {
    # For GLMs
    family <- object$family
    
    if (type == "response") {
      # Response residuals (observed - fitted)
      return(y - fitted_values)
    } else if (type == "pearson") {
      # Pearson residuals
      var_func <- family$variance
      variance <- var_func(fitted_values)
      return((y - fitted_values) / sqrt(variance))
    } else if (type == "deviance") {
      # Deviance residuals
      dev_resids <- family$dev.resids(y, fitted_values, rep(1, length(y)))
      return(sign(y - fitted_values) * sqrt(dev_resids))
    } else if (type == "working") {
      # Working residuals
      eta <- predict(object, type = "link")
      mu.eta <- family$mu.eta
      d.eta <- mu.eta(eta)
      return((y - fitted_values) / d.eta)
    }
  }
}

#' Summarize a glmsusie Model Fit
#'
#' @description
#' Produces a summary of a glmsusie model fit, including coefficients,
#' confidence sets, and model fit statistics.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param ... Additional arguments (not used).
#'
#' @return A list with class "summary.glmsusie" containing various summary statistics.
#'
#' @export
summary.glmsusie <- function(object, ...) {
  # Extract coefficients
  theta <- coef(object)
  
  # Calculate inclusion probabilities
  marginal_probs <- object$pip
  
  # Prepare coefficient table
  coef_table <- data.frame(
    Estimate = theta,
    PIP = marginal_probs,
    row.names = names(theta)
  )
  
  # Sort by absolute coefficient value
  coef_table <- coef_table[order(coef_table$PIP, decreasing = TRUE), ]
  
  # Information about convergence
  converged <- object$niter < object$max_iter
  iter_msg <- if(converged) "Blockwise coordinate algorithm converged" else "Blockwise coordinate algorithm  did not converge"
  
  # Create summary object
  result <- list(
    call = object$call,
    family = object$family,
    coefficients = coef_table,
    confidence_sets = object$cs$sets,
    coverage = object$coverage,
    claimed = object$cs$claimed,
    effects_kept = object$kept,
    dispersion = object$dispersion,
    iterations = object$niter,
    converged = converged,
    elapsed = object$elapsed
  )
  
  class(result) <- "summary.glmsusie"
  return(result)
}

#' Print a glmsusie Model Summary
#'
#' @param x A summary.glmsusie object created by \code{summary.glmsusie()}.
#' @param digits Number of significant digits to use for printing.
#' @param max_sets Maximum number of confidence sets to display (default: 5).
#' @param ... Additional arguments (not used).
#'
#' @return The summary object invisibly.
#'
#' @export
print.summary.glmsusie <- function(x, digits = max(3, getOption("digits") - 3), 
                               max_sets = 5, ...) {
  cat("\nCall:\n")
  print(x$call)
  
  cat("\nFamily:", if(is.list(x$family)) x$family$family else x$family$family(), "\n")
  
  cat("\nCoefficients: (sorted by PIP)\n")
  coef_table <- x$coefficients
  coef_table$Estimate <- coef_table$Estimate
  coef_table$PIP <- round(coef_table$PIP, digits)
  print(head(coef_table, 10))
  if (nrow(coef_table) > 10) {
    cat("... (", nrow(coef_table) - 10, " more coefficients not shown)\n", sep = "")
  }
  
  if (!is.null(x$confidence_sets) && length(x$confidence_sets) > 0) {
    cat(paste0("\n", round(x$coverage * 100), "% Confidence Sets:\n"))
    n_sets <- length(x$confidence_sets)
    sets_to_show <- min(max_sets, n_sets)
    
    # Create a data frame for confidence sets
    cs_table <- data.frame(
      Set = sapply(x$confidence_sets, function(cs) paste0("{", paste(cs, collapse = ", "), "}")),
      Coverage = round(x$claimed, digits)
    )
    
    # Display the first max_sets rows
    print(head(cs_table, sets_to_show))
    
    if (n_sets > max_sets) {
      cat("... (", n_sets - max_sets, " more confidence sets not shown)\n", sep = "")
    }
  }
  
  if (!is.null(x$dispersion) && x$dispersion != 1) {
    cat("\nDispersion parameter:", round(x$dispersion, digits), "\n")
  }
  
  cat("\nModel converged after", x$iterations, "iterations.\n")
  cat("Computation time:", round(x$elapsed, 2), "seconds.\n")
  
  invisible(x)
}

#' Print a glmsusie Model Object
#'
#' @param x A glmsusie object.
#' @param ... Additional arguments passed to print.
#'
#' @return The object invisibly.
#'
#' @export
print.glmsusie <- function(x, ...) {
  cat("glmsusie (Generalized Linear Model with Confidence Sets)\n")
  
  cat("\nCall:\n")
  print(x$call)
  
  family_name <- if(is.list(x$family) && !is.null(x$family$family)) {
    x$family$family
  } else {
    x$family$family()
  }
  cat("\nFamily:", family_name, "\n")
  
  cat("\nNumber of observations:", nrow(x$X), "\n")
  cat("Number of predictors:", ncol(x$X), "\n")
  cat("Number of effects retained:", sum(x$kept), "out of", length(x$kept), "\n")
  
  non_zero <- sum(abs(coef(x)) > 0)
  cat("Non-zero coefficients:", non_zero, "\n")
  
  if (!is.null(x$cs$sets)) {
    cat("Number of confidence sets:", length(x$cs$sets), "\n")
  }
  
  cat("\nUse summary() for more detailed information.\n")
  
  invisible(x)
}

#' Plot inclusion probabilities for a set of "credible sets"
#'
#' @param prob_mat Numeric matrix (r × p) of probabilities (each row sums to 1).
#' @param row_labels Optional character vector of length r. Default `paste0("CS", 1:r)`.
#' @param var_labels Optional character or numeric vector of length p. Default 1:p.
#' @param alpha Coverage level for credible sets. Default 0.05 (95% coverage).
#' @return A ggplot2 object showing credible sets.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @export
plot_cs_matrix <- function(prob_mat,
                           row_labels = paste0("CS", seq_len(nrow(prob_mat))),
                           var_labels = seq_len(ncol(prob_mat)),
                           alpha = 0.05) {
  stopifnot(is.matrix(prob_mat))
  r <- nrow(prob_mat); p <- ncol(prob_mat)
  if (length(row_labels) != r) stop("row_labels must match rows")
  if (length(var_labels) != p) stop("var_labels must match cols")
  
  # Find credible sets for each row
  cs_data <- data.frame()
  for (i in 1:r) {
    # Sort probabilities in descending order and get indices
    sorted_indices <- order(prob_mat[i, ], decreasing = TRUE)
    sorted_probs <- prob_mat[i, sorted_indices]
    
    # Find cumulative sum until > 1-alpha
    cumsum_probs <- cumsum(sorted_probs)
    threshold <- 1 - alpha
    
    # Find the smallest set with cumulative probability > threshold
    n_vars <- which(cumsum_probs > threshold)[1]
    if (is.na(n_vars)) n_vars <- p  # if never exceeds threshold, take all
    
    cs_vars <- sorted_indices[1:n_vars]
    
    # Add to data frame
    if (length(cs_vars) > 0) {
      cs_data <- rbind(cs_data, data.frame(
        Row = factor(rep(i, length(cs_vars)), levels = 1:r, labels = row_labels),
        Var = cs_vars
      ))
    }
  }
  
  # Create plot
  if (nrow(cs_data) == 0) {
    stop("No variables in any credible sets")
  }
  
  x_breaks <- sort(unique(cs_data$Var))
  
  p <- ggplot2::ggplot(cs_data, ggplot2::aes(x = .data$Var, y = .data$Row)) +
    ggplot2::geom_point(size = 4, shape = 7) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::labs(
      title = paste0((1-alpha)*100, "% Credible sets"),
      x = "Variable",
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
    return(p)
}

#' Plot a glmsusie Model Object
#'
#' @description
#' Produces various plots for visualizing glmsusie model results.
#'
#' @param x A glmsusie object.
#' @param which Integer or character string specifying which plot to create:
#'        \itemize{
#'          \item 1 or "coefficients": Plot of coefficient estimates
#'          \item 2 or "probabilities": Plot of inclusion probabilities
#'          \item 3 or "sets": Visualization of confidence sets
#'        }
#' @param n_top Integer specifying how many top coefficients to display.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot2 object.
#' @importFrom rlang .data
#' @export
plot.glmsusie <- function(x, which = c("coefficients", "probabilities", "sets"), 
                       n_top = 20, ...) {
  # Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for plotting")
  }
  
  # Match plot type
  if (is.numeric(which)) {
    which <- c("coefficients", "probabilities", "sets")[which]
  } else {
    which <- match.arg(which)
  }
  
  # Get coefficients and probabilities
  theta <- coef(x)
  probs <- x$pip
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Variable = names(theta),
    Coefficient = theta,
    Probability = probs,
    AbsCoef = abs(theta)
  )
  
  # Sort by absolute coefficient value
  plot_data <- plot_data[order(plot_data$AbsCoef, decreasing = TRUE), ]
  
  # Limit to top n
  if (nrow(plot_data) > n_top) {
    plot_data <- plot_data[1:n_top, ]
  }
  
  # Reverse order for plotting (so largest is at the top)
  plot_data$Variable <- factor(plot_data$Variable, 
                              levels = rev(plot_data$Variable))
  
  # Create the requested plot
  if (which == "coefficients") {
    raw_vars <- levels(plot_data$Variable)
    idx      <- as.integer(gsub("^X", "", raw_vars))
    newlv    <- paste0("X[", idx, "]")
    levels(plot_data$Variable) <- newlv
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Coefficient, y = .data$Variable)) +
      ggplot2::geom_point(size = 2.5, shape = 7) +  # Changed to cross shape (4)
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$Coefficient, 
                                        yend = .data$Variable), color = "gray90", linewidth = 0.2) +
      ggplot2::labs(title = "glmsusie Coefficient Estimates",
                   x = "Coefficient Value", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        panel.grid.major.x = ggplot2::element_line(color = "gray90"),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
      scale_y_discrete(labels=parse(text=levels(plot_data$Variable)))
    
  } else if (which == "probabilities") {
    raw_vars <- levels(plot_data$Variable)
    idx      <- as.integer(gsub("^X", "", raw_vars))
    newlv <- paste0("X[", idx, "]")
    levels(plot_data$Variable) <- newlv
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Probability, y = .data$Variable)) +
      ggplot2::geom_point(size = 2.5, shape = 7) +  # Changed to cross shape (4)
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$Probability, 
                                        yend = .data$Variable), color = "gray90", linewidth = 0.2) +
      ggplot2::labs(title = "glmsusie Inclusion Probabilities",
                   x = "Probability", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        panel.grid.major.x = ggplot2::element_line(color = "gray90"),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
      ) + 
      scale_y_discrete(labels=parse(text=levels(plot_data$Variable)))
    
  } else if (which == "sets") {
    # Check if confidence sets exist
    if (is.null(x$cs$sets)) {
      stop("No confidence sets available in the model object")
    }
    
    # Extract posterior model probabilities for kept variables
    if (is.null(x$pmp) || is.null(x$kept)) {
      stop("Posterior model probabilities or kept variables not available")
    }
    
    # Get posterior model probabilities matrix
    pmp_matrix <- t(x$pmp[, x$kept, drop = FALSE])
    p <- plot_cs_matrix(pmp_matrix)
  }
  
  return(p)
}

#' Extract Fitted Values from a glmsusie Object
#'
#' @description
#' Returns fitted values from a glmsusie model object.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param ... Additional arguments (not used).
#'
#' @return A vector of fitted values on the scale of the response variable.
#'
#' @export
fitted.glmsusie <- function(object, ...) {
  # Use predict with type "response" to get fitted values
  predict(object, type = "response")
}

#' Extract Log-Likelihood from a glmsusie Object
#'
#' @description
#' Returns the log-likelihood of a glmsusie model object.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param ... Additional arguments (not used).
#'
#' @return Log-likelihood value.
#'
#' @export
logLik.glmsusie <- function(object, ...) {
  # pull out basics
  fam  <- object$family
  X    <- object$X
  y    <- object$y
  n    <- NROW(X)
  # number of nonzero single effects = df for β
  df_beta <- sum(object$kept)
  # include intercept?
  has_int <- !is.null(object$intercept)
  df <- df_beta + as.integer(has_int)
  
  # compute linear predictor
  eta <- predict(object, type = "link")
  
  if ((is.list(fam) && fam$family == "cox") ||
      (is.character(fam) && fam == "cox")) {
    # --- partial log-likelihood for Cox (Breslow) ---
    times  <- y[,1]
    status <- y[,2]
    exp_lp <- exp(eta)
    # sum over events: lp_i - log(sum_{j: t_j >= t_i} exp(lp_j))
    event_ids <- which(status == 1)
    ll <- sum(vapply(event_ids, function(i) {
      lp_i   <- eta[i]
      denom  <- sum(exp_lp[times >= times[i]])
      lp_i - log(denom)
    }, numeric(1)))
  } else {
    mu <- fam$linkinv(eta)
    phi <- object$dispersion
    ll <- switch(fam$family,
      gaussian = {
        # base R Gaussian log-density
        sd <- sqrt(phi)
        sum(stats::dnorm(object$y, mean = mu, sd = sd, log = TRUE))
      },
      binomial = {
        # assume y is 0/1
        sum(stats::dbinom(object$y, size = 1, prob = mu, log = TRUE))
      },
      poisson = {
        sum(stats::dpois(object$y, lambda = mu, log = TRUE))
      },
      Gamma = {
        # R’s Gamma(shape, scale): shape=1/phi, scale=mu*phi
        sh    <- 1/phi
        sc    <- mu * phi
        sum(stats::dgamma(object$y, shape = sh, scale = sc, log = TRUE))
      },
      {
        # fallback: up to constant, −½ sum deviance residuals
        dev <- fam$dev.resids(object$y, mu, rep(1, n))
        -0.5 * sum(dev)
      }
    )
  }
  
  # attach attributes
  attr(ll, "df")   <- df
  attr(ll, "nobs") <- n
  class(ll)        <- "logLik"
  ll
}

#' Extract AIC from a glmsusie Object
#'
#' @description
#' Returns the Akaike Information Criterion (AIC) for a glmsusie model.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param ... Additional arguments (not used).
#' @param k Numeric, the penalty per parameter to be used; default is 2.
#'
#' @return AIC value.
#'
#' @export
AIC.glmsusie <- function(object, ..., k = 2) {
  # Get log-likelihood
  ll <- logLik(object)
  
  # Calculate AIC
  df <- attr(ll, "df")
  aic <- -2 * as.numeric(ll) + k * df
  
  return(aic)
}


#' Extract BIC from a glmsusie Object
#'
#' @description
#' Returns the Bayesian Information Criterion (BIC) for a glmsusie model.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param ... Additional arguments (not used).
#'
#' @return BIC value.
#'
#' @export
BIC.glmsusie <- function(object, ...) {
  # Get log-likelihood
  ll <- logLik(object)
  
  # Calculate BIC
  df <- attr(ll, "df")
  n <- attr(ll, "nobs")
  bic <- -2 * as.numeric(ll) + log(n) * df
  
  return(bic)
}

#' Calculate Deviance for a glmsusie Object
#'
#' @description
#' Returns the deviance of a glmsusie model.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param ... Additional arguments (not used).
#'
#' @return Deviance value.
#'
#' @export
deviance.glmsusie <- function(object, ...) {
  # Extract components needed for deviance calculation
  y <- object$y
  fitted_values <- fitted(object)
  weights <- if (!is.null(object$weights)) object$weights else rep(1, length(y))
  family <- object$family
  
  # Calculate deviance
  dev <- sum(family$dev.resids(y, fitted_values, weights))
  
  return(dev)
}

#' Extract Dispersion Parameter from a glmsusie Object
#'
#' @description
#' Returns the estimated dispersion parameter from a glmsusie model.
#'
#' @param object A fitted model object of class "glmsusie".
#' @param ... Additional arguments (not used).
#'
#' @return Dispersion parameter value.
#'
#' @export
dispersion.glmsusie <- function(object, ...) {
  # Return dispersion if available, otherwise calculate
  if (!is.null(object$dispersion)) {
    return(object$dispersion)
  } else {
    # For GLMs with fixed dispersion (e.g., binomial, poisson)
    if (object$family$family %in% c("binomial", "poisson")) {
      return(1)
    }
    
    # For other GLMs, calculate dispersion from deviance
    dev <- deviance(object)
    df <- df.residual(object)
    disp <- dev / df
    return(disp)
  }
}