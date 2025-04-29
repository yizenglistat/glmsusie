#include "core.h"
#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

/**
 * @title Fit Single-Effect Regression (SER) Across Predictors
 * @description
 * Performs single-effect regression for each predictor variable against the response,
 * supporting both Generalized Linear Models (GLMs) and Cox Proportional Hazards models.
 *
 * @details
 * For each predictor, this function:
 * 1. Fits a model with just that predictor (and optional offset)
 * 2. Computes coefficient estimate, standard error, and p-value
 * 3. Calculates model fit statistics (log-likelihood, BIC)
 * 4. Derives Bayes factors and posterior probabilities
 *
 * All models are compared against a null model with only the offset.
 *
 * @param X NumericMatrix (n × p). Predictor matrix with:
 *        n = number of observations
 *        p = number of predictors
 * @param y Response variable:
 *        - For GLMs: numeric vector of length n
 *        - For Cox: matrix with 2 columns (time, status)
 * @param family Model family specification:
 *        - String: "cox", "gaussian", "binomial", "poisson", etc.
 *        - List with element "family" or "name" containing family string
 * @param offset NumericVector. Known offset term for the linear predictor:
 *        - Length 1: Same offset applied to all observations
 *        - Length n: Observation-specific offsets
 * @param standardize Logical. Whether to standardize predictors before fitting (default: TRUE).
 * @param ties String. Method for handling ties in Cox regression:
 *        - "efron": Efron approximation (default, more accurate)
 *        - "breslow": Breslow approximation (faster)
 *
 * @return DataFrame with columns:
 *   - theta: Estimated coefficient for each predictor
 *   - se: Standard error of the coefficient estimate
 *   - p_value: Two-sided p-value for coefficient significance
 *   - logLik: Log-likelihood of the model
 *   - BIC: Bayesian Information Criterion
 *   - delta_BIC: Difference in BIC compared to minimum BIC (stabilized)
 *   - bf: Bayes factor derived from delta_BIC
 *   - posterior: Posterior probability for variable importance
 *
 * @examples
 * \dontrun{
 * # Example for Gaussian regression
 * X <- matrix(rnorm(100 * 10), 100, 10)
 * y <- rnorm(100)
 * results <- get_ser_fit(X, y, "gaussian", rep(0, 100))
 * 
 * # Example for Cox regression
 * time <- rexp(100)
 * status <- rbinom(100, 1, 0.5)
 * y_cox <- cbind(time, status)
 * results_cox <- get_ser_fit(X, y_cox, "cox", rep(0, 100), ties = "efron")
 * }
 */
// [[Rcpp::export]]
DataFrame get_ser_fit(NumericMatrix X,
                      SEXP y,
                      SEXP family,
                      NumericVector offset,
                      bool standardize = true,
                      std::string ties = "efron") {

  // --- Input dimensions ---
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  const int n = Xmat.n_rows;
  const int p = Xmat.n_cols;

  // --- Convert offset ---
  arma::vec offset_vec(offset.begin(), offset.size());
  if (offset_vec.n_elem == 1) {
    offset_vec = arma::vec(n, arma::fill::value(offset_vec[0]));
  } else if ((int)offset_vec.n_elem != n) {
    stop("`offset` must be length 1 or length n");
  }

  // --- Extract family name ---
  std::string fam_name;
  if (TYPEOF(family) == STRSXP) {
    fam_name = as<std::string>(family);
  } else if (TYPEOF(family) == VECSXP) {
    List fam_list(family);
    if (fam_list.containsElementNamed("family")) {
      fam_name = as<std::string>(fam_list["family"]);
    } else if (fam_list.containsElementNamed("name")) {
      fam_name = as<std::string>(fam_list["name"]);
    } else {
      stop("Family object must contain 'family' or 'name'.");
    }
  } else {
    stop("family must be a List or a string");
  }

  const bool is_cox = (fam_name == "cox");

  // --- Fit null model (intercept + offset only) ---
  double BIC_null;
  if (is_cox) {
    arma::mat y_mat = as<arma::mat>(y);
    List null_fit = null_cox_fit(y_mat, offset_vec, ties);
    BIC_null = as<double>(null_fit["BIC"]);
  } else {
    arma::vec y_vec = as<arma::vec>(y);
    List null_fit = null_glm_fit(y_vec, family, offset_vec);
    BIC_null = as<double>(null_fit["BIC"]);
  }

  // --- Prepare result containers ---
  arma::vec theta(p, arma::fill::zeros);
  arma::vec se(p, arma::fill::zeros);
  arma::vec p_value(p, arma::fill::ones);
  arma::vec logLik(p, arma::fill::zeros);
  arma::vec BIC(p, arma::fill::zeros);

  // --- Loop through each predictor ---
  #pragma omp parallel for if(p > 100)
  for (int j = 0; j < p; ++j) {
    arma::vec xj = Xmat.col(j);

    List res;
    if (is_cox) {
      arma::mat y_mat = as<arma::mat>(y);
      res = get_cox_fit(xj, y_mat, offset_vec, standardize, ties);
    } else {
      arma::vec y_vec = as<arma::vec>(y);
      res = get_glm_fit(xj, y_vec, family, offset_vec, standardize);
    }

    theta[j]   = as<double>(res["theta"]);
    se[j]      = as<double>(res["se"]);
    p_value[j] = as<double>(res["p_value"]);
    logLik[j]  = as<double>(res["logLik"]);
    BIC[j]     = as<double>(res["BIC"]);
  }

  // --- Compute model comparison metrics ---
  arma::vec delta_BIC = BIC - BIC_null;
  delta_BIC -= delta_BIC.min();  // stabilization for numerical precision

  arma::vec bf = arma::exp(-0.5 * delta_BIC);
  arma::vec posterior = bf / arma::accu(bf);

  // --- Build output DataFrame ---
  return DataFrame::create(
    _["theta"]     = theta,
    _["se"]        = se,
    _["p_value"]   = p_value,
    _["logLik"]    = logLik,
    _["BIC"]       = BIC,
    _["delta_BIC"] = delta_BIC,
    _["bf"]        = bf,
    _["posterior"] = posterior
  );
}

/**
 * @title Fit LASER Model (L-Effect Sparse Regression) with Early Stopping
 * 
 * @description
 * Fits the LASER model, which represents a response as a sum of L sparse effects,
 * using block coordinate ascent. The model supports both GLM families and Cox 
 * regression and includes early stopping for efficient computation.
 *
 * @details
 * The LASER model represents the linear predictor as:
 * 
 *    η = intercept + X * Σ(θᵢ)
 * 
 * where θᵢ are L sparse vectors of coefficients. The model is fitted using a block 
 * coordinate ascent algorithm that iteratively updates each effect while holding 
 * others constant. Three update strategies are available:
 * 
 * 1. "cyclic": Update each effect in order (default)
 * 2. "shuffle": Update effects in random order each iteration
 * 3. "greedy": Select the single effect update giving maximum improvement
 * 
 * For each update, the function uses single-effect regression (SER) to compute posterior
 * probabilities and coefficient estimates. The dispersion parameter and intercept are
 * updated after each complete iteration.
 *
 * @param X NumericMatrix (n × p) predictor matrix, where:
 *        n = number of observations
 *        p = number of predictors
 * @param y Response variable:
 *        - For GLMs: numeric vector of length n
 *        - For Cox: matrix with 2 columns (time, status)
 * @param family Model family specification:
 *        - String: "cox", "gaussian", "binomial", "poisson", etc.
 *        - List with element "family" or "name" containing family string
 * @param L Integer. Number of latent effects to include (default = 10).
 *                  If ≤ 0 or > p, automatically set to min(10, p).
 * @param standardize Logical. Whether to standardize predictors (default: TRUE).
 * @param ties String. Method for handling ties in Cox regression (default: "efron").
 * @param max_iter Integer. Maximum number of coordinate ascent iterations (default: 100).
 * @param seed Nullable<int>. Random seed for reproducibility (default: NULL).
 *             If NULL or NA, no seed is set.
 * @param step_size Double. Step size multiplier for updates (default: 1.0).
 *                  Values < 1 provide more conservative updates.
 * @param tol Double. Convergence tolerance on log-likelihood change (default: 1e-6).
 * @param method String. Update strategy: "cyclic", "shuffle", or "greedy" (default: "greedy").
 *
 * @return A list containing:
 *   - total_elapsed: Total execution time (seconds)
 *   - elapsed_per_ser: Average time per SER computation
 *   - n_iter: Number of iterations performed
 *   - final_loglike: Final model log-likelihood
 *   - intercept: Estimated intercept term
 *   - dispersion: Estimated dispersion parameter (for applicable GLMs)
 *   - theta: Matrix (p×L) of estimated coefficients for each effect
 *   - posterior: Matrix (p×L) of posterior probabilities
 *   - p_values: Matrix (p×L) of p-values
 *   - theta_hat: Vector (length p) of summed coefficients across effects
 *
 * @examples
 * \dontrun{
 * # Gaussian example with 10 latent effects
 * X <- matrix(rnorm(100 * 20), 100, 20)
 * y <- rnorm(100)
 * laser_fit <- get_laser_fit(X, y, "gaussian", L = 5, method = "shuffle")
 * 
 * # Binomial example with fewer effects
 * y_bin <- rbinom(100, 1, 0.5)
 * laser_bin <- get_laser_fit(X, y_bin, "binomial", L = 3, seed = 123)
 * 
 * # Cox regression example
 * time <- rexp(100)
 * status <- rbinom(100, 1, 0.5)
 * y_cox <- cbind(time, status)
 * laser_cox <- get_laser_fit(X, y_cox, "cox", L = 5, method = "greedy")
 * }
 * 
 * @note
 * The "greedy" method may be slower but often produces higher-quality models,
 * especially for smaller datasets. For large datasets, "cyclic" offers a good
 * balance of speed and quality.
 */
// [[Rcpp::export]]
List get_laser_fit(NumericMatrix   X,
                   SEXP            y,
                   SEXP            family,
                   int             L            = 10,
                   bool            standardize  = true,
                   std::string     ties         = "efron",
                   std::string     method       = "greedy",
                   int             max_iter     = 100,
                   double          step_size    = 1.0,
                   double          tol          = 1e-6,
                   Nullable<int>   seed         = R_NilValue)
{
  // Set random seed if provided
  if (seed.isNotNull()) {
    Function set_seed("set.seed");
    set_seed(as<int>(seed));
  }

  // Convert R matrix to Armadillo for efficient computation
  arma::mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
  [[maybe_unused]] const int n = Xmat.n_rows;
  const int p = Xmat.n_cols;

  // Validate and adjust number of effects
  if (L <= 0 || L > p) {
    L = std::min(10, p);
  }

  // Extract family name from input
  std::string fam_name;
  if (TYPEOF(family) == STRSXP) {
    fam_name = as<std::string>(family);
  } else if (TYPEOF(family) == VECSXP) {
    List fam_list(family);
    if (fam_list.containsElementNamed("family")) {
      fam_name = as<std::string>(fam_list["family"]);
    } else if (fam_list.containsElementNamed("name")) {
      fam_name = as<std::string>(fam_list["name"]);
    } else {
      stop("Family object must contain 'family' or 'name' element.");
    }
  } else {
    stop("family must be a List or a string.");
  }
  const bool is_cox = (fam_name == "cox");

  // Initialize parameter matrices
  arma::mat theta     = arma::zeros<arma::mat>(p, L);
  arma::mat posterior = arma::zeros<arma::mat>(p, L);
  arma::mat p_values  = arma::zeros<arma::mat>(p, L);

  // Initialize model parameters
  double intercept   = 0.0;
  double dispersion  = 1.0;
  double prev_ll     = R_NegInf;
  double curr_ll     = 0.0;
  int    iter        = 0;

  // Initialize timing variables
  NumericVector t0 = Function("proc.time")();
  const double t_start = t0[2];
  double it_start = 0.0, it_end = 0.0;

  // Prepare indices for effect updates
  std::vector<int> idx(L);
  std::iota(idx.begin(), idx.end(), 0);  // Fill with 0, 1, ..., L-1

  // Main coordinate ascent loop
  for (; iter < max_iter; ++iter) {
    
    // Refresh GLM family if needed (for updated dispersion)
    SEXP reg_fam = family;
    if (!is_cox) {
      Function get_family("get_family");
      reg_fam = get_family(family, Named("dispersion") = dispersion);
    }
    
    // Compute current log-likelihood
    curr_ll = get_loglike(X, y, reg_fam, theta, intercept, dispersion);

    // Start timer for this iteration
    NumericVector t1 = Function("proc.time")();
    it_start = t1[2];

    // Calculate sum of all effects for offset calculation
    arma::vec sum_all = arma::sum(theta, 1);

    // === Greedy update strategy ===
    if (method == "greedy") {
      // Find best single-effect update
      double best_gain = 0.0;
      int    best_l    = -1;
      arma::vec best_col(p);
      NumericVector best_post, best_pv;

      // Try updating each effect and keep the best
      for (int l = 0; l < L; ++l) {
        // Calculate offset without current effect
        arma::vec other = sum_all - theta.col(l);
        arma::vec off   = intercept + Xmat * other;

        // Fit single-effect regression
        List ser = get_ser_fit(X, y, reg_fam, wrap(off),
                               standardize, ties);

        // Extract results
        NumericVector th   = ser["theta"];
        NumericVector post = ser["posterior"];
        NumericVector pv   = ser["p_value"];

        // Apply posterior shrinkage
        arma::vec cand(p);
        for (int j = 0; j < p; ++j) {
          cand[j] = step_size * th[j] * post[j];
        }

        // Evaluate new log-likelihood with this update
        arma::mat theta_tmp = theta;
        theta_tmp.col(l) = cand;
        double ll_new = get_loglike(X, y, family, theta_tmp, intercept, dispersion);

        // Check if this is the best update so far
        double gain = ll_new - curr_ll;
        if (gain > best_gain) {
          best_gain   = gain;
          best_l      = l;
          best_col    = cand;
          best_post   = post;
          best_pv     = pv;
        }
      }

      // Stop if no improvement found
      if (best_gain < tol) {
        break;
      }
            
      // Apply the best update
      theta.col(best_l)     = best_col;
      posterior.col(best_l) = as<arma::vec>(best_post);
      p_values.col(best_l)  = as<arma::vec>(best_pv);
      curr_ll += best_gain;
    }
    // === Cyclic or shuffled update strategy ===
    else {
      // Shuffle effect order if requested
      if (method == "shuffle") {
        std::shuffle(idx.begin(), idx.end(), std::mt19937{std::random_device{}()});
      }
      
      // Update each effect in sequence
      for (int ii = 0; ii < L; ++ii) {
        int l = idx[ii];
        
        // Calculate offset without current effect
        arma::vec other = sum_all - theta.col(l);
        arma::vec off   = intercept + Xmat * other;

        // Fit single-effect regression
        List ser = get_ser_fit(X, y, reg_fam, wrap(off),
                               standardize, ties);

        // Extract and apply results with posterior shrinkage
        NumericVector th   = ser["theta"];
        NumericVector post = ser["posterior"];
        NumericVector pv   = ser["p_value"];

        for (int j = 0; j < p; ++j) {
          theta(j, l)     = step_size * th[j] * post[j];
          posterior(j, l) = post[j];
          p_values(j, l)  = pv[j];
        }
        
        // Update sum for next effect
        sum_all = other + theta.col(l);
      }
      
      // Recompute log-likelihood after all updates
      curr_ll = get_loglike(X, y, family, theta, intercept, dispersion);
    }

    // Stop timing for this iteration
    NumericVector t2 = Function("proc.time")();
    it_end = t2[2];

    // Update intercept & dispersion after all effects
    arma::vec off_no_int = Xmat * arma::sum(theta, 1);
    intercept  = update_intercept(y, family, off_no_int);
    arma::vec full_off  = intercept + off_no_int;
    dispersion = update_dispersion(y, family, full_off);

    // Check convergence
    if (std::abs(curr_ll - prev_ll) < tol) {
      break;
    }
    prev_ll = curr_ll;
  }

  // Calculate total elapsed time
  NumericVector t3 = Function("proc.time")();
  const double t_end = t3[2];

  // Return results
  return List::create(
    _["total_elapsed"]   = t_end - t_start,
    _["elapsed_per_ser"] = (it_end - it_start) / L,
    _["n_iter"]          = iter + 1,
    _["final_loglike"]   = curr_ll,
    _["intercept"]       = intercept,
    _["dispersion"]      = dispersion,
    _["theta"]           = theta,
    _["posterior"]       = posterior,
    _["p_values"]        = p_values,
    _["theta_hat"]       = arma::sum(theta, 1)
  );
}

RCPP_MODULE(model) {
  Rcpp::function("get_laser_fit", &get_laser_fit);
  Rcpp::function("get_ser_fit",   &get_ser_fit);
}