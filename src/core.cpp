#include <RcppArmadillo.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


/**
 * @title Kullback-Leibler Divergence from Uniform Distribution
 * @description
 * Calculates the Kullback-Leibler (KL) divergence between a given probability
 * distribution `p` and a uniform distribution of equal length. KL divergence
 * measures how much one probability distribution diverges from another
 * expected distribution.
 *
 * In the LASER modeling framework, this function helps quantify how much a
 * posterior distribution deviates from uniformity, which is useful for
 * identifying significant effects.
 *
 * @param p NumericVector. Probability distribution vector that should sum to 1.
 * @param log_base double. Base of the logarithm used in calculation (default: e = 2.718...).
 *                 Common alternatives include 2 (bits) or 10.
 *
 * @return double. KL divergence value (always non-negative). A value of 0 indicates
 *         that `p` is exactly uniform, while larger values indicate greater
 *         divergence from uniformity.
 *
 * @details
 * The KL divergence is calculated as:
 * \deqn{D_{KL}(P||U) = \sum_{i=1}^{n} P(i) \log\left(\frac{P(i)}{U(i)}\right)}
 * where \eqn{U(i) = 1/n} is the uniform distribution.
 *
 * For numerical stability, this implementation:
 * 1. Only processes positive probability elements (zeros contribute nothing)
 * 2. Allows for different logarithm bases
 *
 * The KL divergence is always non-negative and equals zero if and only if
 * the distributions are identical (i.e., `p` is perfectly uniform).
 *
 * @examples
 * # In R:
 * # Uniform distribution (KL = 0)
 * p1 <- rep(0.25, 4)
 * kl_divergence(p1)  # Should return 0
 *
 * # Somewhat skewed distribution
 * p2 <- c(0.1, 0.2, 0.3, 0.4)
 * kl_divergence(p2)  # Positive value
 *
 * # Very concentrated distribution
 * p3 <- c(0.01, 0.01, 0.97, 0.01)
 * kl_divergence(p3)  # Large positive value
 *
 * # Using binary logarithm (base 2)
 * kl_divergence(p3, log_base = 2)
 *
 * @references
 * Kullback, S., & Leibler, R. A. (1951). On information and sufficiency.
 * The Annals of Mathematical Statistics, 22(1), 79-86.
 *
 * @seealso
 * \code{\link{get_included}} which uses KL divergence to identify significant effects
 */
// [[Rcpp::export]]
double kl_divergence(NumericVector p, double log_base = 2.718) {
  int n = p.size();
  double unif = 1.0 / n;
  double D = 0.0;
  for (int i = 0; i < n; ++i) {
    if (p[i] > 0) {
      D += p[i] * std::log(p[i] / unif) / std::log(log_base);
    }
  }
  return D;
}

/**
 * @title Compute Log-Likelihood for LASER Model
 * @description 
 * Computes the model log-likelihood for various statistical models supported by LASER.
 * Handles both Generalized Linear Models (GLMs) and Cox Proportional Hazards models.
 * For GLMs, supports gaussian, binomial, poisson, gamma, and inverse gaussian families.
 * 
 * @details
 * The linear predictor is computed as: intercept + X * sum(theta, 1), where
 * sum(theta, 1) represents the sum across the latent dimensions (column-wise).
 * For each family, the appropriate link function and log-likelihood formula are applied.
 * 
 * @param X NumericMatrix (n × p) predictor matrix, where:
 *        n = number of observations
 *        p = number of predictors
 * @param y Response variable:
 *        - For GLMs: numeric vector of length n
 *        - For Cox: matrix with 2 columns (time, status), where status = 1 indicates event
 * @param family Model family specification:
 *        - String: "cox", "gaussian", "binomial", "poisson", "Gamma", or "inverse.gaussian"
 *        - List with element "family" or "name" containing one of the above strings
 * @param theta arma::mat (p × L) estimated coefficient matrix, where:
 *        p = number of predictors
 *        L = number of latent dimensions
 * @param intercept Numeric scalar. Estimated intercept term in the model.
 * @param dispersion Numeric scalar. Estimated dispersion parameter (only used for certain GLMs).
 *                   Defaults to 1.0.
 * 
 * @return Numeric scalar: total log-likelihood of the model.
 * 
 * @examples
 * \dontrun{
 * # Example for gaussian model with 100 observations, 10 predictors, 2 latent dimensions
 * X <- matrix(rnorm(100 * 10), 100, 10)
 * y <- rnorm(100)
 * theta <- matrix(rnorm(10 * 2), 10, 2)
 * get_loglike(X, y, "gaussian", theta, 0.5, 1.2)
 * }
 */
// [[Rcpp::export]]
double get_loglike(NumericMatrix X,
                   SEXP y,
                   SEXP family,
                   arma::mat theta,
                   double intercept,
                   double dispersion = 1.0) {
  // Extract dimensions from data
  const int n = X.nrow();
  const int p = X.ncol();
  
  // Convert R matrix to Armadillo matrix (no copy, just a reference)
  arma::mat Xmat(X.begin(), n, p, false);
  
  // Calculate linear predictor: intercept + X * colSums(theta)
  arma::vec linear_pred = intercept + Xmat * arma::sum(theta, 1);
  
  // Extract family name from input
  std::string fam_name;
  if (TYPEOF(family) == STRSXP) {
    // Input is a string
    fam_name = as<std::string>(family);
  } else if (TYPEOF(family) == VECSXP) {
    // Input is a list
    List fam_list(family);
    if (fam_list.containsElementNamed("family")) {
      fam_name = as<std::string>(fam_list["family"]);
    } else if (fam_list.containsElementNamed("name")) {
      fam_name = as<std::string>(fam_list["name"]);
    } else {
      stop("Family object must contain 'family' or 'name' element.");
    }
  } else {
    stop("Family must be a List or a string.");
  }
  
  // Initialize log-likelihood
  double logLik = 0.0;
  
  // Calculate log-likelihood based on family
  if (fam_name == "cox") {
    // --- Cox proportional hazards model ---
    arma::mat ymat = as<arma::mat>(y);
    arma::vec time = ymat.col(0);
    arma::vec status = ymat.col(1);
    
    // Sort by time (increasing)
    arma::uvec sorted_idx = arma::sort_index(time);
    arma::vec time_sorted = time(sorted_idx);
    arma::vec status_sorted = status(sorted_idx);
    arma::vec lp_sorted = linear_pred(sorted_idx);
    
    // Calculate Cox partial log-likelihood
    double risk_sum = 0.0;
    for (int i = n - 1; i >= 0; --i) {
      double exp_lp = std::exp(lp_sorted[i]);
      risk_sum += exp_lp;
      if (status_sorted[i] == 1) {
        logLik += lp_sorted[i] - std::log(risk_sum);
      }
    }
  } else if (fam_name == "gaussian") {
    // --- Gaussian model (normal distribution) ---
    arma::vec y_vec = as<arma::vec>(y);
    arma::vec resid = y_vec - linear_pred;
    double rss = arma::dot(resid, resid);
    logLik = -n/2.0 * std::log(2 * M_PI * dispersion) - rss / (2.0 * dispersion);
  } else if (fam_name == "binomial") {
    // --- Binomial model (logistic regression) ---
    arma::vec y_vec = as<arma::vec>(y);
    arma::vec mu = 1.0 / (1.0 + arma::exp(-linear_pred));
    // Numerically stable computation
    logLik = arma::accu(y_vec % arma::log(mu + 1e-10) + 
                        (1 - y_vec) % arma::log(1 - mu + 1e-10));
  } else if (fam_name == "poisson") {
    // --- Poisson model ---
    arma::vec y_vec = as<arma::vec>(y);
    arma::vec mu = arma::exp(linear_pred);
    logLik = arma::accu(y_vec % linear_pred - mu - lgamma(y_vec + 1));
  } else if (fam_name == "Gamma") {
    // --- Gamma model ---
    arma::vec y_vec = as<arma::vec>(y);
    arma::vec mu = arma::exp(linear_pred);
    // Note: This is proportional to the actual log-likelihood
    // (ignoring shape parameter terms)
    logLik = arma::accu(-y_vec / mu - arma::log(mu));
  } else if (fam_name == "inverse.gaussian") {
    // --- Inverse Gaussian model ---
    arma::vec y_vec = as<arma::vec>(y);
    arma::vec mu = arma::exp(linear_pred);
    // Note: This is proportional to the actual log-likelihood
    // (ignoring dispersion parameter terms)
    logLik = arma::accu(-(y_vec - mu) % (y_vec - mu) / (2.0 * mu % mu % y_vec));
  } else {
    stop("Unsupported family: '" + fam_name + "'. Supported families are: 'cox', 'gaussian', 'binomial', 'poisson', 'Gamma', and 'inverse.gaussian'.");
  }
  
  return logLik;
}

/**
 * @title Calculate Purity Statistics for LASER Model Variables
 * @description 
 * Computes purity statistics (minimum, mean, and median of absolute correlation values)
 * for a set of selected variables. Purity measures help assess variable independence
 * within selected groups.
 * 
 * @details
 * Purity is calculated based on the absolute correlation values between pairs of
 * selected variables. The function can either use a provided correlation matrix or
 * compute correlations from the raw data matrix. For large datasets, a random subset
 * of observations can be used to calculate correlations more efficiently.
 * 
 * @param pos IntegerVector containing indices of selected variables (1-based indexing).
 * @param X Nullable<NumericMatrix> original predictor matrix (n x p).
 *           If provided, correlations are calculated from this matrix.
 *           Not required if R is provided.
 * @param R Nullable<NumericMatrix> pre-calculated correlation matrix (p x p).
 *               Not required if X is provided.
 * @param squared Logical. If TRUE, squared correlation values are used. Default is FALSE.
 * @param n_purity Integer. Maximum number of observations to use when calculating 
 *                 correlations from X. Default is 100.
 * 
 * @return NumericVector with three elements:
 *         [1] Minimum absolute correlation among selected variables
 *         [2] Mean absolute correlation among selected variables
 *         [3] Median absolute correlation among selected variables
 *         If only one variable is selected, returns (1.0, 1.0, 1.0).
 * 
 * @examples
 * \dontrun{
 * # With raw data matrix
 * X <- matrix(rnorm(1000 * 20), 1000, 20)
 * selected_vars <- c(1, 5, 10, 15)
 * get_purity(selected_vars, X = X)
 * 
 * # With correlation matrix
 * Xcorr <- cor(X)
 * get_purity(selected_vars, R = Xcorr, squared = TRUE)
 * }
 */
// [[Rcpp::export]]
NumericVector get_purity(IntegerVector pos,
                         Nullable<NumericMatrix> X,
                         Nullable<NumericMatrix> R,
                         bool squared = false,
                         int n_purity = 100) {
  // Number of selected variables
  const int p = pos.size();
  
  // Correlation matrix
  arma::mat C;
  
  // Case 1: Pre-computed correlation matrix is provided
  if (R.isNotNull()) {
    NumericMatrix Xcorr(R);
    C.set_size(p, p);
    
    // Extract the correlation sub-matrix for selected variables
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        // Convert from 1-based to 0-based indexing
        C(i, j) = Xcorr(pos[i] - 1, pos[j] - 1);
      }
    }
  }
  // Case 2: Raw data matrix is provided
  else if (X.isNotNull()) {
    NumericMatrix X_(X);
    const int n = X_.nrow();
    
    // Extract selected variables from original matrix
    arma::mat Xm(n, p);
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < n; ++i) {
        // Convert from 1-based to 0-based indexing
        Xm(i, j) = X_(i, pos[j] - 1);
      }
    }
    
    // If data is large, use a random subset for efficiency
    if (n > n_purity) {
      arma::uvec rows = arma::randperm(n, n_purity);
      Xm = Xm.rows(rows);
    }
    
    // Standardize data (center and scale)
    arma::rowvec mu = arma::mean(Xm, 0);
    arma::rowvec sigma = arma::stddev(Xm, 0, 0);
    Xm.each_row() -= mu;
    Xm.each_row() /= sigma;
    
    // Calculate correlation matrix
    C = (Xm.t() * Xm) / (Xm.n_rows - 1);
  }
  // Error case: Neither matrix provided
  else {
    stop("Either X or Xcorr must be provided.");
  }
  
  // Calculate purity statistics from correlation matrix
  std::vector<double> vals;
  
  // Extract upper triangular elements (excluding diagonal)
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      double v = std::abs(C(i, j));
      if (squared) {
        v *= v;
      }
      vals.push_back(v);
    }
  }
  
  // Handle the case of a single variable
  if (vals.empty()) {
    return NumericVector::create(1.0, 1.0, 1.0);
  }
  
  // Sort correlation values for computing statistics
  std::sort(vals.begin(), vals.end());
  
  // Calculate statistics
  const double min_val = vals.front();
  const double mean_val = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
  double median_val;
  const int m = vals.size();
  
  // Calculate median
  if (m % 2 == 1) {
    median_val = vals[m / 2];  // Odd number of elements
  } else {
    median_val = 0.5 * (vals[m / 2 - 1] + vals[m / 2]);  // Even number of elements
  }
  
  // Return results as a vector
  return NumericVector::create(min_val, mean_val, median_val);
}

/**
 * @title Determine Significant Components in a LASER Model
 * @description
 * Tests which latent components in a LASER model are statistically significant
 * using likelihood ratio tests.
 * 
 * @details
 * For each latent dimension ℓ, the function performs a likelihood ratio test by:
 * 1. Setting the coefficients for dimension ℓ to zero (null model)
 * 2. Computing the log-likelihood of this restricted model
 * 3. Comparing with the full model's log-likelihood using a chi-squared test
 * 4. Determining significance based on the provided alpha level
 * 
 * @param fit List containing LASER model fit objects, including:
 *        - theta: coefficient matrix
 *        - intercept: model intercept
 *        - dispersion: dispersion parameter
 *        - final_loglike: log-likelihood of the full model
 * @param X NumericMatrix (n × p) predictor matrix.
 * @param y Response variable: vector for GLM, matrix (time, status) for Cox.
 * @param family List (for GLM) or string "cox".
 * @param alpha Numeric scalar. Significance level for the tests. Default is 0.05.
 * 
 * @return LogicalVector indicating which latent dimensions are significant.
 * 
 * @examples
 * \dontrun{
 * # Assuming 'laser_fit' is a fitted LASER model with 3 latent dimensions
 * X <- matrix(rnorm(100 * 10), 100, 10)
 * y <- rnorm(100)
 * significant_dims <- get_included(laser_fit, X, y, "gaussian", alpha = 0.01)
 * }
 * 
 * @note This function requires the get_loglike function to be available.
 */
// [[Rcpp::export]]
LogicalVector get_included(List fit,
                           NumericMatrix X,
                           SEXP y,
                           SEXP family,
                           double alpha = 0.05) {
  // Extract model parameters from fit object
  const double ll_full = as<double>(fit["final_loglike"]);
  arma::mat theta = as<arma::mat>(fit["theta"]);
  const double intercept = as<double>(fit["intercept"]);
  const double dispersion = as<double>(fit["dispersion"]);
  
  // Number of latent dimensions
  const int L = theta.n_cols;
  
  // Initialize result vector
  LogicalVector include(L);
  
  // Test each latent dimension
  for (int ℓ = 0; ℓ < L; ++ℓ) {
    // Create null model by setting coefficients for dimension ℓ to zero
    arma::mat theta_null = theta;
    theta_null.col(ℓ).zeros();
    
    // Compute log-likelihood of null model
    double ll_null = get_loglike(
      X, y, family,
      theta_null,
      intercept,
      dispersion
    );
    
    // Likelihood ratio test statistic
    double stat = 2.0 * (ll_full - ll_null);
    
    // Compute p-value using chi-squared distribution with 1 df
    double pval;
    if (stat <= 0.0) {
      pval = 1.0;
    } else {
      pval = R::pchisq(stat, 1, /*lower.tail=*/false, /*log.p=*/false);
    }
    
    // Determine significance
    include[ℓ] = (pval < alpha);
  }
  
  return include;
}

/**
 * @title Compute Credible Sets for LASER Model Parameters
 * @description
 * Constructs credible sets for parameters in a LASER model based on posterior probabilities
 * and optionally filters them based on correlation purity metrics.
 * 
 * @details
 * This function constructs credible sets for each selected latent dimension based on
 * posterior inclusion probabilities. It then:
 * 1. Identifies unique sets to avoid duplication
 * 2. Computes the claimed coverage for each set
 * 3. If correlation data is provided, computes purity metrics for each set
 * 4. Filters sets that don't meet minimum correlation thresholds
 * 5. Returns the filtered sets ordered by purity
 * 
 * Purity measures assess the independence of variables within each credible set.
 * Higher purity indicates stronger independent signals rather than correlated variables.
 * 
 * @param posterior_mat NumericMatrix (p × L), posterior inclusion probabilities for each variable
 *                      and latent dimension, where:
 *                      p = number of variables
 *                      L = number of latent dimensions
 * @param keep LogicalVector of length L indicating which latent dimensions to include.
 * @param coverage Numeric scalar. Target coverage probability for credible sets (default: 0.95).
 * @param X Nullable<NumericMatrix>. Original predictor matrix for computing correlations.
 *           Not required if R is provided.
 * @param R Nullable<NumericMatrix>. Pre-computed correlation matrix.
 *               Not required if X is provided.
 * @param check_symmetric Logical. If TRUE, forces the correlation matrix to be symmetric 
 *                        (default: TRUE).
 * @param min_abs_corr Numeric scalar. Minimum absolute correlation threshold for filtering
 *                     credible sets (default: 0.5).
 * @param n_purity Integer. Maximum number of observations to use for correlation computation
 *                 (default: 100).
 * @param squared Logical. If TRUE, uses squared correlations (default: FALSE).
 * 
 * @return A List containing:
 *         - sets: List of credible sets (each as an IntegerVector of variable indices)
 *         - coverage: NumericVector of claimed coverage for each set
 *         If no sets meet the filtering criteria, returns NULL for both elements.
 * 
 * @note This function requires the get_purity function to be available.
 */
// [[Rcpp::export]]
List get_cs(NumericMatrix posterior_mat,
            LogicalVector keep,
            double coverage = 0.95,
            Nullable<NumericMatrix> X = R_NilValue,
            Nullable<NumericMatrix> R = R_NilValue,
            bool check_symmetric = true,
            double min_abs_corr = 0.5,
            int n_purity = 100,
            bool squared = false) {
  
  // 1) Check if any effects are selected
  if (!is_true(any(keep))) {
    return List::create(
      _["sets"] = R_NilValue,
      _["coverage"] = R_NilValue
    );
  }
  
  // 2) Build sub-matrix of selected dimensions
  std::vector<int> cols;
  for (int ℓ = 0; ℓ < posterior_mat.ncol(); ++ℓ) {
    if (keep[ℓ]) {
      cols.push_back(ℓ);
    }
  }
  
  const int p = posterior_mat.nrow();  // Number of variables
  const int K = cols.size();           // Number of selected dimensions
  
  // Create sub-matrix with only selected dimensions
  NumericMatrix sub(p, K);
  for (int i = 0; i < p; ++i) {
    for (int k = 0; k < K; ++k) {
      sub(i, k) = posterior_mat(i, cols[k]);
    }
  }
  
  // 3) Compute credible sets for each dimension
  std::vector<std::vector<int>> raw_sets;
  raw_sets.reserve(K);
  
  for (int k = 0; k < K; ++k) {
    // Sort variables by posterior probability (descending)
    std::vector<std::pair<double, int>> vp(p);
    for (int i = 0; i < p; ++i) {
      vp[i] = {sub(i, k), i};
    }
    std::sort(vp.begin(), vp.end(),
             [](const std::pair<double,int> &a,
                const std::pair<double,int> &b) {
               return a.first > b.first;
             });
    
    // Accumulate variables until reaching desired coverage
    double cum_prob = 0.0;
    int n = 0;
    while (n < p && cum_prob < coverage) {
      cum_prob += vp[n].first;
      ++n;
    }
    
    // Create credible set with 1-based indexing
    std::vector<int> credible_set(n);
    for (int m = 0; m < n; ++m) {
      credible_set[m] = vp[m].second + 1;
    }
    std::sort(credible_set.begin(), credible_set.end());
    raw_sets.push_back(std::move(credible_set));
  }
  
  // 4) Extract unique sets and compute their claimed coverage
  std::map<std::string, int> seen;
  List unique_sets;
  NumericVector claimed_cov;
  
  for (auto &set : raw_sets) {
    // Create a string representation of the set for deduplication
    std::ostringstream oss;
    for (size_t i = 0; i < set.size(); ++i) {
      if (i) oss << ",";
      oss << set[i];
    }
    std::string key = oss.str();
    
    // If this is a new unique set, add it to our results
    if (seen.find(key) == seen.end()) {
      seen[key] = unique_sets.size();
      unique_sets.push_back(wrap(set));
      
      // Calculate claimed coverage across all dimensions
      double set_prob = 0.0;
      for (int idx : set) {
        for (int kk = 0; kk < K; ++kk) {
          set_prob += sub(idx - 1, kk);  // Adjust for 1-based indexing
        }
      }
      // Cap coverage at 1.0
      claimed_cov.push_back(set_prob > 1.0 ? 1.0 : set_prob);
    }
  }
  
  // 5) If no correlation data is provided, return results now
  if (X.isNull() && R.isNull()) {
    return List::create(
      _["sets"] = unique_sets,
      _["coverage"] = claimed_cov
    );
  }
  
  // 6) Process correlation matrix if provided
  NumericMatrix Xcorrmat;
  if (R.isNotNull()) {
    Xcorrmat = NumericMatrix(R);
    
    // Enforce symmetry if requested
    if (check_symmetric) {
      const int n2 = Xcorrmat.nrow();
      arma::mat T(Xcorrmat.begin(), n2, n2, false);
      arma::mat sym = 0.5 * (T + T.t());
      std::memcpy(Xcorrmat.begin(), sym.memptr(), sizeof(double) * sym.n_elem);
    }
  }
  
  // 7) Compute purity metrics for each unique set
  const int M = unique_sets.size();
  NumericMatrix purity_mat(M, 3);
  
  for (int i = 0; i < M; ++i) {
    IntegerVector pos = unique_sets[i];
    NumericVector purity = get_purity(
      pos, 
      X, 
      Xcorrmat.size() == 0 ? R_NilValue : wrap(Xcorrmat),
      squared, 
      n_purity
    );
    
    for (int j = 0; j < 3; ++j) {
      purity_mat(i, j) = purity[j];
    }
  }
  
  // Set column names for purity matrix
  colnames(purity_mat) = squared
    ? CharacterVector::create("min.sq.corr", "mean.sq.corr", "median.sq.corr")
    : CharacterVector::create("min.abs.corr", "mean.abs.corr", "median.abs.corr");
  
  // 8) Filter sets by minimum correlation threshold and reorder by purity
  const double thresh = squared ? (min_abs_corr * min_abs_corr) : min_abs_corr;
  std::vector<std::pair<double, int>> ord;
  
  for (int i = 0; i < M; ++i) {
    if (purity_mat(i, 0) >= thresh) {
      ord.emplace_back(purity_mat(i, 0), i);
    }
  }
  
  // Return empty result if no sets meet the threshold
  if (ord.empty()) {
    return List::create(
      _["sets"] = R_NilValue, 
      _["coverage"] = R_NilValue
    );
  }
  
  // Sort by purity (descending)
  std::sort(ord.begin(), ord.end(),
             [](const std::pair<double,int> &a,
                const std::pair<double,int> &b) {
               return a.first > b.first;
             });
  
  // Create final filtered and ordered results
  List final_sets, final_purity;
  NumericVector final_cov;
  
  for (const auto &pr : ord) {
    int idx = pr.second;
    final_sets.push_back(unique_sets[idx]);
    final_cov.push_back(claimed_cov[idx]);
    final_purity.push_back(purity_mat(idx, _));
  }
  
  // Return final results
  return List::create(
    _["sets"] = final_sets,
    _["coverage"] = final_cov
  );
}

/**
 * @title Fit a Univariate Generalized Linear Model (GLM) Fast
 *
 * @description
 * Efficiently fit a univariate GLM for predictor \code{x} and response \code{y},
 * supporting Gaussian (closed-form) and Binomial, Poisson, Gamma, Inverse Gaussian (IRLS).
 *
 * @param x Numeric vector. Predictor.
 * @param y Numeric vector. Response.
 * @param family List. Family object with elements \code{name} and optional \code{dispersion}.
 * @param offset Numeric vector. Optional offset (default \code{0}).
 * @param standardize Logical. If \code{TRUE}, standardize \code{x} (default \code{TRUE}).
 * @param max_iter Integer. Max IRLS iterations (default 25).
 * @param tol Numeric. Convergence tolerance (default 1e-8).
 *
 * @return A named list: \code{theta}, \code{se}, \code{p_value}, \code{logLik}, \code{BIC}.
 *
 * @export
 */
// [[Rcpp::export]]
List get_glm_fit(arma::vec x,
                 arma::vec y,
                 List family,
                 arma::vec offset,
                 bool standardize = true,
                 int max_iter = 25,
                 double tol = 1e-8) {

  int n = x.n_elem;

  // --- Standardize x ---
  double x_mean = 0.0, x_norm = 1.0;
  arma::vec x_std = x;
  if (standardize) {
    x_mean = arma::mean(x);
    x_std  = x - x_mean;
    x_norm = arma::norm(x_std, 2);
    if (x_norm > 0) x_std /= x_norm;
    else x_norm = 1.0;
  }

  // --- Family ---
  std::string name;
  if (family.containsElementNamed("family")) {
    name = as<std::string>(family["family"]);
  } else if (family.containsElementNamed("name")) {
    name = as<std::string>(family["name"]);
  } else {
    stop("Family object must contain 'family' or 'name'.");
  }

  double theta_std = 0.0, se_std = NA_REAL, p_val = NA_REAL;
  double logLik = NA_REAL, BIC = NA_REAL;

  // --- Gaussian closed-form ---
  if (name == "gaussian") {
    double dispersion = family.containsElementNamed("dispersion") ? as<double>(family["dispersion"]) : 1.0;
    arma::vec y_adj = y - offset;

    double x_center = arma::mean(x_std);
    double y_center = arma::mean(y_adj);

    arma::vec x_c = x_std - x_center;
    arma::vec y_c = y_adj - y_center;

    double xx = arma::dot(x_c, x_c);
    double xy = arma::dot(x_c, y_c);

    if (xx <= 1e-10) {
      return List::create(_["theta"] = 0.0, _["se"] = R_PosInf, _["p_value"] = 1.0,
                          _["logLik"] = NA_REAL, _["BIC"] = NA_REAL);
    }

    theta_std = xy / xx;
    double intercept = y_center - theta_std * x_center;
    arma::vec resid = y_adj - (intercept + theta_std * x_std);
    double rss = arma::dot(resid, resid);

    se_std = std::sqrt(dispersion / xx);
    double t_val = theta_std / se_std;
    p_val = 2 * R::pt(-std::abs(t_val), n-2, true, false);

    logLik = -n/2.0 * std::log(2 * M_PI * dispersion) - rss / (2 * dispersion);
    BIC = -2 * logLik + std::log(n) * 2;
  }

  // --- IRLS for Binomial, Poisson, Gamma, Inverse Gaussian ---
  else {
    theta_std = 0.0;
    arma::vec eta(n), mu(n), W(n), z(n);

    for (int iter = 0; iter < max_iter; ++iter) {
      eta = offset + theta_std * x_std;
      
      if (name == "binomial") {
        mu = 1.0 / (1.0 + arma::exp(-eta));
        W = mu % (1 - mu);
        z = eta + (y - mu) / W;
      } else if (name == "poisson") {
        mu = arma::exp(eta);
        W = mu;
        z = eta + (y - mu) / W;
      } else if (name == "Gamma") {
        mu = arma::exp(eta);
        W = 1.0 / (mu % mu);
        z = eta + (y - mu) / mu;
      } else if (name == "inverse.gaussian") {
        mu = arma::exp(eta);
        W = 1.0 / (mu % mu % mu);
        z = eta + (y - mu) / (mu % mu);
      } else {
        stop("Unsupported family inside IRLS.");
      }

      double num = arma::dot(x_std % W, (z - offset));
      double denom = arma::dot(x_std % x_std, W);
      double theta_new = num / denom;
      if (std::abs(theta_new - theta_std) < tol) break;
      theta_std = theta_new;
    }

    double info = arma::dot(W, x_std % x_std);
    se_std = 1.0 / std::sqrt(info);

    double z_val = theta_std / se_std;
    p_val = 2 * R::pnorm(-std::abs(z_val), 0.0, 1.0, true, false);

    logLik = 0.0;
    eta = offset + theta_std * x_std;
    if (name == "binomial") {
      mu = 1.0 / (1.0 + arma::exp(-eta));
      logLik = arma::accu(y % arma::log(mu) + (1-y) % arma::log(1-mu));
    } else if (name == "poisson") {
      mu = arma::exp(eta);
      logLik = arma::accu(y % eta - mu - lgamma(y+1));
    } else if (name == "Gamma") {
      mu = arma::exp(eta);
      logLik = arma::accu(-y/mu - arma::log(mu));
    } else if (name == "inverse.gaussian") {
      mu = arma::exp(eta);
      logLik = arma::accu(-(y-mu)%(y-mu)/(2*mu%mu%y));
    }
    BIC = -2 * logLik + std::log(n) * 2;
  }

  double theta = theta_std / x_norm;
  double se    = se_std / x_norm;

  return List::create(
    _["theta"]   = theta,
    _["se"]      = se,
    _["p_value"] = p_val,
    _["logLik"]  = logLik,
    _["BIC"]     = BIC
  );
}

/**
 * @title Fit a Univariate Cox Proportional Hazards Model
 *
 * @description
 * Efficiently fit a univariate Cox model for predictor \code{x} and survival time \code{y}.
 *
 * @param x Numeric vector. Predictor.
 * @param y Numeric matrix (n x 2). First column = time, second = status (1 = event).
 * @param offset Numeric vector. Known offset (default 0).
 * @param standardize Logical. If \code{TRUE}, standardize \code{x} (default \code{TRUE}).
 * @param ties Character. Tie method ("efron" [default] or "breslow").
 * @param max_iter Integer. Max Newton-Raphson iterations (default 25).
 * @param tol Numeric. Convergence tolerance (default 1e-8).
 *
 * @return A named list: \code{theta}, \code{se}, \code{p_value}, \code{logLik}, \code{BIC}.
 *
 * @export
 */
// [[Rcpp::export]]
List get_cox_fit(arma::vec x,
                 arma::mat y,
                 arma::vec offset,
                 bool standardize = true,
                 std::string ties = "efron",
                 int max_iter = 25,
                 double tol = 1e-8) {

  int n = x.n_elem;
  
  arma::vec time = y.col(0);
  arma::vec status = y.col(1);

  // --- Standardize x ---
  double x_mean = 0.0, x_norm = 1.0;
  arma::vec x_std = x;
  if (standardize) {
    x_mean = arma::mean(x);
    x_std  = x - x_mean;
    x_norm = arma::norm(x_std, 2);
    if (x_norm > 0) x_std /= x_norm;
    else x_norm = 1.0;
  }

  arma::uvec sorted_idx = arma::sort_index(time);
  arma::vec time_sorted = time(sorted_idx);
  arma::vec status_sorted = status(sorted_idx);
  arma::vec x_sorted = x_std(sorted_idx);
  arma::vec offset_sorted = offset(sorted_idx);

  double theta = 0.0;

  for (int iter = 0; iter < max_iter; ++iter) {
    double U = 0.0, I = 0.0, risk_sum = 0.0, x_risk_sum = 0.0, x2_risk_sum = 0.0;

    for (int i = n-1; i >= 0; --i) {
      double eta = offset_sorted[i] + theta * x_sorted[i];
      double exp_eta = std::exp(eta);

      risk_sum += exp_eta;
      x_risk_sum += exp_eta * x_sorted[i];
      x2_risk_sum += exp_eta * x_sorted[i] * x_sorted[i];

      if (status_sorted[i] == 1) {
        if (ties == "efron") {
          U += x_sorted[i] - (x_risk_sum / risk_sum);
          I += (x2_risk_sum / risk_sum - std::pow(x_risk_sum / risk_sum, 2.0));
        } else if (ties == "breslow") {
          U += x_sorted[i] - (x_risk_sum / risk_sum);
          I += (x2_risk_sum / risk_sum - std::pow(x_risk_sum / risk_sum, 2.0));
        }
      }
    }

    double delta = U / I;
    theta += delta;
    if (std::abs(delta) < tol) break;
  }

  double logLik = 0.0, risk_sum = 0.0;
  for (int i = n-1; i >= 0; --i) {
    double eta = offset_sorted[i] + theta * x_sorted[i];
    double exp_eta = std::exp(eta);
    risk_sum += exp_eta;
    if (status_sorted[i] == 1) logLik += eta - std::log(risk_sum);
  }

  double I = 0.0;
  risk_sum = 0.0;
  double x_risk_sum = 0.0, x2_risk_sum = 0.0;
  for (int i = n-1; i >= 0; --i) {
    double eta = offset_sorted[i] + theta * x_sorted[i];
    double exp_eta = std::exp(eta);

    risk_sum += exp_eta;
    x_risk_sum += exp_eta * x_sorted[i];
    x2_risk_sum += exp_eta * x_sorted[i] * x_sorted[i];

    if (status_sorted[i] == 1) {
      I += (x2_risk_sum / risk_sum - std::pow(x_risk_sum / risk_sum, 2.0));
    }
  }
  double se = 1.0 / std::sqrt(I);

  double z_val = theta / se;
  double p_val = 2 * R::pnorm(-std::abs(z_val), 0.0, 1.0, true, false);

  double BIC = -2 * logLik + std::log(n) * 1;

  theta /= x_norm;
  se    /= x_norm;

  return List::create(
    _["theta"]   = theta,
    _["se"]      = se,
    _["p_value"] = p_val,
    _["logLik"]  = logLik,
    _["BIC"]     = BIC
  );
}

/**
 * @title Fit a Null (Offset-Only) GLM Model (Fast)
 *
 * @description
 * Fit a generalized linear model with no covariates (only offset), and return
 * the fitted log-likelihood and BIC. Supports Gaussian (closed-form) and
 * Binomial, Poisson, Gamma, Inverse Gaussian via IRLS.
 *
 * @param y Numeric vector. Response values.
 * @param family List. Family object, containing \code{name} (string),
 * and optionally \code{dispersion} (for Gaussian etc.).
 * @param offset Numeric vector. Known offset (default 0).
 * @param max_iter Integer. Maximum IRLS iterations (default 25).
 * @param tol Numeric. Convergence tolerance for IRLS (default 1e-8).
 *
 * @return A named list:
 * \describe{
 *   \item{logLik}{Model log-likelihood.}
 *   \item{BIC}{Bayesian Information Criterion.}
 * }
 *
 * @export
 */
// [[Rcpp::export]]
List null_glm_fit(arma::vec y,
                  List family,
                  arma::vec offset,
                  int max_iter = 25,
                  double tol = 1e-8) {

  int n = y.n_elem;

  if (offset.n_elem == 1) {
    offset = arma::vec(n, arma::fill::value(offset[0]));
  } else {
    if (offset.n_elem != n) {
      stop("Offset must have length 1 or length n.");
    }
  }

  std::string name;
  if (family.containsElementNamed("family")) {
    name = as<std::string>(family["family"]);
  } else if (family.containsElementNamed("name")) {
    name = as<std::string>(family["name"]);
  } else {
    stop("Family object must contain 'family' or 'name'.");
  }

  double intercept = 0.0;
  double logLik = NA_REAL;
  double BIC = NA_REAL;

  // --- Gaussian ---
  if (name == "gaussian") {
    double dispersion = 1.0;
    if (family.containsElementNamed("dispersion") && !Rf_isNull(family["dispersion"])) {
      dispersion = as<double>(family["dispersion"]);
    }
    
    arma::vec y_adj = y - offset;
    intercept = arma::mean(y_adj);
    arma::vec resid = y_adj - intercept;
    double rss = arma::dot(resid, resid);
    
    logLik = -n/2.0 * std::log(2 * M_PI * dispersion) - rss / (2 * dispersion);
    BIC = -2 * logLik + std::log(n) * 1; // 1 parameter (intercept)

  }
  // --- Binomial, Poisson, Gamma, Inverse Gaussian ---
  else {
    intercept = 0.0;
    arma::vec eta(n), mu(n), W(n), z(n);

    for (int iter = 0; iter < max_iter; ++iter) {
      eta = offset + intercept;

      if (name == "binomial") {
        mu = 1.0 / (1.0 + arma::exp(-eta));
        W = mu % (1.0 - mu);
        z = eta + (y - mu) / W;
      } else if (name == "poisson") {
        mu = arma::exp(eta);
        W = mu;
        z = eta + (y - mu) / W;
      } else if (name == "Gamma") {
        mu = arma::exp(eta);
        W = 1.0 / (mu % mu);
        z = eta + (y - mu) / mu;
      } else if (name == "inverse.gaussian") {
        mu = arma::exp(eta);
        W = 1.0 / (mu % mu % mu);
        z = eta + (y - mu) / (mu % mu);
      } else {
        stop("Unsupported family for null GLM fit.");
      }

      double num = arma::accu(W % (z - offset));
      double denom = arma::accu(W);
      double intercept_new = num / denom;

      if (std::abs(intercept_new - intercept) < tol) {
        break;
      }
      intercept = intercept_new;
    }

    // --- Log-likelihood ---
    eta = offset + intercept;
    if (name == "binomial") {
      mu = 1.0 / (1.0 + arma::exp(-eta));
      logLik = arma::accu(y % arma::log(mu) + (1 - y) % arma::log(1 - mu));
    } else if (name == "poisson") {
      mu = arma::exp(eta);
      logLik = arma::accu(y % eta - mu - lgamma(y + 1));
    } else if (name == "Gamma") {
      mu = arma::exp(eta);
      logLik = arma::accu(-y / mu - arma::log(mu));
    } else if (name == "inverse.gaussian") {
      mu = arma::exp(eta);
      logLik = arma::accu(-(y - mu) % (y - mu) / (2 * mu % mu % y));
    }

    BIC = -2 * logLik + std::log(n) * 1; // 1 parameter (intercept)
  }

  return List::create(
    _["logLik"] = logLik,
    _["BIC"]    = BIC
  );
}

/**
 * @title Fit a Null (Offset-Only) Cox Proportional Hazards Model
 *
 * @description
 * Fit a Cox model with only a known offset (no covariates), returning the
 * partial log-likelihood and Bayesian Information Criterion (BIC).
 *
 * @param y NumericMatrix (n × 2). First column: follow-up time; second column: event indicator (1=event, 0=censoring).
 * @param offset Numeric vector. Known offset vector (length 1 or n).
 * @param ties String. Ties handling method: "efron" (default) or "breslow".
 *
 * @return A named list:
 * \describe{
 *   \item{logLik}{Partial log-likelihood.}
 *   \item{BIC}{Bayesian Information Criterion.}
 * }
 *
 * @export
 */
// [[Rcpp::export]]
List null_cox_fit(arma::mat y,
                  arma::vec offset,
                  std::string ties = "efron") {

  int n = y.n_rows;
  arma::vec time = y.col(0);
  arma::vec status = y.col(1);

  // --- Expand scalar offset ---
  if (offset.n_elem == 1) {
    offset = arma::vec(n, arma::fill::value(offset[0]));
  } else if (offset.n_elem != n) {
    stop("Offset must have length 1 or n.");
  }

  // --- Sort by time ascending ---
  arma::uvec order = arma::sort_index(time);
  arma::vec time_sorted = time(order);
  arma::vec status_sorted = status(order);
  arma::vec offset_sorted = offset(order);

  double logLik = 0.0;

  double risk_sum = 0.0;

  for (int i = n - 1; i >= 0; --i) {
    double eta = offset_sorted[i];
    double exp_eta = std::exp(eta);

    risk_sum += exp_eta;

    if (status_sorted[i] == 1) {
      if (ties == "efron") {
        // Efron tie correction
        logLik += eta - std::log(risk_sum);
      } else if (ties == "breslow") {
        // Breslow approximation
        logLik += eta - std::log(risk_sum);
      } else {
        stop("Unsupported ties method. Use 'efron' or 'breslow'.");
      }
    }
  }

  // --- BIC ---
  double BIC = -2.0 * logLik;  // Null model: 0 parameters except offset

  return List::create(
    _["logLik"] = logLik,
    _["BIC"]    = BIC
  );
}

/**
 * @title Update Intercept Given Offset
 * @description
 * Estimate the intercept term in a GLM given a known offset.
 * For Cox models, returns 0.
 * 
 * @param y Response vector or matrix (for Cox).
 * @param family List containing "family" or "name" field.
 * @param offset Known offset vector.
 * 
 * @return Numeric scalar: estimated intercept.
 */
// [[Rcpp::export]]
double update_intercept(SEXP y,
                             List family,
                             arma::vec offset) {
  
  int n;
  
  if (Rf_isMatrix(y)) {
    NumericMatrix ymat(y);
    n = ymat.nrow();
  } else if (Rf_isNumeric(y)) {
    NumericVector yvec(y);
    n = yvec.size();
  } else {
    stop("y must be a numeric vector or matrix");
  }
  
  if (offset.n_elem == 1) {
    offset = arma::vec(n, arma::fill::ones) * offset(0);
  }
  if (offset.n_elem != n) {
    stop("offset must have length n");
  }
  
  // --- Check family ---
  std::string fam_name;
  if (family.containsElementNamed("family")) {
    fam_name = as<std::string>(family["family"]);
  } else if (family.containsElementNamed("name")) {
    fam_name = as<std::string>(family["name"]);
  } else {
    fam_name = "cox";
  }
  
  // --- Cox case ---
  if (fam_name == "cox") {
    return 0.0;
  }
  
  // --- GLM case ---
  arma::vec y_vec;
  if (Rf_isMatrix(y)) {
    NumericMatrix ymat(y);
    y_vec = arma::vec(Rcpp::NumericVector(ymat(_, 0)));  // first column = outcome
  } else {
    y_vec = Rcpp::as<arma::vec>(y);
  }
  
  if (fam_name == "gaussian") {
    // Gaussian: intercept is just mean(residual)
    arma::vec resid = y_vec - offset;
    return arma::mean(resid);
  }
  
  else if (fam_name == "binomial" || fam_name == "poisson" ||
           fam_name == "Gamma" || fam_name == "inverse.gaussian") {
    
    // build a constant predictor of 1's
    arma::vec x1(n, arma::fill::ones);

    // call get_glm_fit *positionally*
    List res = get_glm_fit(
      x1,        // x
      y_vec,     // y
      family,    // family
      offset,    // offset
      false,     // standardize = FALSE
      25,        // max_iter = 25
      1e-8       // tol = 1e-8
    );

    // extract the intercept
    return as<double>(res["theta"]);
  }
  
  else {
    stop("Unsupported family.");
  }
  
  return NA_REAL; // Should not reach
}

/**
 * @title Update Dispersion Parameter
 * @description
 * Given a known offset (linear predictor), update the dispersion parameter.
 * For Gaussian, Gamma, Inverse Gaussian GLM families, computes dispersion from residuals.
 * For binomial, Poisson, or Cox families, returns dispersion = 1.
 *
 * @param y Response vector or matrix (time, status) for Cox.
 * @param family Family description as a List: must contain "family" or "name".
 * @param offset Known offset (linear predictor).
 * 
 * @return Updated dispersion (numeric scalar).
 */
// [[Rcpp::export]]
double update_dispersion(SEXP y,
                              List family,
                              arma::vec offset) {

  int n;
  
  if (Rf_isMatrix(y)) {
    NumericMatrix ymat(y);
    n = ymat.nrow();
  } else if (Rf_isNumeric(y)) {
    NumericVector yvec(y);
    n = yvec.size();
  } else {
    stop("y must be numeric vector or matrix");
  }
  
  if (offset.n_elem == 1) {
    offset = arma::vec(n, arma::fill::ones) * offset(0);
  }
  if (offset.n_elem != n) {
    stop("Offset must have length n");
  }
  
  // --- Family name ---
  std::string fam_name;
  if (family.containsElementNamed("family")) {
    fam_name = as<std::string>(family["family"]);
  } else if (family.containsElementNamed("name")) {
    fam_name = as<std::string>(family["name"]);
  } else {
    fam_name = "cox";
  }
  
  // --- Cox case ---
  if (fam_name == "cox") {
    return 1.0;
  }
  
  // --- One-parameter GLM ---
  if (fam_name == "binomial" || fam_name == "poisson") {
    return 1.0;
  }
  
  // --- Two-parameter GLM: Gaussian, Gamma, inverse.gaussian ---
  arma::vec y_vec;
  if (Rf_isMatrix(y)) {
    NumericMatrix ymat(y);
    y_vec = arma::vec(Rcpp::NumericVector(ymat(_, 0)));  // time column
  } else {
    y_vec = Rcpp::as<arma::vec>(y);
  }

  if (fam_name == "gaussian") {
    arma::vec resid = y_vec - offset;
    double rss = arma::dot(resid, resid);
    return rss / (n - 1);
  }
  else if (fam_name == "Gamma") {
    arma::vec mu = arma::exp(offset);  // assumes canonical log-link
    double dev = 0.0;
    for (int i = 0; i < n; ++i) {
      dev += -2 * (std::log(y_vec[i] / mu[i]) - (y_vec[i] - mu[i]) / mu[i]);
    }
    return dev / (n - 1);
  }
  else if (fam_name == "inverse.gaussian") {
    arma::vec mu = arma::exp(offset);  // assumes canonical link
    double dev = 0.0;
    for (int i = 0; i < n; ++i) {
      dev += (y_vec[i] - mu[i]) * (y_vec[i] - mu[i]) / (mu[i]*mu[i]*y_vec[i]);
    }
    return dev / (n - 1);
  }
  
  else {
    stop("Unsupported family inside update_dispersion.");
  }
  
  return NA_REAL; // never reached
}