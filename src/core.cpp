#include "core.h"
#include <Rcpp.h>
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double update_dispersion(const arma::vec& y,
                         SEXP family,
                         arma::vec offset,
                         std::string approach = "pearson") {
  
  int n = y.n_elem;
  int p = 0;

  // Expand offset if scalar or missing
  if (offset.n_elem == 0) {
    offset = arma::zeros(n);
  } else if (offset.n_elem == 1 && n > 1) {
    offset = arma::vec(n).fill(offset(0));
  } else if (offset.n_elem != n) {
    stop("offset must have length 1 or same length as y");
  }

  // Validate family object
  if (TYPEOF(family) != VECSXP) {
    stop("family must be a list (i.e., GLM family)");
  }

  List fam = as<List>(family);
  Function linkinv = fam["linkinv"];
  NumericVector mu_r = linkinv(wrap(offset));
  arma::vec mu = as<arma::vec>(mu_r);

  double dispersion = NA_REAL;

  if (approach == "pearson") {
    Function variance = fam["variance"];
    NumericVector var_r = variance(mu_r);
    arma::vec var_mu = as<arma::vec>(var_r);

    arma::vec residuals = (y - mu) / arma::sqrt(var_mu + 1e-8);
    dispersion = arma::accu(residuals % residuals) / (n - p);

  } else if (approach == "deviance") {
    Function dev_resids = fam["dev.resids"];
    NumericVector dev = dev_resids(wrap(y), wrap(mu), NumericVector(n, 1.0));
    dispersion = std::accumulate(dev.begin(), dev.end(), 0.0) / (n - p);

  } else {
    stop("approach must be 'pearson' or 'deviance'");
  }

  return dispersion;
}

// [[Rcpp::export]]
double univariate_loglik_cox(
    const arma::vec& x,
    const arma::mat& y,
    double theta,
    arma::vec offset,
    std::string ties = "efron"
) {
  // Validate dimensions
  int n = x.n_elem;
  if (y.n_rows != n || y.n_cols != 2) {
    Rcpp::stop("y must be an n x 2 matrix (time, status)");
  }

  // Extract time and status
  arma::vec time   = y.col(0);
  arma::vec status = y.col(1);

  // Prepare offset
  if (offset.n_elem == 1) {
    offset = arma::vec(n).fill(offset(0));
  } else if ((int)offset.n_elem != n) {
    Rcpp::stop("offset must be length 1 or length n");
  }

  // Linear predictor and exp(lp)
  arma::vec lp     = offset + theta * x;
  arma::vec exp_lp = arma::exp(lp);

  double logLik = 0.0;

  if (ties == "breslow") {
    // Breslow approximation: sum over each failure
    for (int i = 0; i < n; ++i) {
      if (status(i) == 1) {
        // risk set: times >= time[i]
        arma::uvec R = arma::find(time >= time(i));
        double denom = arma::sum(exp_lp.elem(R));
        logLik += lp(i) - std::log(denom);
      }
    }
  }
  else if (ties == "efron") {
    // Efron approximation: group failures by unique times
    std::set<double> utimes;
    for (int i = 0; i < n; ++i) if (status(i) == 1) utimes.insert(time(i));
    for (double t : utimes) {
      // indices of events at time t
      arma::uvec E = arma::find((time == t) && (status == 1));
      int d = E.n_elem;
      if (d == 0) continue;
      // sum of lp for events
      for (auto idx : E) logLik += lp(idx);
      // risk set at time t
      arma::uvec R = arma::find(time >= t);
      double Rsum  = arma::sum(exp_lp.elem(R));
      double Esum  = arma::sum(exp_lp.elem(E));
      // contribution from denominator
      for (int l = 0; l < d; ++l) {
        double adj = Rsum - (double)l / d * Esum;
        logLik -= std::log(adj);
      }
    }
  }
  else {
    Rcpp::stop("Unsupported ties method: '" + ties + "'. Use 'breslow' or 'efron'.");
  }

  return logLik;
}

// [[Rcpp::export]]
double univariate_loglik_glm(
    const arma::vec& x,
    const arma::vec& y,  
    SEXP family,
    double theta,
    const arma::vec& offset
) {
  int n = x.n_elem;
  
  // Check if y has correct dimensions
  if (y.n_elem != n) {
    stop("Length of 'y' must match length of 'x'");
  }
  
  // Check if offset is empty, if so, create a zero vector
  arma::vec offset_vec;
  if (offset.n_elem == 0) {
    offset_vec = arma::zeros(n);
  } else if (offset.n_elem != n) {
    stop("'offset' must have length equal to x");
  } else {
    offset_vec = offset;
  }
  
  // Convert family to List to properly access elements
  List fam(family);
  
  // Get family functions for link and deviance
  Function linkinv = fam["linkinv"];
  Function dev_resids = fam["dev.resids"];
  
  // Calculate linear predictor and fitted values
  arma::vec eta = offset_vec + theta * x;
  NumericVector eta_r = wrap(eta);
  NumericVector mu_r = linkinv(eta_r);
  
  // Use unit weights
  NumericVector w(n, 1.0);
  
  // Calculate deviance residuals
  NumericVector y_r = wrap(y);
  NumericVector dev = dev_resids(y_r, mu_r, w);
  
  // Convert to log-likelihood: loglik = -sum(dev)/2
  double sum_dev = sum(NumericVector(dev));
  
  // For families with estimated dispersion (gaussian, Gamma, inverse.gaussian),
  // we should profile out the dispersion parameter
  std::string family_name;
  try {
    family_name = as<std::string>(fam["family"]);
  } catch(...) {
    // If we can't extract the family name, default to fixed dispersion
    return -sum_dev / 2.0;
  }
  
  if (family_name == "gaussian") {
    // For Gaussian, profile out sigma^2
    double sigma2 = sum_dev / n;
    return -0.5 * n * log(2.0 * M_PI * sigma2) - 0.5 * n;
  } 
  else if (family_name == "Gamma") {
    // For Gamma, profile out the shape parameter
    double shape = n / sum_dev; 
    return n * shape * log(shape) - n * lgamma(shape) - n * shape;
  }
  else if (family_name == "inverse.gaussian") {
    // For inverse.gaussian, profile out lambda
    double lambda = n / sum_dev;
    return 0.5 * n * log(lambda / (2.0 * M_PI)) - 
           1.5 * sum(log(y_r)) - 0.5 * lambda * sum_dev;
  }
  else {
    // For families with fixed dispersion (binomial, poisson)
    return -sum_dev / 2.0;
  }
}


// [[Rcpp::export]]
double univariate_loglik(
    const arma::vec&   x,
    SEXP               y,        // vector for GLM or matrix for Cox
    SEXP               family,  
    double             theta,  
    const arma::vec&   offset,  
    std::string        ties     
) {
  int n = x.n_elem;
  // --- offset must be length n ---
  if ((int)offset.n_elem != n) {
    stop("`offset` must be length(x)");
  }
  
  // --- extract family name ---
  // Convert family to List first
  List fam(family);
  std::string famname = as<std::string>(fam["family"]);
  
  if (famname == "cox") {
    // Cox: y must be an n×2 matrix
    NumericMatrix Ym(y);
    arma::mat ymat = as<arma::mat>(Ym);
    if ((int)ymat.n_rows != n || (int)ymat.n_cols != 2) {
      stop("For Cox, y must be an n×2 matrix (time,status)");
    }
    return univariate_loglik_cox(x, ymat, theta, offset, ties);
  }
  else {
    NumericVector y_r(y);
    arma::vec y_arma = as<arma::vec>(y_r);
    return univariate_loglik_glm(x, y_arma, family, theta, offset);
  }
}

// [[Rcpp::export]]
double univariate_irls_cox(arma::vec x, 
                          arma::mat y,
                          arma::vec offset,
                          std::string ties = "efron",
                          double lambda = 0.0,       // Penalty strength parameter
                          double tau = 1e-5,         // Truncation parameter for TLP
                          int max_iter = 25,
                          double tol = 1e-8) {
  
  // Extract time and status from y matrix
  int n = x.n_elem;
  if (y.n_cols != 2 || y.n_rows != n) {
    stop("y must be a matrix with time and status columns and same number of rows as x");
  }
  
  arma::vec time = y.col(0);
  arma::vec status = y.col(1);
  
  // Initialize offset if not provided
  if (offset.n_elem == 0) {
    offset = arma::zeros(n);
  } else if (offset.n_elem == 1) {
    double offset_val = offset(0);
    offset = arma::ones(n) * offset_val;
  } else if (offset.n_elem != n) {
    stop("offset must have length 1 or same length as x");
  }
  
  // Check penalty parameters
  if (tau <= 0) {
    stop("tau must be positive");
  }
  
  // Sort data by time
  arma::uvec ord = arma::stable_sort_index(time);
  arma::vec time_sorted = time.elem(ord);
  arma::vec status_sorted = status.elem(ord);
  arma::vec x_sorted = x.elem(ord);
  arma::vec offset_sorted = offset.elem(ord);
  
  // IRLS loop for Cox regression with TLP
  double theta = 0.0;
  bool converged = false;
  double c = lambda / tau;
  
  for (int iter = 0; iter < max_iter; iter++) {
    // Linear predictor
    arma::vec eta = offset_sorted + theta * x_sorted;
    arma::vec exp_eta = arma::exp(eta);
    
    // Score and information calculation
    double score = 0.0;
    double info = 0.0;
    
    if (ties == "breslow") {
      // Breslow method for tied events (simpler)
      for (int i = 0; i < n; i++) {
        if (status_sorted(i) == 1) {  // If event occurred
          // Find risk set (subjects still at risk at time t_i)
          arma::uvec risk_set = arma::find(time_sorted >= time_sorted(i));
          arma::vec risk_x = x_sorted.elem(risk_set);
          arma::vec risk_exp_eta = exp_eta.elem(risk_set);
          
          // Sum of exp(eta) in risk set
          double sum_risk_exp_eta = arma::sum(risk_exp_eta);
          
          // Weighted average of x in risk set
          double weighted_x = arma::sum(risk_x % risk_exp_eta) / sum_risk_exp_eta;
          
          // Update score (first derivative)
          score += (x_sorted(i) - weighted_x);
          
          // Update information (negative second derivative)
          double weighted_x2 = arma::sum(arma::square(risk_x) % risk_exp_eta) / sum_risk_exp_eta;
          info += (weighted_x2 - weighted_x * weighted_x);
        }
      }
    } 
    else if (ties == "efron") {
      // Efron method for tied events (more accurate)
      
      // Find unique event times
      std::vector<double> unique_times;
      for (int i = 0; i < n; i++) {
        if (status_sorted(i) == 1) {
          // Check if this time is already in our vector
          bool found = false;
          for (size_t j = 0; j < unique_times.size(); j++) {
            if (std::abs(time_sorted(i) - unique_times[j]) < 1e-8) {
              found = true;
              break;
            }
          }
          if (!found) {
            unique_times.push_back(time_sorted(i));
          }
        }
      }
      
      // Process each unique event time
      for (size_t t_idx = 0; t_idx < unique_times.size(); t_idx++) {
        double t = unique_times[t_idx];
        
        // Find events at this time
        arma::uvec event_indices = arma::find((time_sorted == t) && (status_sorted == 1));
        int d = event_indices.n_elem;  // Number of tied events
        
        if (d == 0) continue;
        
        // Find risk set (all subjects at risk at time t)
        arma::uvec risk_set = arma::find(time_sorted >= t);
        arma::vec risk_x = x_sorted.elem(risk_set);
        arma::vec risk_exp_eta = exp_eta.elem(risk_set);
        
        // Event set variables
        arma::vec event_x = x_sorted.elem(event_indices);
        arma::vec event_exp_eta = exp_eta.elem(event_indices);
        
        // Sums for the risk set and event set
        double sum_risk_exp_eta = arma::sum(risk_exp_eta);
        double sum_event_exp_eta = arma::sum(event_exp_eta);
        arma::vec weighted_risk_x = risk_x % risk_exp_eta;
        double sum_weighted_risk_x = arma::sum(weighted_risk_x);
        arma::vec weighted_event_x = event_x % event_exp_eta;
        double sum_weighted_event_x = arma::sum(weighted_event_x);
        
        // Sum of x for all events at this time
        [[maybe_unused]] double sum_event_x = arma::sum(event_x);
        
        // Efron adjustment for score
        for (int j = 0; j < d; j++) {
          double fraction = j / static_cast<double>(d);
          double denom = sum_risk_exp_eta - fraction * sum_event_exp_eta;
          double weighted_mean = (sum_weighted_risk_x - fraction * sum_weighted_event_x) / denom;
          
          score += (event_x(j) - weighted_mean);
        }
        
        // Efron adjustment for information
        for (int j = 0; j < d; j++) {
          double fraction = j / static_cast<double>(d);
          double denom = sum_risk_exp_eta - fraction * sum_event_exp_eta;
          
          arma::vec sq_risk_x = arma::square(risk_x);
          arma::vec weighted_sq_risk_x = sq_risk_x % risk_exp_eta;
          double sum_weighted_sq_risk_x = arma::sum(weighted_sq_risk_x);
          
          arma::vec sq_event_x = arma::square(event_x);
          arma::vec weighted_sq_event_x = sq_event_x % event_exp_eta;
          double sum_weighted_sq_event_x = arma::sum(weighted_sq_event_x);
          
          double weighted_mean_x = (sum_weighted_risk_x - fraction * sum_weighted_event_x) / denom;
          double weighted_mean_x2 = (sum_weighted_sq_risk_x - fraction * sum_weighted_sq_event_x) / denom;
          
          info += (weighted_mean_x2 - std::pow(weighted_mean_x, 2));
        }
      }
    }
    else {
      stop("Unsupported ties method: '" + ties + "'. Use 'breslow' or 'efron'.");
    }
    
    // Check for numerical issues
    if (info <= 0 || !std::isfinite(score/info)) {
      Rcpp::warning("Non-positive or non-finite information matrix in Cox IRLS; returning current slope.");
      break;
    }
    
    // Add TLP gradient and hessian adjustment
    double theta_tlp = 0.0;
    
    if (lambda > 0) {
      // Apply TLP using closed-form solution
      // For Cox regression, the quadratic approximation gives us
      // score = b and info = a in the notation from univariate_glm_truncLasso
      
      double a = info;
      double b = score;
      
      // Ensure a is positive
      if (a <= 0) a = 1e-8;
      
      // Capped-ℓ₁ closed form solution
      if (std::abs(b) <= c) {
        theta_tlp = 0.0;
      } else if (std::abs(b) < a * tau + c) {
        theta_tlp = (b - c * sgn(b)) / a;
      } else {
        theta_tlp = b / a;
      }
    } else {
      // If no penalty (lambda = 0), use standard update
      theta_tlp = score / info;
    }
    
    // Apply the update
    double theta_new = theta + theta_tlp;
    
    // Convergence check
    if (std::abs(theta_new - theta) < tol) {
      theta = theta_new;
      converged = true;
      break;
    }
    
    theta = theta_new;
  }
  
  if (!converged) {
    Rcpp::warning("Cox IRLS algorithm with TLP did not converge in %d iterations", max_iter);
  }
  
  return theta;
}


// [[Rcpp::export]]
double univariate_irls_glm(const arma::vec&   x,
                           const arma::vec&   y,
                           SEXP               family,
                           arma::vec          offset,
                           double             lambda     = 0.0,
                           double             tau        = 1e-5,
                           int                max_iter   = 25,
                           double             tol        = 1e-8) {
  int n = x.n_elem;
  if ((int)y.n_elem != n)   stop("x and y must have same length");
  if (offset.n_elem == 0)   offset = arma::zeros<arma::vec>(n);
  if (offset.n_elem == 1 && n>1) offset = arma::vec(n).fill(offset(0));
  if ((int)offset.n_elem != n) stop("offset must be length 1 or n");
  if (tau <= 0)             stop("tau must be > 0");
  if (TYPEOF(family) != VECSXP) stop("family must be a stats::family object");


  // Extract family functions
  List fam = as<List>(family);
  Function linkinv = fam["linkinv"];
  Function varfun = fam["variance"];
  Function mu_eta = fam["mu.eta"];
  double dispersion = fam["dispersion"];
  CharacterVector family_name = fam["family"];
  bool is_poisson = (as<std::string>(family_name[0]) == "poisson");
  
  // Initialize
  arma::vec eta = offset;
  // Bound eta for numerical stability (especially for Poisson)
  if (is_poisson) {
    eta = arma::clamp(eta, -20.0, 20.0);  // Prevent overflow in exp()
  }
  
  NumericVector mu_r = linkinv(wrap(eta));
  arma::vec mu = as<arma::vec>(mu_r);
  
  // Ensure mu is valid (especially for Poisson where mu > 0)
  if (is_poisson) {
    for (int i = 0; i < n; i++) {
      if (mu[i] <= 0) mu[i] = 1e-8;
    }
  }
  
  double theta = 0.0, theta_new = 0.0;
  double c = lambda / tau;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    // Calculate IRLS weights & working response
    NumericVector var_r = varfun(wrap(mu));
    NumericVector gprime = mu_eta(wrap(eta));  // Using eta instead of mu_r for derivative
    
    // Handle numerical issues
    for (int i = 0; i < n; i++) {
      // For Poisson with log link, gprime = mu, and var_mu = mu
      // For small mu, we can get unstable values - set a floor
      if (R_IsNaN(var_r[i]) || var_r[i] <= 0) var_r[i] = 1e-8;
      if (R_IsNaN(gprime[i]) || gprime[i] <= 0) gprime[i] = 1e-8;
    }
    
    arma::vec var_mu = as<arma::vec>(var_r);
    arma::vec gprime_vec = as<arma::vec>(gprime);
    
    // Calculate weights (w = gprime² / (var_mu * dispersion))
    arma::vec w = arma::square(gprime_vec) / (var_mu * dispersion);
    
    // Handle weights that are too large or NaN
    for (int i = 0; i < n; i++) {
      if (!arma::is_finite(w[i]) || w[i] <= 0) w[i] = 1e-8;
      else if (w[i] > 1e8) w[i] = 1e8;  // Cap large weights
    }
    
    // Calculate working response
    arma::vec z(n);
    for (int i = 0; i < n; i++) {
      double diff = y[i] - mu[i];
      z[i] = eta[i] + diff / gprime_vec[i];
      // Handle extreme values
      if (!arma::is_finite(z[i])) z[i] = eta[i];
    }
    
    // Weighted LS subproblem
    arma::vec z0 = z - offset;
    double a = arma::dot(w % x, x);
    double b = arma::dot(w % x, z0);
    
    // Ensure a is positive
    if (a <= 0) a = 1e-8;
    
    // Capped-ℓ₁ closed form solution
    if (std::abs(b) <= c) {
      theta_new = 0.0;
    } else if (std::abs(b) < a * tau + c) {
      theta_new = (b - c * sgn(b)) / a;
    } else {
      theta_new = b / a;
    }
    
    // Check convergence
    if (std::abs(theta_new - theta) < tol) {
      theta = theta_new;
      break;
    }
    theta = theta_new;
    
    // Update eta and mu for next IRLS iteration
    eta = offset + theta * x;
    
    // Bound eta for numerical stability (especially for Poisson)
    if (is_poisson) {
      eta = arma::clamp(eta, -20.0, 20.0);  // Prevent overflow in exp()
    }
    
    mu_r = linkinv(wrap(eta));
    mu = as<arma::vec>(mu_r);
    
    // Ensure mu is valid (especially for Poisson where mu > 0)
    if (is_poisson) {
      for (int i = 0; i < n; i++) {
        if (mu[i] <= 0) mu[i] = 1e-8;
      }
    }
  }
  
  return theta;
}

// [[Rcpp::export]]
List univariate_fit(
    const arma::vec& x, 
    SEXP y, 
    SEXP family,
    arma::vec offset,
    bool standardize = true,
    std::string ties = "efron",
    double lambda = 0.0,
    double tau = 0.5,
    double null_threshold = 1e-6)
{
  // Get dimensions and extract family type
  int n = x.n_elem;
  bool is_cox = TYPEOF(family) == STRSXP ? as<std::string>(family) == "cox" : 
                 as<std::string>(as<List>(family)["family"]) == "cox";
  
  // Handle offset
  if (offset.n_elem == 1 && n > 1) offset = arma::vec(n, arma::fill::value(offset(0)));
  
  // Standardize predictor (only if non-constant)
  double x_mean = 0.0, x_norm = 1.0;
  arma::vec x_std = x;
  bool all_zero = arma::all(arma::abs(x) <= 1e-10);
  
  if (standardize && !all_zero) {
    x_mean = arma::mean(x);
    x_std = x - x_mean;
    x_norm = arma::norm(x_std, 2);
    if (x_norm > 1e-10) x_std /= x_norm;
    else { x_std = x; x_mean = 0.0; x_norm = 1.0; }
  }
  
  // Fit model using appropriate function
  double theta_std = 0.0;
  
  if (!all_zero) {
    if (is_cox) {
      // Convert y to appropriate format for Cox
      arma::mat y_mat = TYPEOF(y) == VECSXP ? 
        arma::mat(n, 2) : as<arma::mat>(y);
      
      if (TYPEOF(y) == VECSXP) {
        NumericVector time = as<NumericVector>(as<List>(y)[0]);
        NumericVector status = as<NumericVector>(as<List>(y)[1]);
        y_mat.col(0) = as<arma::vec>(time);
        y_mat.col(1) = as<arma::vec>(status);
      }
      
      theta_std = univariate_irls_cox(x_std, y_mat, offset, ties, lambda, tau);
    } else {
      // Convert y for GLM
      arma::vec y_vec = TYPEOF(y) == VECSXP ? 
        as<arma::vec>(as<NumericVector>(as<List>(y)[0])) : 
        as<arma::vec>(as<NumericVector>(y));
      
      theta_std = univariate_irls_glm(x_std, y_vec, family, offset, lambda, tau);
    }
  }
  
  // Transform coefficient, apply threshold, and compute metrics
  double theta = theta_std / x_norm;
  if (std::abs(theta) <= null_threshold) theta = 0.0;
  
  double loglik = univariate_loglik(x, y, family, theta, offset, ties);
  double bic = -2.0 * loglik + std::log(n) * (theta != 0.0 ? 2.0 : 0.0);
  
  return List::create(
    Named("theta") = theta,
    Named("loglik") = loglik,
    Named("bic") = bic
  );
}


// [[Rcpp::export]]
Rcpp::List single_effect_fit(
    const arma::mat&   X,
    SEXP               y,
    SEXP               family,
    arma::vec          offset,
    bool               standardize = true,
    std::string        ties = "efron",
    double             lambda = 0.0,
    double             tau = 0.5,
    double             null_threshold = 1e-6
) {
  // Get dimensions
  int n = X.n_rows;
  int p = X.n_cols;
  
  // Initialize result vectors
  arma::vec theta(p, arma::fill::zeros);
  arma::vec loglik(p, arma::fill::zeros);
  arma::vec bic(p, arma::fill::zeros);
  arma::vec bic_diff(p, arma::fill::zeros);
  
  // Expand offset if needed
  if (offset.n_elem == 1 && n > 1) {
    offset = arma::vec(n, arma::fill::value(offset(0)));
  }
  
  // Fit null model
  List res_null = univariate_fit(
    arma::zeros<arma::vec>(n),
    y,
    family,
    offset,
    standardize,
    ties,
    lambda,
    tau,
    null_threshold
  );
  
  double null_bic = as<double>(res_null["bic"]);
  
  // Fit univariate models for each predictor
  for (int j = 0; j < p; j++) {
    arma::vec x_j = X.col(j);
    
    List res = univariate_fit(
      x_j,
      y,
      family,
      offset,
      standardize,
      ties,
      lambda,
      tau,
      null_threshold
    );
    
    theta[j] = as<double>(res["theta"]);
    loglik[j] = as<double>(res["loglik"]);
    bic[j] = as<double>(res["bic"]);
    bic_diff[j] = bic[j] - null_bic;
  }
  
  // Shift BIC differences so minimum is 0
  double min_bic_diff = bic_diff.min();
  bic_diff = bic_diff - min_bic_diff;
  
  // Calculate Bayes factors and posterior model probabilities
  arma::vec bf = arma::exp(-0.5 * bic_diff);
  double sum_bf = arma::sum(bf);
  arma::vec pmp = bf / sum_bf;
  
  // Threshold small values
  for (int j = 0; j < p; j++) {
    if (pmp[j] <= null_threshold) {
      theta[j] = 0.0;
      pmp[j] = 0.0;
    }
  }
  
  // Calculate PMP-weighted expectations
  arma::vec expect_theta = pmp % theta;
  double mu1 = arma::dot(pmp, theta);
  double mu2 = arma::dot(pmp, theta % theta);
  double expect_variance = mu2 - mu1*mu1;

  
  // Threshold small expected values
  for (int j = 0; j < p; j++) {
    if (expect_theta[j] <= null_threshold) {
      expect_theta[j] = 0.0;
    }
  }
  
  if (expect_variance <= null_threshold) {
    expect_variance = 0.0;
  }
  
  // Return results as a list
  return List::create(
    Named("loglik") = loglik,
    Named("bic") = bic,
    Named("bic_diff") = bic_diff,
    Named("bf") = bf,
    Named("pmp") = pmp,
    Named("theta") = theta,
    Named("expect_theta") = expect_theta,
    Named("expect_variance") = expect_variance
  );
}

// [[Rcpp::export]]
List additive_effect_fit(
    const arma::mat& X,
    SEXP y,
    int L,
    SEXP family,
    bool standardize = true,
    std::string ties = "efron",
    double lambda = 0.0,
    double tau = 0.5,
    double null_threshold = 1e-6,
    double tol = 5e-2,
    double eps = 1e-5,
    int max_iter = 100)
{
  // Start timing using standard C++ time
  clock_t start_time = clock();
  
  // Get dimensions
  int n = X.n_rows;
  int p = X.n_cols;
  L = std::min(10, p);
  
  // Initialize family and check if we need to estimate dispersion
  List fam;
  bool is_cox = false;
  bool estimate_dispersion = false;
  
  if (TYPEOF(family) == STRSXP) {
    std::string fam_str = as<std::string>(family);
    is_cox = (fam_str == "cox");
    if (is_cox) {
      fam = List::create(
        Named("family") = "cox",
        Named("link") = "log",
        Named("dispersion") = 1.0
      );
    }
  } else if (TYPEOF(family) == VECSXP) {
    fam = as<List>(family);
    if (fam.containsElementNamed("family")) {
      std::string fam_str = as<std::string>(fam["family"]);
      is_cox = (fam_str == "cox");
      
      if (!is_cox && fam.containsElementNamed("dispersion")) {
        fam["dispersion"] = 1.0;
      }
      
      // Check if dispersion should be estimated
      if (!is_cox) {
        std::string fam_str = as<std::string>(fam["family"]);
        estimate_dispersion = (
          fam_str == "gaussian" || 
          fam_str == "Gamma" || 
          fam_str == "inverse.gaussian" || 
          fam_str == "quasibinomial" || 
          fam_str == "quasipoisson"
        );
      }
    }
  } else {
    Rcpp::stop("family must be a family object; e.g., gaussian() or binomial()");
  }
  
  // Initialize parameters
  double intercept = 0.0;
  arma::mat theta(p, L, arma::fill::zeros);
  
  // Initialize result matrices
  arma::mat loglik(p, L, arma::fill::zeros);
  arma::mat bic(p, L, arma::fill::zeros);
  arma::mat bic_diff(p, L, arma::fill::zeros);
  arma::mat bf(p, L, arma::fill::zeros);
  arma::mat pmp(p, L, arma::fill::zeros);
  arma::vec expect_variance(L, arma::fill::zeros);
  
  // Log-likelihood tracking
  arma::vec expect_loglik(max_iter, arma::fill::zeros);
  
  // Main iterations
  int iter;
  for (iter = 0; iter < max_iter; iter++) {
    // Calculate current linear predictor
    arma::vec linear_predictor(n, arma::fill::zeros);
    linear_predictor.fill(intercept);
    linear_predictor += X * theta * arma::ones<arma::vec>(L);
    
    // Update each single effect
    for (int l = 0; l < L; l++) {
      // Update offset by removing current effect
      arma::vec offset = linear_predictor - X * theta.col(l);
      
      // Fit single effect
      List res = single_effect_fit(
        X,
        y,
        fam,
        offset,
        standardize,
        ties,
        lambda,
        tau,
        null_threshold
      );
      
      // Extract results
      arma::vec res_theta = as<arma::vec>(res["theta"]);
      arma::vec res_loglik = as<arma::vec>(res["loglik"]);
      arma::vec res_bic = as<arma::vec>(res["bic"]);
      arma::vec res_bic_diff = as<arma::vec>(res["bic_diff"]);
      arma::vec res_bf = as<arma::vec>(res["bf"]);
      arma::vec res_pmp = as<arma::vec>(res["pmp"]);
      
      // Apply thresholding to expect_theta
      arma::vec res_expect_theta = as<arma::vec>(res["expect_theta"]);
      for (int j = 0; j < p; j++) {
        if (res_expect_theta(j) <= null_threshold) {
          res_expect_theta(j) = 0.0;
        }
      }
      
      // Update parameters
      theta.col(l) = res_expect_theta;
      loglik.col(l) = res_loglik;
      bic.col(l) = res_bic;
      bic_diff.col(l) = res_bic_diff;
      bf.col(l) = res_bf;
      pmp.col(l) = res_pmp;
      expect_variance(l) = as<double>(res["expect_variance"]);
      
      // Update linear predictor for next effect
      linear_predictor = offset + X * theta.col(l);
    }
    
    // Update intercept (GLM only)
    if (!is_cox) {
      arma::vec ones(n, arma::fill::ones);
      arma::vec current_offset = X * (theta * arma::ones<arma::vec>(L));
      
      // Convert y to vector for GLM
      arma::vec y_vec = TYPEOF(y) == VECSXP ? 
          as<arma::vec>(as<NumericVector>(as<List>(y)[0])) : 
          as<arma::vec>(as<NumericVector>(y));
      
      intercept = univariate_irls_glm(
        ones,
        y_vec,
        fam,
        current_offset,
        0.0,  // No penalty on intercept
        tau,
        100,
        1e-8
      );
    }
    
    // Update dispersion if needed
    if (estimate_dispersion) {
      arma::vec full_pred = arma::ones(n) * intercept + X * (theta * arma::ones<arma::vec>(L));
      
      // Convert y to vector for GLM
      arma::vec y_vec = TYPEOF(y) == VECSXP ? 
          as<arma::vec>(as<NumericVector>(as<List>(y)[0])) : 
          as<arma::vec>(as<NumericVector>(y));
      
      double new_dispersion = update_dispersion(
        y_vec,
        fam,
        full_pred,
        "pearson"
      );
      
      fam["dispersion"] = new_dispersion;
    }
    
    // Compute expected log-likelihood
    double ell = 0.0;
    for (int l = 0; l < L; l++) {
      for (int j = 0; j < p; j++) {
        ell += pmp(j, l) * loglik(j, l);
      }
    }
    expect_loglik(iter) = ell / (p * L);
    
    // Check convergence
    if (iter > 0 && std::abs(expect_loglik(iter) - expect_loglik(iter-1)) < tol) {
      break;
    }
  }

  int last;
  if (iter == max_iter) {
    last = max_iter - 1;
  } else {
    last = iter;
  }
  
  // Identify kept effects
  // arma::uvec kept = (expect_variance > tol);
  std::vector<bool> kept(expect_variance.n_elem);
  for (size_t i = 0; i < expect_variance.n_elem; ++i) {
    kept[i] = (expect_variance[i] < eps);
  }

  // Calculate elapsed time
  clock_t end_time = clock();
  double elapsed_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
  
  // Return results
  return List::create(
    Named("niter") = last + 1,
    Named("loglik") = loglik,
    Named("expect_loglik") = expect_loglik.subvec(0, last),
    Named("final_loglik") = expect_loglik(last),
    Named("intercept") = intercept,
    Named("dispersion") = as<double>(fam["dispersion"]),
    Named("theta") = theta,
    Named("pmp") = pmp,
    Named("bic") = bic,
    Named("bic_diff") = bic_diff,
    Named("bf") = bf,
    Named("expect_variance") = expect_variance,
    Named("kept") = kept,
    Named("elapsed_time") = elapsed_time
  );
}