#include "core.h"
#include <Rcpp.h>
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat decompose_theta(const arma::mat& theta, int L) {
  int p = theta.n_rows;

  // Step 1: compute row sums
  arma::vec row_sums = arma::sum(theta, 1);  // p × 1

  // Step 2: initialize new theta matrix
  arma::mat theta_new(p, L, arma::fill::zeros);

  // Step 3: redistribute row sums into L columns
  for (int i = 0; i < p; ++i) {
    int col = i % L;  // assign each row sum to one column (round-robin or diagonal-style)
    theta_new(i, col) = row_sums(i);
  }

  return theta_new;

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
    const arma::vec& offset,
    double intercept
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
  
  // Get family functions
  Function linkinv = fam["linkinv"];
  Function dev_resids = fam["dev.resids"];
  Function variance = fam["variance"];
  
  // Calculate linear predictor with intercept and fitted values
  arma::vec eta = intercept + offset_vec + theta * x;
  NumericVector eta_r = wrap(eta);
  NumericVector mu_r = linkinv(eta_r);
  
  // Use unit weights
  NumericVector w(n, 1.0);
  
  // Get the family name
  std::string family_name;
  try {
    family_name = as<std::string>(fam["family"]);
  } catch(...) {
    stop("Could not extract family name");
  }
  
  // Calculate log-likelihood based on family
  double loglik = 0.0;
  
  if (family_name == "gaussian") {
    // For Gaussian: -n/2*log(2*pi*sigma^2) - sum((y-mu)^2)/(2*sigma^2)
    NumericVector y_r = wrap(y);
    NumericVector resid = y_r - mu_r;
    double RSS = sum(resid * resid);
    double sigma2 = RSS / n;  // ML estimate of sigma^2
    loglik = -0.5 * n * log(2.0 * M_PI * sigma2) - 0.5 * RSS / sigma2;
  } 
  else if (family_name == "binomial") {
    // For binomial: sum(y*log(mu/(1-mu)) + log(1-mu))
    NumericVector y_r = wrap(y);
    loglik = 0.0;
    for (int i = 0; i < n; i++) {
      if (mu_r[i] <= 0.0 || mu_r[i] >= 1.0) {
        stop("Fitted probabilities numerically 0 or 1 occurred");
      }
      double p = mu_r[i];
      loglik += y_r[i] * log(p / (1.0 - p)) + log(1.0 - p);
    }
  }
  else if (family_name == "poisson") {
    // For Poisson: sum(y*log(mu) - mu - log(y!))
    NumericVector y_r = wrap(y);
    loglik = 0.0;
    for (int i = 0; i < n; i++) {
      if (mu_r[i] <= 0.0) {
        stop("Non-positive fitted means in Poisson family");
      }
      // Omit factorial term for constants
      loglik += y_r[i] * log(mu_r[i]) - mu_r[i];
      // Subtract log(y!) only for non-zero y
      if (y_r[i] > 0) {
        loglik -= lgamma(y_r[i] + 1.0);
      }
    }
  }
  else if (family_name == "Gamma") {
    // For Gamma with shape parameter alpha:
    // sum(alpha*log(alpha*y/mu) - alpha*y/mu - log(y) - lgamma(alpha))
    NumericVector y_r = wrap(y);
    NumericVector dev = dev_resids(y_r, mu_r, w);
    double sum_dev = sum(NumericVector(dev));
    double shape = n / sum_dev;  // ML estimate of shape parameter
    
    loglik = 0.0;
    for (int i = 0; i < n; i++) {
      if (y_r[i] <= 0.0 || mu_r[i] <= 0.0) {
        stop("Non-positive responses or fitted values in Gamma family");
      }
      loglik += shape * log(shape * y_r[i] / mu_r[i]) - 
                shape * y_r[i] / mu_r[i] - 
                log(y_r[i]) - 
                lgamma(shape);
    }
  }
  else if (family_name == "inverse.gaussian") {
    // For inverse.gaussian:
    // sum(0.5*log(lambda/(2*pi*y^3)) - 0.5*lambda*(y-mu)^2/(y*mu^2))
    NumericVector y_r = wrap(y);
    NumericVector dev = dev_resids(y_r, mu_r, w);
    double sum_dev = sum(NumericVector(dev));
    double lambda = n / sum_dev;  // ML estimate of dispersion parameter
    
    loglik = 0.0;
    for (int i = 0; i < n; i++) {
      if (y_r[i] <= 0.0 || mu_r[i] <= 0.0) {
        stop("Non-positive responses or fitted values in inverse.gaussian family");
      }
      loglik += 0.5 * log(lambda / (2.0 * M_PI * pow(y_r[i], 3))) - 
                0.5 * lambda * pow(y_r[i] - mu_r[i], 2) / (y_r[i] * pow(mu_r[i], 2));
    }
  }
  else if (family_name == "quasi" || family_name == "quasibinomial" || 
           family_name == "quasipoisson") {
    // For quasi-families, use deviance-based approximation
    NumericVector y_r = wrap(y);
    NumericVector dev = dev_resids(y_r, mu_r, w);
    double sum_dev = sum(NumericVector(dev));
    loglik = -sum_dev / 2.0;
  }
  else {
    // For any other family, use deviance-based approximation
    NumericVector y_r = wrap(y);
    NumericVector dev = dev_resids(y_r, mu_r, w);
    double sum_dev = sum(NumericVector(dev));
    loglik = -sum_dev / 2.0;
  }
  
  return loglik;
}


// [[Rcpp::export]]
double univariate_loglik(
    const arma::vec&   x,
    SEXP               y,        // vector for GLM or matrix for Cox
    SEXP               family,  
    double             theta,  
    const arma::vec&   offset,  
    double             intercept = 0.0, 
    std::string        ties = "efron"
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
    return univariate_loglik_glm(x, y_arma, family, theta, offset, intercept);
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
double univariate_irls_glm_no_intercept(const arma::vec&   x,
                                        const arma::vec&   y,
                                        SEXP               family,
                                        arma::vec          offset,
                                        int                max_iter   = 25,
                                        double             tol        = 1e-8) {
  int n = x.n_elem;
  if ((int)y.n_elem != n)   stop("x and y must have same length");
  if (offset.n_elem == 0)   offset = arma::zeros<arma::vec>(n);
  if (offset.n_elem == 1 && n > 1) offset = arma::vec(n).fill(offset(0));
  if ((int)offset.n_elem != n) stop("offset must be length 1 or n");
  if (TYPEOF(family) != VECSXP) stop("family must be a stats::family object");
  
  // Extract family functions
  List fam = as<List>(family);
  Function linkinv = fam["linkinv"];
  Function varfun = fam["variance"];
  Function mu_eta = fam["mu.eta"];
  
  // Handle dispersion - checking if NULL first
  double dispersion = 1.0;  // Default to 1.0 (appropriate for binomial/Poisson)

  // Check if "dispersion" exists in the family object
  if (fam.containsElementNamed("dispersion")) {
    // Additional check if it's not NULL
    if (!Rf_isNull(fam["dispersion"])) {
      dispersion = as<double>(fam["dispersion"]);
    }
  }
  
  CharacterVector family_name = fam["family"];
  bool is_poisson = (as<std::string>(family_name[0]) == "poisson");
  
  // Initialize slope parameter
  double theta = 0.0;
  
  // Initialize linear predictor
  arma::vec eta = offset;
  if (is_poisson) eta = arma::clamp(eta, -20.0, 20.0);
  
  NumericVector mu_r = linkinv(wrap(eta));
  arma::vec mu = as<arma::vec>(mu_r);
  
  if (is_poisson) {
    for (int i = 0; i < n; i++) if (mu[i] <= 0) mu[i] = 1e-8;
  }
  
  for (int iter = 0; iter < max_iter; ++iter) {
    NumericVector var_r = varfun(wrap(mu));
    NumericVector gprime = mu_eta(wrap(eta));
    
    for (int i = 0; i < n; i++) {
      if (R_IsNaN(var_r[i]) || var_r[i] <= 0) var_r[i] = 1e-8;
      if (R_IsNaN(gprime[i]) || gprime[i] <= 0) gprime[i] = 1e-8;
    }
    
    arma::vec var_mu = as<arma::vec>(var_r);
    arma::vec gprime_vec = as<arma::vec>(gprime);
    
    arma::vec w = arma::square(gprime_vec) / (var_mu * dispersion);
    
    for (int i = 0; i < n; i++) {
      if (!std::isfinite(w[i]) || w[i] <= 0) w[i] = 1e-8;
      else if (w[i] > 1e8) w[i] = 1e8;
    }
    
    arma::vec z(n);
    for (int i = 0; i < n; i++) {
      double diff = y[i] - mu[i];
      z[i] = eta[i] + diff / gprime_vec[i];
      if (!std::isfinite(z[i])) z[i] = eta[i];
    }
    
    arma::vec z0 = z - offset;
    arma::vec xw = x % sqrt(w);
    arma::vec zw = z0 % sqrt(w);
    
    // Handle potential numerical issues with XtWX
    double XtWX = dot(xw, xw);
    if (XtWX < 1e-8) XtWX = 1e-8;  // Prevent division by very small numbers
    
    double XtWz = dot(xw, zw);
    double theta_new = XtWz / XtWX;
    
    if (std::abs(theta_new - theta) < tol) {
      theta = theta_new;
      break;
    }
    
    theta = theta_new;
    
    eta = offset + theta * x;
    if (is_poisson) eta = arma::clamp(eta, -20.0, 20.0);
    
    mu_r = linkinv(wrap(eta));
    mu = as<arma::vec>(mu_r);
    
    if (is_poisson) {
      for (int i = 0; i < n; i++) if (mu[i] <= 0) mu[i] = 1e-8;
    }
  }
  
  return theta;
}

// [[Rcpp::export]]
Rcpp::List univariate_irls_glm(const arma::vec&   x,
                               const arma::vec&   y,
                               SEXP               family,
                               arma::vec          offset,
                               int                max_iter   = 25,
                               double             tol        = 1e-8) {
  int n = x.n_elem;
  if ((int)y.n_elem != n)   stop("x and y must have same length");
  if (offset.n_elem == 0)   offset = arma::zeros<arma::vec>(n);
  if (offset.n_elem == 1 && n>1) offset = arma::vec(n).fill(offset(0));
  if ((int)offset.n_elem != n) stop("offset must be length 1 or n");
  if (TYPEOF(family) != VECSXP) stop("family must be a stats::family object");

  // Extract family functions
  List fam = as<List>(family);
  Function linkinv = fam["linkinv"];
  Function varfun = fam["variance"];
  Function mu_eta = fam["mu.eta"];
  
  // Handle dispersion - checking if NULL first
  double dispersion = 1.0;  // Default to 1.0 (appropriate for binomial/Poisson)

  // Check if "dispersion" exists in the family object
  if (fam.containsElementNamed("dispersion")) {
    // Additional check if it's not NULL
    if (!Rf_isNull(fam["dispersion"])) {
      dispersion = as<double>(fam["dispersion"]);
    }
  }
  
  CharacterVector family_name = fam["family"];
  bool is_poisson = (as<std::string>(family_name[0]) == "poisson");
  
  // Initialize parameters
  double intercept = 0.0, intercept_new = 0.0;
  double theta = 0.0, theta_new = 0.0;
  
  // Initialize linear predictor
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
      if (!std::isfinite(w[i]) || w[i] <= 0) w[i] = 1e-8;
      else if (w[i] > 1e8) w[i] = 1e8;  // Cap large weights
    }
    
    // Calculate working response
    arma::vec z(n);
    for (int i = 0; i < n; i++) {
      double diff = y[i] - mu[i];
      z[i] = eta[i] + diff / gprime_vec[i];
      // Handle extreme values
      if (!std::isfinite(z[i])) z[i] = eta[i];
    }
    
    // Weighted LS subproblem for intercept + theta*x
    arma::vec z0 = z - offset;
    
    // Create design matrix with intercept (ones) and x
    arma::mat X(n, 2);
    X.col(0) = arma::ones<arma::vec>(n); // Intercept column
    X.col(1) = x;                        // x column
    
    // Create weighted X
    arma::mat Xw = X;
    for (int i = 0; i < n; i++) {
      Xw.row(i) *= std::sqrt(w[i]);
    }
    
    // Create weighted z
    arma::vec zw = z0 % arma::sqrt(w);
    
    // Normal equations: (X'WX)b = X'Wz
    arma::mat XtWX = Xw.t() * Xw;
    arma::vec XtWz = Xw.t() * zw;
    
    // Solve for intercept and theta
    arma::vec params;
    bool solved = arma::solve(params, XtWX, XtWz);
    
    if (!solved) {
      // Fallback to a more stable approach if standard solve fails
      arma::mat XtWX_reg = XtWX;
      XtWX_reg.diag() += 1e-6; // Add small regularization
      solved = arma::solve(params, XtWX_reg, XtWz);
      
      if (!solved) {
        // If still fails, use pseudoinverse
        params = arma::pinv(XtWX) * XtWz;
      }
    }
    
    intercept_new = params(0);
    theta_new = params(1);
    
    // Check convergence
    if (std::abs(intercept_new - intercept) < tol && 
        std::abs(theta_new - theta) < tol) {
      intercept = intercept_new;
      theta = theta_new;
      break;
    }
    
    intercept = intercept_new;
    theta = theta_new;
    
    // Update eta and mu for next IRLS iteration
    eta = offset + intercept + theta * x;
    
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
  
  // Return both parameters in a List
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("intercept") = intercept,
    Rcpp::Named("theta") = theta
  );
  
  return result;
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
    double tau = 0.5)
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
  
  if (standardize) {
    x_mean = arma::mean(x);
    x_std = x - x_mean;
    x_norm = arma::norm(x_std, 2);
    if (x_norm > 1e-10) x_std /= x_norm;
    else { x_std = x; x_mean = 0.0; x_norm = 1.0; }

  }
  
  // Fit model using appropriate function
  double intercept_std = 0.0;
  double theta_std = 0.0;
  
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
    
    Rcpp::List res = univariate_irls_glm(x_std, y_vec, family, offset);

    intercept_std = as<double>(res["intercept"]);
    theta_std     = as<double>(res["theta"]);

  }
  
  // Transform coefficient, apply threshold, and compute metrics
  double theta = theta_std / x_norm;
  double intercept = intercept_std - theta * x_mean;

  double loglik = univariate_loglik(x, y, family, theta, offset, intercept, ties);
  double loglik0 = univariate_loglik(x, y, family, 0.0, offset, intercept, ties);
  double lrt = 2.0 * (loglik - loglik0);
  double pval = (lrt > 0) ? R::pchisq(lrt, 1.0, false, false) : 1.0;
  double bic = -2.0 * (loglik - loglik0) + std::log(n);

  return List::create(
    Named("intercept") = intercept,
    Named("theta") = theta,
    Named("loglik") = loglik,
    Named("bic") = bic,
    Named("pval") = pval
  );
}

// [[Rcpp::export]]
List univariate_glm(const arma::vec& x,
                    const arma::vec& y,
                    SEXP family,
                    Nullable<NumericVector> offset) {
  int n = y.n_elem;
  if ((int)x.n_elem != n) stop("x and y must have the same length.");

  // Build data.frame (convert arma::vec -> NumericVector via wrap)
  bool has_off = offset.isNotNull();
  DataFrame df;
  if (has_off) {
    NumericVector off(offset);
    if (off.size() != n && off.size() != 1) stop("offset must be length n or 1.");
    if (off.size() == 1) off = NumericVector(n, off[0]);
    df = DataFrame::create(_["y"]=wrap(y), _["x"]=wrap(x), _["off"]=off);
  } else {
    df = DataFrame::create(_["y"]=wrap(y), _["x"]=wrap(x));
  }

  // stats namespace & concrete methods
  Environment stats      = Environment::namespace_env("stats");
  Function     glm       = stats["glm"];
  Function     logLik    = stats["logLik"];
  Function     BICf      = stats["BIC"];
  Function     pchisq    = stats["pchisq"];
  Function     as_formula= stats["as.formula"];
  Function     summaryGLM= stats["summary.glm"]; // concrete method

  // Formulas
  SEXP f1 = has_off ? as_formula("y ~ x + offset(off)") : as_formula("y ~ x");
  SEXP f0 = has_off ? as_formula("y ~ 1 + offset(off)") : as_formula("y ~ 1");

  // Fits
  List fit1 = glm(_["formula"]=f1, _["data"]=df, _["family"]=family);
  List fit0 = glm(_["formula"]=f0, _["data"]=df, _["family"]=family);

  // Likelihoods & BICs
  double ll1  = as<double>(logLik(fit1));
  double ll0  = as<double>(logLik(fit0));
  double bic1 = as<double>(BICf(fit1));
  double bic0 = as<double>(BICf(fit0));

  // Coef/SE via summary.glm
  List s1 = summaryGLM(fit1);
  NumericMatrix coefs = s1["coefficients"];  // rows: (Intercept), x
  if (coefs.nrow() < 2) stop("summary.glm did not return expected coefficient matrix.");
  double intercept = coefs(0,0);
  double beta      = coefs(1,0);
  double se_beta   = coefs(1,1);
  double wald_p    = coefs(1,3);

  // LRT vs null
  double lrt   = 2.0 * (ll1 - ll0);
  double lrt_p = (lrt > 0.0) ? as<double>(pchisq(lrt, 1.0, false, false)) : 1.0;

  // BIC-based deltas / BF
  double deltaBIC = bic1 - bic0;   // negative favors variable model
  double twoLogBF = bic0 - bic1;   // ≈ 2*log BF_{M1 vs M0}

  return List::create(
    _["intercept"] = intercept,
    _["beta"]      = beta,
    _["se"]        = se_beta,
    _["wald_p"]    = wald_p,
    _["logLik1"]   = ll1,
    _["logLik0"]   = ll0,
    _["BIC1"]      = bic1,
    _["BIC0"]      = bic0,
    _["LRT"]       = lrt,
    _["LRT_p"]     = lrt_p,
    _["deltaBIC"]  = deltaBIC,
    _["twoLogBF"]  = twoLogBF
  );
}


// [[Rcpp::export]]
Rcpp::List univariate_cox(SEXP y,
                          const arma::vec& x,
                          Rcpp::Nullable<Rcpp::NumericVector> offset = R_NilValue,
                          std::string ties = "efron") {
  // ---- parse y = (time, status) ----
  NumericVector time, status;
  if (TYPEOF(y) == VECSXP) {
    List yl = y;
    if (yl.size() < 2) stop("y as list must contain (time, status).");
    time   = as<NumericVector>(yl[0]);
    status = as<NumericVector>(yl[1]);
  } else if (TYPEOF(y) == REALSXP || TYPEOF(y) == INTSXP) {
    NumericMatrix ym = as<NumericMatrix>(y);
    if (ym.ncol() < 2) stop("y as matrix must be n x 2: (time, status).");
    time   = ym(_, 0);
    status = ym(_, 1);
  } else {
    stop("y must be list(time, status) or an n x 2 numeric matrix.");
  }

  int n = time.size();
  if (status.size() != n || static_cast<int>(x.n_elem) != n)
    stop("Lengths of time, status, and x must match.");

  // ---- handle offset (optional; scalar expands) ----
  bool has_off = offset.isNotNull();
  NumericVector off;
  if (has_off) {
    off = NumericVector(offset);
    if (off.size() == 1) off = NumericVector(n, off[0]);
    if (off.size() != n) stop("offset must be length n or length 1.");
  }

  // ---- build data.frame ----
  // (wrap(x) converts arma::vec -> NumericVector)
  DataFrame df = has_off
    ? DataFrame::create(_["time"]=time, _["status"]=status, _["x"]=wrap(x), _["off"]=off)
    : DataFrame::create(_["time"]=time, _["status"]=status, _["x"]=wrap(x));

  // ---- namespaces & functions ----
  Environment survival = Environment::namespace_env("survival");
  Environment stats    = Environment::namespace_env("stats");
  Function coxph       = survival["coxph"];
  Function summ_coxph  = survival["summary.coxph"];   // concrete method
  Function as_formula  = stats["as.formula"];
  Function pchisq      = stats["pchisq"];

  // ---- formula: Surv(time, status) ~ x (+ offset(off)) ----
  SEXP f1 = has_off
    ? as_formula("Surv(time, status) ~ x + offset(off)")
    : as_formula("Surv(time, status) ~ x");

  // ---- fit model ----
  List fit1 = coxph(_["formula"]=f1, _["data"]=df, _["ties"]=ties);

  // ---- coef / se / z / p via summary.coxph ----
  List s1 = summ_coxph(fit1);
  NumericMatrix coefmat = s1["coefficients"]; // cols: coef, exp(coef), se(coef), z, p
  if (coefmat.nrow() < 1 || coefmat.ncol() < 5)
    stop("summary.coxph did not return expected coefficients matrix.");

  double beta    = coefmat(0, 0);
  double se_beta = coefmat(0, 2);
  double wald_z  = coefmat(0, 3);
  double wald_p  = coefmat(0, 4);

  // ---- log-likelihoods: fit1$loglik = c(ll0, ll1) ----
  NumericVector llv = fit1["loglik"];
  if (llv.size() < 2) stop("coxph did not return both null and fitted log-likelihoods.");
  double ll0 = llv[0];
  double ll1 = llv[1];

  // ---- LRT vs null ----
  double LRT   = 2.0 * (ll1 - ll0);
  double LRT_p = (LRT > 0.0) ? as<double>(pchisq(LRT, 1.0, false, false)) : 1.0;

  // ---- BIC (partial log-likelihood), k1=1, k0=0 ----
  double BIC1     = -2.0 * ll1 + std::log(static_cast<double>(n));
  double BIC0     = -2.0 * ll0;
  double deltaBIC = BIC1 - BIC0;       // negative favors model with x
  double twoLogBF = BIC0 - BIC1;       // ≈ 2*log BF_{x vs null}

  return List::create(
    _["beta"]       = beta,
    _["se"]         = se_beta,
    _["wald_z"]     = wald_z,
    _["wald_p"]     = wald_p,
    _["logLik0"]    = ll0,
    _["logLik1"]    = ll1,
    _["LRT"]        = LRT,
    _["LRT_p"]      = LRT_p,
    _["BIC0"]       = BIC0,
    _["BIC1"]       = BIC1,
    _["deltaBIC"]   = deltaBIC,
    _["twoLogBF"]   = twoLogBF,
    _["ties"]       = ties,
    _["has_offset"] = has_off
  );
}

// [[Rcpp::export]]
Rcpp::List single_effect_fit(
    const arma::mat&   X,
    SEXP               y,
    SEXP               family,
    arma::vec          offset,
    bool               standardize = true,
    bool               shrinkage = true,
    std::string        ties = "efron",
    double             lambda = 0.0,
    double             tau = 0.5,
    double             alpha = 0.05
) {
  // Get dimensions
  int n = X.n_rows;
  int p = X.n_cols;
  
  // Initialize result vectors
  arma::vec intercept(p, arma::fill::zeros);
  arma::vec theta(p, arma::fill::zeros);
  arma::vec loglik(p, arma::fill::zeros);
  arma::vec bic(p, arma::fill::zeros);
  arma::vec bic_diff(p, arma::fill::zeros);
  arma::vec pval_raw(p, arma::fill::zeros);
  arma::vec pval_intercept(p, arma::fill::zeros);
  arma::vec pval_theta(p, arma::fill::zeros);
  arma::vec evidence(p, arma::fill::zeros);
  arma::vec evidence_raw(p, arma::fill::zeros);

  // Expand offset if needed
  if (offset.n_elem == 1 && n > 1) {
    offset = arma::vec(n, arma::fill::value(offset(0)));
  }
  
  double null_bic = 0.0;
  
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
      tau
    );
    
    intercept[j] = as<double>(res["intercept"]);
    theta[j] = as<double>(res["theta"]);
    loglik[j] = as<double>(res["loglik"]);
    bic[j] = as<double>(res["bic"]);
    bic_diff[j] = bic[j] - null_bic;
    pval_raw[j] = as<double>(res["pval"]);
    evidence_raw[j] = -0.5*bic[j];
  }
  
  // Shift BIC differences so minimum is 0
  double min_bic_diff = bic_diff.min();
  bic_diff = bic_diff - min_bic_diff;
  
  // Calculate Bayes factors and posterior model probabilities
  arma::vec bf = arma::exp(-0.5 * bic_diff);
  double sum_bf = arma::sum(bf);
  arma::vec pmp = bf / sum_bf;
  
  // Calculate PMP-weighted expectations
  arma::vec expect_intercept = pmp % intercept;
  arma::vec expect_theta = pmp % theta;
  double mu1 = arma::dot(pmp, theta);
  double mu2 = arma::dot(pmp, theta % theta);
  double expect_variance = mu2 - mu1*mu1;
  
  for (int j = 0; j < p; j++) {
    arma::vec x_j = X.col(j);
    // Compute full model log-likelihood
    double ll1 = univariate_loglik(x_j, y, family, expect_theta[j], offset, expect_intercept[j]);

    // === Intercept Test ===
    // Null model: intercept = 0, theta fixed
    double ll0_intercept = univariate_loglik(x_j, y, family, expect_theta[j], offset, 0.0);
    double lrt_intercept = 2.0 * (ll1 - ll0_intercept);
    pval_intercept[j] = R::pchisq(lrt_intercept, 1.0, false, false);
    if (shrinkage && pval_intercept[j] > alpha) expect_intercept[j] = 0.0;
    
    // === Slope Test ===
    // Null model: theta = 0, intercept fixed
    double ll0_theta = univariate_loglik(x_j, y, family, 0.0, offset, expect_intercept[j]);
    double lrt_theta = 2.0 * (ll1 - ll0_theta);
    pval_theta[j] = R::pchisq(lrt_theta, 1.0, false, false);
    if (shrinkage && pval_theta[j] > alpha) expect_theta[j] = 0.0;

    // === Evidence ===
    evidence[j] = 2*(ll1 - ll0_theta) - std::log(n);
  }

  // Return results as a list
  return List::create(
    Named("loglik") = loglik,
    Named("bic") = bic,
    Named("bic_diff") = bic_diff,
    Named("bf") = bf,
    Named("pmp") = pmp,
    Named("intercept") = intercept,
    Named("theta") = theta,
    Named("pval_raw") = pval_raw,
    Named("pval_intercept") = pval_intercept,
    Named("pval_theta") = pval_theta,
    Named("evidence") = evidence,
    Named("evidence_raw") = evidence_raw,
    Named("expect_intercept") = expect_intercept,
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
    bool decompose = true,
    bool shrinkage = true,
    double alpha = 0.05,
    double tol = 5e-2,
    int max_iter = 100)
{
  // Start timing using standard C++ time
  clock_t start_time = clock();
  
  // Get dimensions
  int n = X.n_rows;
  int p = X.n_cols;

  List fam;
  bool is_cox = false;
  
  if (TYPEOF(family) == STRSXP) {
    std::string fam_str = as<std::string>(family);
    is_cox = (fam_str == "cox");
    if (is_cox) {
      fam = List::create(
        Named("family") = "cox",
        Named("link") = "log"
      );
    }
  } else if (TYPEOF(family) == VECSXP) {
    fam = as<List>(family);
    if (fam.containsElementNamed("family")) {
      std::string fam_str = as<std::string>(fam["family"]);
      is_cox = (fam_str == "cox");
    }
  } else {
    Rcpp::stop("family must be a family object; e.g., gaussian() or binomial()");
  }
  
  // Initialize parameters
  arma::mat intercept(p, L, arma::fill::zeros);
  arma::mat theta(p, L, arma::fill::zeros);

  // P-values
  arma::mat pval_raw(p, L, arma::fill::zeros);
  arma::mat pval_intercept(p, L, arma::fill::zeros);
  arma::mat pval_theta(p, L, arma::fill::zeros);

  // Evidence
  arma::mat evidence(p, L, arma::fill::zeros);
  arma::mat evidence_raw(p, L, arma::fill::zeros);
  
  // Initialize result matrices
  arma::mat loglik(p, L, arma::fill::zeros);
  arma::mat bic(p, L, arma::fill::zeros);
  arma::mat bic_diff(p, L, arma::fill::zeros);
  arma::mat bf(p, L, arma::fill::zeros);
  arma::mat pmp(p, L, arma::fill::zeros);
  arma::vec expect_variance(L, arma::fill::zeros);
  
  arma::mat ONES(n, p, arma::fill::ones);

  // Log-likelihood tracking
  arma::vec expect_loglik(max_iter, arma::fill::zeros);
  
  // Main iterations
  int iter;
  for (iter = 0; iter < max_iter; iter++) {
    if(iter > 1 && decompose) theta = decompose_theta(theta, L);// theta = arma::diagmat(arma::sum(theta, 1));

    // Calculate current linear predictor
    arma::vec linear_predictor(n, arma::fill::zeros);
    linear_predictor += ONES * intercept * arma::ones<arma::vec>(L);
    linear_predictor += X * theta * arma::ones<arma::vec>(L);
    
    // Update each single effect
    for (int l = 0; l < L; l++) {
      // Update offset by removing current effect
      arma::vec offset = linear_predictor - ONES * intercept.col(l) - X * theta.col(l);
      
      // Fit single effect
      List res = single_effect_fit(
        X,
        y,
        fam,
        offset,
        standardize,
        shrinkage,
        ties,
        lambda,
        tau,
        alpha
      );
      
      // Extract results
      arma::vec res_loglik = as<arma::vec>(res["loglik"]);
      arma::vec res_bic = as<arma::vec>(res["bic"]);
      arma::vec res_bic_diff = as<arma::vec>(res["bic_diff"]);
      arma::vec res_bf = as<arma::vec>(res["bf"]);
      arma::vec res_pmp = as<arma::vec>(res["pmp"]);
      
      arma::vec res_pval_intercept = as<arma::vec>(res["pval_intercept"]);
      arma::vec res_pval_theta = as<arma::vec>(res["pval_theta"]);
      arma::vec res_pval_raw = as<arma::vec>(res["pval_raw"]);

      arma::vec res_evidence = as<arma::vec>(res["evidence"]);
      arma::vec res_evidence_raw = as<arma::vec>(res["evidence_raw"]);

      // Apply thresholding to expect_theta
      arma::vec res_expect_intercept = as<arma::vec>(res["expect_intercept"]);
      arma::vec res_expect_theta = as<arma::vec>(res["expect_theta"]);
      
      // Update parameters
      intercept.col(l) = res_expect_intercept;
      theta.col(l) = res_expect_theta;
      pval_intercept.col(l) = res_pval_intercept;
      pval_theta.col(l) = res_pval_theta;
      pval_raw.col(l) = res_pval_raw;
      evidence.col(l) = res_evidence;
      evidence_raw.col(l) = res_evidence_raw;
      loglik.col(l) = res_loglik;
      bic.col(l) = res_bic;
      bic_diff.col(l) = res_bic_diff;
      bf.col(l) = res_bf;
      pmp.col(l) = res_pmp;
      expect_variance(l) = as<double>(res["expect_variance"]);
      
      // Update linear predictor for next effect
      linear_predictor = offset + ONES * intercept.col(l) + X * theta.col(l);
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
    Named("theta") = theta,
    Named("pval_intercept") = pval_intercept,
    Named("pval_theta") = pval_theta,
    Named("pval_raw") = pval_raw,
    Named("evidence") = evidence,
    Named("evidence_raw") = evidence_raw,
    Named("pmp") = pmp,
    Named("bic") = bic,
    Named("bic_diff") = bic_diff,
    Named("bf") = bf,
    Named("expect_variance") = expect_variance,
    Named("elapsed_time") = elapsed_time
  );
}