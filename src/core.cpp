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
    

    // Solve for intercept and theta (silence Armadillo's "approx solution" warning)
    arma::vec params;
    bool solved = arma::solve(params, XtWX, XtWz, arma::solve_opts::no_approx);
    
    if (!solved) {
      params = arma::pinv(XtWX) * XtWz;
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

// // [[Rcpp::export]]
// List univariate_fit(
//     const arma::vec& x, 
//     SEXP y, 
//     SEXP family,
//     arma::vec offset,
//     bool standardize = true,
//     std::string ties = "efron",
//     double lambda = 0.0,
//     double tau = 0.5)
// {
//   // Get dimensions and extract family type
//   int n = x.n_elem;
//   bool is_cox = TYPEOF(family) == STRSXP ? as<std::string>(family) == "cox" : 
//                  as<std::string>(as<List>(family)["family"]) == "cox";
  
//   // Handle offset
//   if (offset.n_elem == 1 && n > 1) offset = arma::vec(n, arma::fill::value(offset(0)));
  
//   // Standardize predictor (only if non-constant)
//   double x_mean = 0.0, x_norm = 1.0;
//   arma::vec x_std = x;
  
//   if (standardize) {
//     x_mean = arma::mean(x);
//     x_std = x - x_mean;
//     x_norm = arma::norm(x_std, 2);
//     if (x_norm > 1e-10) x_std /= x_norm;
//     else { x_std = x; x_mean = 0.0; x_norm = 1.0; }

//   }
  
//   // Fit model using appropriate function
//   double intercept_std = 0.0;
//   double theta_std = 0.0;
  
//   if (is_cox) {
//     // Convert y to appropriate format for Cox
//     arma::mat y_mat = TYPEOF(y) == VECSXP ? 
//       arma::mat(n, 2) : as<arma::mat>(y);
    
//     if (TYPEOF(y) == VECSXP) {
//       NumericVector time = as<NumericVector>(as<List>(y)[0]);
//       NumericVector status = as<NumericVector>(as<List>(y)[1]);
//       y_mat.col(0) = as<arma::vec>(time);
//       y_mat.col(1) = as<arma::vec>(status);
//     }
    
//     theta_std = univariate_irls_cox(x_std, y_mat, offset, ties, lambda, tau);
//   } else {
//     // Convert y for GLM
//     arma::vec y_vec = TYPEOF(y) == VECSXP ? 
//       as<arma::vec>(as<NumericVector>(as<List>(y)[0])) : 
//       as<arma::vec>(as<NumericVector>(y));
    
//     Rcpp::List res = univariate_irls_glm(x_std, y_vec, family, offset);

//     intercept_std = as<double>(res["intercept"]);
//     theta_std     = as<double>(res["theta"]);

//   }
  
//   // Transform coefficient, apply threshold, and compute metrics
//   double theta = theta_std / x_norm;
//   double intercept = intercept_std - theta * x_mean;

//   double loglik = univariate_loglik(x, y, family, theta, offset, intercept, ties);
//   double loglik0 = univariate_loglik(x, y, family, 0.0, offset, intercept, ties);
//   double lrt = 2.0 * (loglik - loglik0);
//   double pval = (lrt > 0) ? R::pchisq(lrt, 1.0, false, true) : 1.0;
//   double bic = -2.0 * loglik + std::log(n) * (theta != 0.0 ? 1.0 : 0.0);

//   return List::create(
//     Named("intercept") = intercept,
//     Named("theta") = theta,
//     Named("loglik") = loglik,
//     Named("bic") = bic,
//     Named("pval") = pval
//   );
// }

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
  const int n = x.n_elem;

  // Identify family type
  const bool is_cox = TYPEOF(family) == STRSXP
    ? as<std::string>(family) == "cox"
    : as<std::string>(as<List>(family)["family"]) == "cox";

  // Broadcast offset if needed
  if (offset.n_elem == 1 && n > 1)
    offset = arma::vec(n, arma::fill::value(offset(0)));

  // Prepare response in the desired shape
  arma::vec y_vec;      // for GLM
  arma::mat y_mat;      // for Cox (time, status)
  if (is_cox) {
    if (TYPEOF(y) == VECSXP) {
      NumericVector time   = as<NumericVector>(as<List>(y)[0]);
      NumericVector status = as<NumericVector>(as<List>(y)[1]);
      y_mat.set_size(n, 2);
      y_mat.col(0) = as<arma::vec>(time);
      y_mat.col(1) = as<arma::vec>(status);
    } else {
      y_mat = as<arma::mat>(y);
    }
  } else {
    y_vec = (TYPEOF(y) == VECSXP)
      ? as<arma::vec>(as<NumericVector>(as<List>(y)[0]))
      : as<arma::vec>(as<NumericVector>(y));
  }

  // Standardize predictor (if non-constant)
  double x_mean = 0.0, x_norm = 1.0;
  arma::vec x_std = x;
  if (standardize) {
    x_mean = arma::mean(x);
    x_std  = x - x_mean;
    x_norm = arma::norm(x_std, 2);
    if (x_norm > 1e-10) {
      x_std /= x_norm;
    } else {
      x_std = x; x_mean = 0.0; x_norm = 1.0;
    }
  }

  // ---------- Fit alternative model (θ free) ----------
  double intercept_std = 0.0;
  double theta_std     = 0.0;

  if (is_cox) {
    // Cox: no intercept in partial likelihood
    theta_std = univariate_irls_cox(x_std, y_mat, offset, ties, lambda, tau);
  } else {
    Rcpp::List res = univariate_irls_glm(x_std, y_vec, family, offset);
    intercept_std = as<double>(res["intercept"]);
    theta_std     = as<double>(res["theta"]);
  }

  // Back-transform to original x-scale
  const double theta     = theta_std / x_norm;
  const double intercept = intercept_std - theta * x_mean;

  // Alt log-likelihood at its own MLEs
  const double loglik = univariate_loglik(x, y, family, theta, offset, intercept, ties);

  // ---------- Fit null model properly (θ = 0 with its own MLE intercept) ----------
  double intercept0 = 0.0;
  if (!is_cox) {
    // GLM: intercept-only fit with same offset
    Rcpp::List res0 = univariate_irls_glm(arma::zeros<arma::vec>(n), y_vec, family, offset);
    intercept0 = as<double>(res0["intercept"]);
  } // Cox: no intercept; β=0 is the null

  const double loglik0 = univariate_loglik(x, y, family, 0.0, offset, intercept0, ties);

  // Likelihood ratio test, p-value, and BIC/evidence
  const double lrt   = 2.0 * (loglik - loglik0);
  const double lrt_p = (lrt > 0.0) ? R::pchisq(lrt, 1.0, /*lower_tail=*/false, /*log_p=*/false) : 1.0;

  // BIC for the alternative model (k = 1 additional slope parameter; intercept cancels in ΔBIC)
  const double bic = -2.0 * loglik + std::log((double)n) * (theta != 0.0 ? 1.0 : 0.0);

  // Evidence on the "2*log BF" scale: evidence = 2*log(BF_10) = -ΔBIC = 2(ℓ1-ℓ0) - log n
  const double evidence = lrt - std::log((double)n);

  // ---------- Profiled SE via finite-diff Hessian in θ ----------
  // Use the *profile* log-likelihood: at θ±h, re-maximize the intercept (GLM).
  const double h = std::max(1e-6, 1e-4 * (1.0 + std::abs(theta)));

  auto intercept_only_glm = [&](double theta_fix) {
    // absorb θ * x into the offset and refit an intercept-only GLM
    arma::vec off = offset + theta_fix * x; // NOTE: theta on original x-scale
    Rcpp::List r = univariate_irls_glm(arma::zeros<arma::vec>(n), y_vec, family, off);
    return as<double>(r["intercept"]);
  };

  double ll_plus  = NA_REAL, ll_minus = NA_REAL;

  if (is_cox) {
    // Cox: no intercept to profile
    ll_plus  = univariate_loglik(x, y, family, theta + h, offset, 0.0, ties);
    ll_minus = univariate_loglik(x, y, family, theta - h, offset, 0.0, ties);
  } else {
    const double intercept_plus  = intercept_only_glm(theta + h);
    const double intercept_minus = intercept_only_glm(theta - h);
    ll_plus  = univariate_loglik(x, y, family, theta + h, offset, intercept_plus,  ties);
    ll_minus = univariate_loglik(x, y, family, theta - h, offset, intercept_minus, ties);
  }

  const double d2   = (ll_plus - 2.0*loglik + ll_minus) / (h*h); // second derivative
  const double info = -d2;                                       // observed information

  double std_err = NA_REAL, z = NA_REAL, p_wald = NA_REAL;
  if (std::isfinite(info) && info > 0.0) {
    std_err = std::sqrt(1.0 / info);
    if (R_finite(std_err) && std_err > 0.0) {
      z = theta / std_err;
      p_wald = 2.0 * R::pnorm(std::fabs(z), 0.0, 1.0, /*lower_tail=*/false, /*log_p=*/false);
    }
  }

  return List::create(
    Named("intercept")    = intercept,
    Named("theta")        = theta,
    Named("std_err")      = std_err,        // profiled SE
    Named("p_wald")       = p_wald,         // Wald p-value
    Named("lrt")          = lrt,
    Named("lrt_p_value")  = lrt_p,
    Named("loglik")       = loglik,
    Named("loglik0")      = loglik0,        // (added for debugging/verification)
    Named("bic")          = bic,
    Named("evidence")     = evidence        // = 2*log BF_10 = -ΔBIC
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
  arma::vec pval_intercept(p, arma::fill::zeros);
  arma::vec pval_theta(p, arma::fill::zeros);
  arma::vec evidence(p, arma::fill::zeros);
  arma::vec pval_wald(p, arma::fill::zeros);
  arma::vec std_err(p, arma::fill::zeros);

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
    evidence[j] = as<double>(res["evidence"]);
    pval_wald[j] = as<double>(res["p_wald"]);
    std_err[j] = as<double>(res["std_err"]);

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
    double ll1 = univariate_loglik(x_j, y, family, theta[j], offset, intercept[j]);

    // === Intercept Test ===
    // Null model: intercept = 0, theta fixed
    double ll0_intercept = univariate_loglik(x_j, y, family, theta[j], offset, 0.0);
    double lrt_intercept = 2.0 * (ll1 - ll0_intercept);
    pval_intercept[j] = R::pchisq(lrt_intercept, 1.0, false, false);
    if (shrinkage && pval_intercept[j] > alpha) expect_intercept[j] = 0.0;
    
    // === Slope Test ===
    // Null model: theta = 0, intercept fixed
    double ll0_theta = univariate_loglik(x_j, y, family, 0.0, offset, intercept[j]);
    double lrt_theta = 2.0 * (ll1 - ll0_theta);
    pval_theta[j] = R::pchisq(lrt_theta, 1.0, false, false);
    if (shrinkage && pval_theta[j] > alpha) expect_theta[j] = 0.0;
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
    Named("pval_intercept") = pval_intercept,
    Named("pval_theta") = pval_theta,
    Named("evidence") = evidence,
    Named("pval_wald") = pval_wald,
    Named("std_err") = std_err,
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
  arma::mat pval_intercept(p, L, arma::fill::zeros);
  arma::mat pval_theta(p, L, arma::fill::zeros);

  // Evidence
  arma::mat evidence(p, L, arma::fill::zeros);
  arma::mat pval_wald(p, L, arma::fill::zeros);
  arma::mat std_err(p, L, arma::fill::zeros);
  
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

      arma::vec res_evidence = as<arma::vec>(res["evidence"]);
      arma::vec res_pval_wald = as<arma::vec>(res["pval_wald"]);
      arma::vec res_std_err = as<arma::vec>(res["std_err"]);

      // Apply thresholding to expect_theta
      arma::vec res_expect_intercept = as<arma::vec>(res["expect_intercept"]);
      arma::vec res_expect_theta = as<arma::vec>(res["expect_theta"]);
      
      // Update parameters
      intercept.col(l) = res_expect_intercept;
      theta.col(l) = res_expect_theta;
      pval_intercept.col(l) = res_pval_intercept;
      pval_theta.col(l) = res_pval_theta;
      evidence.col(l) = res_evidence;
      pval_wald.col(l) = res_pval_wald;
      std_err.col(l) = res_std_err;
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
    Named("evidence") = evidence,
    Named("pval_wald") = pval_wald,
    Named("std_err") = std_err,
    Named("pmp") = pmp,
    Named("bic") = bic,
    Named("bic_diff") = bic_diff,
    Named("bf") = bf,
    Named("expect_variance") = expect_variance,
    Named("elapsed_time") = elapsed_time
  );
}