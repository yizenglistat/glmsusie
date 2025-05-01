#include "core.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' @title Univariate Cox Partial Log-Likelihood
//' @description
//'   Computes the partial log‐likelihood for a univariate Cox proportional hazards model,
//'   allowing either the Breslow or Efron approximation to handle tied event times.
//'
//' @param x Numeric vector of covariate values (length n).
//' @param y Numeric matrix with two columns:  
//'   * column 1 = event or censoring times  
//'   * column 2 = status (1 = event, 0 = censored)
//'   Must have n rows.
//' @param offset Numeric vector or scalar; fixed component of the linear predictor.
//'   If scalar, it is recycled to length n.
//' @param theta Numeric scalar coefficient for x (default = 0).
//' @param ties Character; method for tied times.  
//'   * `"breslow"`: Breslow approximation  
//'   * `"efron"` (default): Efron approximation
//'
//' @return A single numeric value: the partial log‐likelihood \eqn{\ell(\theta)}.
//'
//' @examples
//' \dontrun{
//' n <- 100
//' x <- rnorm(n)
//' time <- rexp(n, rate = exp(0.3 * x))
//' status <- rbinom(n, 1, 0.7)
//' y <- cbind(time, status)
//'
//' # compute log‐likelihood at theta = 0.5
//' univariate_loglik_cox(x, y, offset = 0, theta = 0.5, ties = "efron")
//' }
//' @seealso \code{\link{univariate_irls_cox}}
//' @export
// [[Rcpp::export]]
double univariate_loglik_cox(
    const arma::vec& x,
    const arma::mat& y,
    arma::vec offset,
    double theta = 0.0,
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


//' @title Fit Univariate Cox via Newton–Raphson (IRLS)
//' @description
//'   Finds the MLE of the coefficient in a univariate Cox model by iteratively
//'   maximizing the partial likelihood using Newton–Raphson.
//'
//' @param x Numeric vector of covariate values (length n).
//' @param y Numeric matrix with two columns: event time and status.
//' @param offset Numeric vector or scalar offset term; recycled if scalar.
//' @param ties Character; how to handle ties, `"breslow"` or `"efron"` (default).
//' @param max_iter Integer; maximum number of NR iterations (default = 25).
//' @param tol Numeric; convergence tolerance on theta (default = 1e-8).
//'
//' @return A single numeric value: the estimated coefficient \eqn{\hat\theta}.
//'
//' @examples
//' \dontrun{
//' n <- 100
//' x <- rnorm(n)
//' y <- cbind(rexp(n, rate=exp(0.5*x)), rbinom(n,1,0.8))
//' theta_hat <- univariate_irls_cox(x, y, offset = 0, ties = "efron")
//' }
//' @seealso \code{\link{univariate_loglik_cox}}
//' @export
// [[Rcpp::export]]
double univariate_irls_cox(arma::vec x, 
                           arma::mat y,
                           arma::vec offset,
                           std::string ties = "efron",
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
  
  // Sort data by time
  arma::uvec ord = arma::stable_sort_index(time);
  arma::vec time_sorted = time.elem(ord);
  arma::vec status_sorted = status.elem(ord);
  arma::vec x_sorted = x.elem(ord);
  arma::vec offset_sorted = offset.elem(ord);
  
  // IRLS loop for Cox regression
  double theta = 0.0;
  bool converged = false;
  
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
        double sum_event_x = arma::sum(event_x);
        
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
    
    // Newton-Raphson update
    double theta_new = theta + score/info;
    
    // Convergence check
    if (std::abs(theta_new - theta) < tol) {
      theta = theta_new;
      converged = true;
      break;
    }
    
    theta = theta_new;
  }
  
  if (!converged) {
    Rcpp::warning("Cox IRLS algorithm did not converge in %d iterations", max_iter);
  }
  
  return theta;
}