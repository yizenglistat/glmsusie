#pragma once

// src/core.h

#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// helper for sign
inline double sgn(double z) { return (z>0) - (z<0); }


/**
 * @brief Estimate GLM dispersion using Pearson or Deviance methods.
 *
 * @param y           Response vector (length n).
 * @param family  GLM family object (must contain "linkinv", "variance", "dev.resids").
 * @param offset      Optional linear predictor vector or scalar. Defaults to zero vector.
 * @param approach    Dispersion estimation method: "pearson" or "deviance".
 *
 * @return Estimated dispersion value (double).
 */
double update_dispersion(const arma::vec& y,
                         SEXP family,
                         arma::vec offset,
                         std::string approach);

/**
 * Calculate univariate Cox regression log-likelihood
 * 
 * @param x Vector of covariate values
 * @param y Matrix with 2 columns (time, status)
 * @param offset Vector of offset values
 * @param theta Coefficient value
 * @param ties Method for handling tied event times ("breslow" or "efron")
 * @return Log-likelihood value
 */
double univariate_loglik_cox(
    const arma::vec& x,
    const arma::mat& y,
    arma::vec offset,
    double theta,
    std::string ties
);

/**
 * Compute univariate log-likelihood for GLM or Cox model.
 *
 * @param x      Numeric covariate vector of length n
 * @param y      SEXP: numeric vector (GLM) or n×2 matrix for Cox
 * @param family Rcpp::List GLM family object or Cox family list
 * @param theta  Numeric coefficient to evaluate log-likelihood
 * @param offset Numeric offset (length 1 or n)
 * @param ties   Method for ties in Cox: "breslow" or "efron"
 * @return       Log-likelihood value
 */
double univariate_loglik(
    const arma::vec& x,
    SEXP             y,
    Rcpp::List       family,
    double           theta,
    const arma::vec& offset,
    std::string      ties
);


/**
 * Fit univariate Cox regression via iteratively reweighted least squares
 * 
 * @param x Vector of covariate values
 * @param y Matrix with 2 columns (time, status)
 * @param offset Vector of offset values
 * @param lambda Penalty weight; if ≤0 defaults to √(2 log n / n).
 * @param tau Truncation parameter; defaults to 1e-5
 * @param max_iter Maximum number of iterations
 * @param tol Convergence tolerance
 * @return Estimated coefficient value
 */
double univariate_irls_cox(
    arma::vec x, 
    arma::mat y,
    arma::vec offset,
    std::string ties,
    double lambda,
    double tau,
    int max_iter,
    double tol
);


/**
 * Calculate truncated‐L1 IRLS estimate for a univariate GLM.
 *
 * Solves the penalized objective
 *   min_theta [ - logLik(theta) + lambda * min(1, |theta|/tau) ]
 * via iteratively reweighted least squares with a closed‐form update for
 * the capped‐L1 penalty at each step.
 *
 * @param x          Numeric covariate vector (length n).
 * @param y          Numeric response vector (length n).
 * @param family     Rcpp::List GLM family object (gaussian, binomial, poisson).
 * @param offset     Offset vector (length 1 or n) for the linear predictor.
 * @param lambda     Penalty weight; if ≤0 defaults to √(2 log n / n).
 * @param tau        Truncation parameter; defaults to 1e-5
 * @param max_iter   Maximum number of IRLS iterations.
 * @param tol        Convergence tolerance on θ updates.
 * @return           Estimated coefficient θ minimizing the penalized objective.
 */
double univariate_irls_glm(
    const arma::vec& x,
    const arma::vec& y,
    SEXP             family,
    arma::vec        offset,
    double           lambda,
    double           tau,
    int              max_iter,
    double           tol
);

/**
 * Univariate fit for GLM or Cox with optional TLP penalty and standardization.
 *
 * @param x              Covariate vector (length n). Scalar expands to zeros.
 * @param y              Response: length-n vector (GLM) or n×2 matrix (Cox).
 * @param family         Rcpp::List: GLM family object or Cox family list.
 * @param offset         Offset vector (length 1 or n).
 * @param standardize    Whether to center and scale x (default: true).
 * @param ties           Cox ties method: "efron" or "breslow" (default: "efron").
 * @param lambda         Penalty weight; NULL defaults to sqrt(2*log(n)/n).
 * @param tau            Truncation param; NULL defaults to 0.5.
 * @param null_threshold Threshold below which theta is set to zero (default: 1e-6).
 *
 * @return List with elements:
 *   - theta:   Estimated coefficient.
 *   - loglik:  Unpenalized log-likelihood at theta.
 *   - bic:     Bayesian Information Criterion.
 */
Rcpp::List univariate_fit(
    const arma::vec&         x,
    SEXP                     y,
    SEXP                     family,
    arma::vec                offset,
    bool                     standardize,
    std::string              ties,
    double                   lambda,
    double                   tau,
    double                   null_threshold
);


/**
 * @brief Fit single-effect regression across all predictors
 * 
 * For each column of the design matrix X, fits a univariate model to
 * the response y, computes BIC differences relative to a null model,
 * converts these into Bayes factors and posterior model probabilities (PMP),
 * and returns both the MAP estimates and the PMP-weighted expectations.
 * 
 * @param X Design matrix (n × p)
 * @param y Response variable (vector for GLMs, matrix for Cox)
 * @param family Family object or "cox"
 * @param offset Offset vector
 * @param standardize Whether to standardize the covariates
 * @param ties Method for handling ties in Cox models
 * @param lambda Penalty strength parameter (NULL for default)
 * @param tau Truncation parameter (NULL for default)
 * @param null_threshold Threshold for setting coefficients to zero
 * 
 * @return List with log-likelihood, BIC, Bayes factors, PMP and coefficient estimates
 */
List single_effect_fit(
    const arma::mat& X,
    SEXP y,
    SEXP family,
    arma::vec offset,
    bool standardize,
    std::string ties,
    double lambda,
    double tau,
    double null_threshold
);

/**
 * @brief Fit Likelihood-based Additive Single-Effect Regression (LASER) Model
 * 
 * Fits the LASER model, representing the coefficient vector as a sum of L
 * sparse "single effects". At each iteration, it cyclically applies
 * single_effect_fit to update one effect while holding the others fixed,
 * then updates the intercept (and dispersion, if applicable), until
 * convergence or max_iter is reached.
 * 
 * @param X Design matrix (n × p)
 * @param y Response variable (vector for GLMs, matrix for Cox)
 * @param L Number of single effects to include (will be set to min(10, L))
 * @param family GLM family object or "cox" family
 * @param standardize Whether to standardize the covariates
 * @param ties Method for handling ties in Cox models
 * @param lambda Penalty strength parameter
 * @param tau Truncation parameter
 * @param null_threshold Threshold for setting values to zero
 * @param tol Convergence tolerance
 * @param max_iter Maximum number of iterations
 * 
 * @return List with fitted model components
 */
List additive_effect_fit(
    const arma::mat& X,
    SEXP y,
    int L,
    SEXP family,
    bool standardize,
    std::string ties,
    double lambda,
    double tau,
    double null_threshold,
    double tol,
    int max_iter
);








