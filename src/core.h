// src/core.h

#ifndef CORE_H
#define CORE_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
 * Fit univariate Cox regression via iteratively reweighted least squares
 * 
 * @param x Vector of covariate values
 * @param y Matrix with 2 columns (time, status)
 * @param offset Vector of offset values
 * @param max_iter Maximum number of iterations
 * @param tol Convergence tolerance
 * @return Estimated coefficient value
 */
double univariate_irls_cox(
    arma::vec x, 
    arma::mat y,
    arma::vec offset,
    std::string ties,
    int max_iter,
    double tol
);

#endif // CORE_H