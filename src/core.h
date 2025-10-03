#pragma once

// src/core.h

#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// //' and /** */ different commenting yields different RcppExport.

// helper for sign
inline double sgn(double z) { return (z>0) - (z<0); }


/**
 * @brief Decompose row sums of a theta matrix into a p x L matrix
 *
 * Given a p × L matrix `theta`, this function computes the row sums
 * and redistributes them into a new p × L matrix such that each row sum
 * is preserved, but concentrated in a single column (round-robin by row index).
 *
 * This is useful for reconstructing a single-effect decomposition where
 * each row's total effect is spread sparsely over L components.
 *
 * @param theta Input matrix of shape p × L
 * @param L     Number of columns (effects) in the output
 * @return      A p × L matrix with same row sums as `theta` and only one non-zero per row
 */
arma::mat decompose_theta(const arma::mat& theta, int L);

/**
 * Calculate univariate Cox regression log-likelihood
 * 
 * @param x Vector of covariate values
 * @param y Matrix with 2 columns (time, status)
 * @param theta Coefficient value
 * @param offset Vector of offset values
 * @param ties Method for handling tied event times ("breslow" or "efron")
 * @return Log-likelihood value
 */
double univariate_loglik_cox(
    const arma::vec& x,
    const arma::mat& y,
    double theta,
    arma::vec offset,
    std::string ties
);

/**
 * Calculate the log-likelihood for a univariate GLM model
 *
 * @param x Vector of predictor values
 * @param y Vector of response values
 * @param family An R family object (as SEXP)
 * @param theta Coefficient value
 * @param offset Vector of offset values (optional)
 * @param intercept
 * @return The calculated log-likelihood value
 */
double univariate_loglik_glm(
    const arma::vec& x,
    const arma::vec& y,
    SEXP family,
    double theta,
    const arma::vec& offset,
    double intercept
);

/**
 * Compute univariate log-likelihood for GLM or Cox model.
 *
 * @param x          Numeric covariate vector of length n
 * @param y          SEXP: numeric vector (GLM) or n×2 matrix for Cox
 * @param family     SEXP: GLM family object or Cox family list
 * @param theta      Numeric coefficient to evaluate log-likelihood
 * @param offset     Numeric offset (length 1 or n)
 * @param intercept  Numeric intercept to evaluate log-likelihood
 * @param ties       Method for ties in Cox: "breslow" or "efron"
 * @return           Log-likelihood value
 */
double univariate_loglik(
    const arma::vec& x,
    SEXP             y,
    SEXP             family,
    double           theta,
    const arma::vec& offset,
    double           intercept,
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
 * @brief Fit a univariate GLM (no intercept)
 * 
 * This function fits a generalized linear model (GLM) with a single predictor `x` and no intercept.
 * It is optimized for speed and numerical stability. The model form is:
 * 
 *     g(E[y]) = offset + theta * x
 * 
 * where g is the link function (e.g., logit, log, identity) determined by the specified GLM family.
 * 
 * @param x         Predictor vector (length n)
 * @param y         Response vector (length n)
 * @param family    R family object (e.g., binomial(), gaussian())
 * @param offset    Optional offset vector (length 1 or n). Use numeric(0) to omit.
 * @param max_iter  Maximum number of IRLS iterations
 * @param tol       Convergence tolerance
 * 
 * @return          Estimated slope (theta) as a double
 */
double univariate_irls_glm_no_intercept(const arma::vec& x,
                                        const arma::vec& y,
                                        SEXP             family,
                                        arma::vec        offset,
                                        int              max_iter,
                                        double           tol);

/**
 * @brief Fit a univariate generalized linear model
 * 
 * This function fits a generalized linear model (GLM) with a single predictor variable
 * plus an intercept term. It is optimized for the univariate case and implements the same
 * iteratively reweighted least squares (IRLS) algorithm as R's glm() function, but more 
 * efficiently for the single-predictor case.
 * 
 * The model form is:
 *    g(E[y]) = offset + intercept + theta * x
 * 
 * where g is the link function determined by the specified family.
 * 
 * @param x           Predictor variable (vector)
 * @param y           Response variable (vector)
 * @param family      R family object (e.g., gaussian(), binomial(), poisson())
 * @param offset      Optional offset term in the linear predictor (use numeric(0) for none)
 * @param max_iter    Maximum number of IRLS iterations
 * @param tol         Convergence tolerance for parameter estimates
 * 
 * @return A list containing:
 *    - intercept: The estimated intercept term
 *    - theta: The estimated coefficient for the predictor variable
 * 
 * @details
 * This function implements the IRLS algorithm for GLMs optimized for the univariate case.
 * It supports all standard GLM families available in R, including:
 *    - gaussian: Linear regression
 *    - binomial: Logistic regression
 *    - poisson: Poisson regression
 *    - Others supported by R's family objects
 * 
 * The implementation includes numerical safeguards for stability:
 *    - Handling of extreme values in the linear predictor
 *    - Prevention of zeros or negative values in variance calculations
 *    - Special handling for Poisson regression to prevent numerical issues
 *    - Multiple fallback methods for matrix solution if standard solve fails
 * 
 * @note
 * Results are identical to R's glm() function for the univariate case.
 * 
 */
Rcpp::List univariate_irls_glm(const arma::vec&   x,
                               const arma::vec&   y,
                               SEXP               family,
                               arma::vec          offset,
                               int                max_iter,
                               double             tol);


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
 *
 * @return List with elements:
 *   - theta:   Estimated intercept.
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
    double                   tau
);

/**
 * Univariate GLM fit using R's glm() interface, with optional offset.
 *
 * Fits a single-predictor generalized linear model (GLM) against the null model
 * and reports coefficients, standard errors, likelihood ratio test, and BIC-based
 * model comparison quantities. Internally calls R's stats::glm(), logLik(), BIC(),
 * and summary() functions.
 *
 * @param x        Covariate vector (length n).
 * @param y        Response vector (length n).
 * @param family   R family object (e.g., binomial(), gaussian(), poisson()).
 * @param offset   Optional offset vector (length 1 or n). If length 1, expanded.
 *
 * @return List with elements:
 *   - intercept: Intercept estimate from glm().
 *   - beta:      Slope coefficient estimate for x.
 *   - se:        Standard error of beta (from summary(glm)).
 *   - wald_p:    Wald test p-value for beta.
 *   - logLik1:   Log-likelihood of single-variable model.
 *   - logLik0:   Log-likelihood of null (intercept-only) model.
 *   - BIC1:      BIC of single-variable model.
 *   - BIC0:      BIC of null model.
 *   - LRT:       Likelihood ratio statistic (2*(ll1 - ll0)).
 *   - LRT_p:     Likelihood ratio test p-value (chi-square, df=1).
 *   - deltaBIC:  Difference BIC1 - BIC0 (negative favors variable model).
 *   - twoLogBF:  Approximate 2*log Bayes factor ≈ BIC0 - BIC1.
 */
Rcpp::List univariate_glm(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y,
    SEXP family,
    Rcpp::Nullable<Rcpp::NumericVector> offset
);

/**
 * Univariate Cox proportional hazards fit using R's survival::coxph(),
 * with optional offset and choice of ties handling.
 *
 * Fits a single-predictor Cox model against the null (no covariates) and
 * reports coefficient, standard error, Wald test, likelihood ratio test,
 * and BIC-based model comparison quantities. Internally calls
 * survival::coxph(), survival::summary.coxph(), and uses the partial
 * log-likelihood for information criteria.
 *
 * @param y        Response data, either an n×2 matrix (time, status) or a
 *                 list of two numeric vectors: (time, status).
 * @param x        Covariate vector (length n).
 * @param offset   Optional offset vector (length 1 or n). If length 1, expanded.
 * @param ties     Method for handling tied event times: "efron" (default) or "breslow".
 *
 * @return List with elements:
 *   - beta:      Estimated regression coefficient for x.
 *   - se:        Standard error of beta (from summary.coxph).
 *   - wald_z:    Wald test statistic.
 *   - wald_p:    Wald test p-value.
 *   - logLik0:   Null-model partial log-likelihood.
 *   - logLik1:   Single-variable model partial log-likelihood.
 *   - LRT:       Likelihood ratio statistic (2*(ll1 - ll0)).
 *   - LRT_p:     Likelihood ratio test p-value (chi-square, df=1).
 *   - BIC1:      BIC of single-variable model (partial likelihood).
 *   - BIC0:      BIC of null model (partial likelihood).
 *   - deltaBIC:  Difference BIC1 - BIC0 (negative favors variable model).
 *   - twoLogBF:  Approximate 2*log Bayes factor ≈ BIC0 - BIC1.
 *   - ties:      Ties handling method used.
 *   - has_offset: Logical indicating if an offset was provided.
 */
Rcpp::List univariate_cox(
    SEXP y,
    const Rcpp::NumericVector& x,
    Rcpp::Nullable<Rcpp::NumericVector> offset,
    std::string ties
);

/**
 * For each column of X, fits a univariate model (GLM or Cox, chosen by 'family'),
 * collects per-variable stats (coef, SE, p-values, BIC, evidence), forms Bayes
 * weights via a stabilized softmax on evidence, and (optionally) performs
 * shrinkage tests using univariate_loglik.
 *
 * @param X            n x p design matrix.
 * @param y            GLM vector or Cox (time,status).
 * @param family       GLM family or "cox".
 * @param offset       Optional offset (length 1 or n).
 * @param standardize  (kept for API compatibility; handled in inner fits as needed)
 * @param shrinkage    Whether to zero out expectations if p > alpha in LRT.
 * @param ties         Cox ties method ("efron" default).
 * @param lambda,tau   (kept for API compatibility; not used in glm/cox wrappers)
 * @param alpha        Significance level for shrinkage tests.
 *
 * @return List with elements:
 *   - loglik, bic, bic_diff, bf, pmp
 *   - intercept (length p; 0 for Cox), theta (length p), se_theta (length p)
 *   - pval_raw, pval_intercept, pval_theta
 *   - evidence, evidence_raw
 *   - expect_intercept, expect_theta, expect_variance
 */
Rcpp::List single_effect_fit(
    const arma::mat&   X,
    SEXP               y,
    SEXP               family,
    arma::vec          offset,
    bool               standardize,
    bool               shrinkage,
    std::string        ties,
    double             lambda,
    double             tau,
    double             alpha
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
 * @param decompose Whether to decompose the theta
 * @param shrinkage Whether to shrinkage parameters using pvals
 * @param alpha Level of significance
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
    bool decompose,
    bool shrinkage,
    double alpha,
    double tol,
    int max_iter
);








