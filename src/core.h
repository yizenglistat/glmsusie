#pragma once

// core.h — Core APIs for LASER model fitting, single‐effect regressions, and utilities

#include <RcppArmadillo.h>
using namespace Rcpp;

/**
 * @brief Fit a univariate GLM by IRLS or closed‐form (Gaussian).
 *
 * @param x              Predictor vector (length n).
 * @param y              Response vector (length n).
 * @param family         GLM family object (must contain “family” or “name”).
 * @param offset         Offset vector (length n) or scalar.
 * @param standardize    Whether to center & scale x before fitting.
 * @param max_iter       Maximum IRLS/Newton iterations.
 * @param tol            Convergence tolerance on coefficient change.
 * @return Named list with components:
 *   - theta: estimated coefficient,
 *   - se: standard error,
 *   - p_value: two‐sided p‐value,
 *   - logLik: log‐likelihood,
 *   - BIC: Bayesian Information Criterion.
 */
List get_glm_fit(arma::vec x,
                 arma::vec y,
                 List    family,
                 arma::vec offset,
                 bool    standardize = true,
                 int     max_iter    = 25,
                 double  tol         = 1e-8);


/**
 * @brief Fit a univariate Cox proportional hazards model.
 *
 * @param x              Predictor vector (length n).
 * @param y              Survival matrix (n x 2: time, status).
 * @param offset         Offset vector (length n) or scalar.
 * @param standardize    Whether to center & scale x before fitting.
 * @param ties    Ties handling ("efron" or "breslow").
 * @param max_iter       Maximum Newton‐Raphson iterations.
 * @param tol            Convergence tolerance on coefficient change.
 * @return Named list with components:
 *   - theta: estimated log‐hazard ratio,
 *   - se: standard error,
 *   - p_value: Wald‐test p‐value,
 *   - logLik: partial log‐likelihood,
 *   - BIC: Bayesian Information Criterion.
 */
List get_cox_fit(arma::vec    x,
                 arma::mat    y,
                 arma::vec    offset,
                 bool         standardize = true,
                 std::string  ties = "efron",
                 int          max_iter    = 25,
                 double       tol         = 1e-8);


/**
 * @brief Fit a null GLM (intercept + offset only).
 *
 * @param y        Response vector (length n).
 * @param family   GLM family object.
 * @param offset   Offset vector (length n) or scalar.
 * @param max_iter Maximum IRLS iterations for one‐parameter families.
 * @param tol      Convergence tolerance.
 * @return Named list with:
 *   - logLik: log‐likelihood,
 *   - BIC: Bayesian Information Criterion.
 */
List null_glm_fit(arma::vec y,
                  List    family,
                  arma::vec offset,
                  int     max_iter = 25,
                  double  tol      = 1e-8);


/**
 * @brief Fit a null Cox model (offset only).
 *
 * @param y            Survival matrix (n x 2: time, status).
 * @param offset       Offset vector (length n) or scalar.
 * @param ties  Ties handling ("efron" or "breslow").
 * @return Named list with:
 *   - logLik: partial log‐likelihood,
 *   - BIC: Bayesian Information Criterion.
 */
List null_cox_fit(arma::mat   y,
                  arma::vec   offset,
                  std::string ties = "efron");


/**
 * @brief Given a known offset, update the intercept parameter.
 *
 * For GLMs, solves the single‐parameter intercept via Newton‐Raphson.
 * For Cox, returns 0.
 *
 * @param y       Response (vector or 2-col matrix for Cox).
 * @param family  Family object with “family” or “name”.
 * @param offset  Offset vector (length n).
 * @return Estimated intercept.
 */
double update_intercept(SEXP  y,
                        List  family,
                        arma::vec offset);


/**
 * @brief Given a known offset, update the dispersion parameter.
 *
 * - Gaussian/Gamma/Inverse‐Gaussian: method‐of‐moments or deviance‐based.
 * - Binomial/Poisson/Cox: returns 1 (fixed).
 *
 * @param y       Response (vector or 2-col matrix for Cox).
 * @param family  Family object with “family” or “name”.
 * @param offset  Offset vector (length n).
 * @return Updated dispersion.
 */
double update_dispersion(SEXP  y,
                         List  family,
                         arma::vec offset);


/**
 * @brief Compute total log-likelihood under current LASER parameters.
 *
 * Supports GLM families and Cox.  Summation over L single‐effect contributions.
 *
 * @param X           Design matrix (n × p).
 * @param y           Response (vector or 2-col matrix for Cox).
 * @param family      GLM family object or string "cox".
 * @param theta       p×L matrix of per‐effect coefficients.
 * @param intercept   Intercept scalar.
 * @param dispersion  Dispersion (for Gaussian/Gamma/IG families).
 * @return Total log-likelihood.
 */
double get_loglike(NumericMatrix X,
                   SEXP          y,
                   SEXP          family,
                   arma::mat     theta,
                   double        intercept,
                   double        dispersion = 1.0);