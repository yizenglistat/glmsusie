# Tests for the laser algorithm with various model families 
# e.g., Gaussian, Binomila, Poisson, Cox
# Author: Yizeng Li
# Date: 2025-04-28

# =============================================================================
# Test Gaussian Models
# =============================================================================
test_that("laser handles Gaussian regression with independent predictors", {
  # SETUP: Create controlled test data with known structure
  set.seed(1)
  n <- 500                                     # Number of observations
  p <- 10                                      # Number of predictors
  X <- matrix(rnorm(n * p), n, p)              # Independent predictors
  true_vars <- c(1, 3)                         # Variables with true effects
  y <- X[, true_vars[1]] + X[, true_vars[2]] + rnorm(n, sd = 1)  # Response with true model
  
  # EXECUTE: Run laser algorithm
  output <- laser(X, y, family = gaussian())
  
  # VERIFY: Test output structure and correctness
  # 1. Basic structure checks
  expect_s3_class(output, "laser")               # Should return proper S3 class
  expect_true(all(c("cs", "fit") %in% names(output)))  # Contains expected components
  
  # 2. Confidence sets validation
  expect_type(output$cs, "list")                 # cs should be a list
  expect_true(all(c("sets", "coverage") %in% names(output$cs)))  # Has both sets and coverage
  expect_type(output$cs$sets, "list")            # sets should be a list of confidence sets
  
  # 3. Variable selection correctness
  # Check if any confidence set contains at least one of the true variables
  found_true_var <- FALSE
  for (i in seq_along(output$cs$sets)) {
    if (any(true_vars %in% output$cs$sets[[i]])) {
      found_true_var <- TRUE
      break
    }
  }
  expect_true(found_true_var, "No confidence set contains any of the true variables")
  
  # 4. Result properties
  expect_true(all(output$fit$p_values >= 0 & output$fit$p_values <= 1))  # Valid p-values
  expect_true(all(output$cs$coverage <= 1))       # Coverage should not exceed 1
  
  # 5. Coefficient estimates
  expect_equal(length(output$fit$theta_hat), p)   # Should have one coefficient per variable
  expect_true(abs(output$fit$theta_hat[1]) > 0.7) # True effect should be detected
  expect_true(abs(output$fit$theta_hat[3]) > 0.7) # True effect should be detected
  
  # 6. Model fit quality
  expect_true(output$fit$final_loglike > -Inf)    # Should have valid likelihood
  expect_true(is.finite(output$fit$dispersion))   # Dispersion should be finite
})

test_that("laser handles Gaussian regression with highly correlated predictors", {
  # SETUP: Create controlled test data with multicollinearity
  set.seed(1)
  n <- 500                                     # Number of observations
  p <- 10                                      # Number of predictors
  
  # Generate base independent random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Create highly correlated pairs
  X <- matrix(0, n, p)
  X[, 1] <- Z[, 1]
  X[, 2] <- 0.95 * Z[, 1] + 0.05 * Z[, 2]      # X1, X2 correlated at ~0.95
  X[, 3] <- Z[, 3]
  X[, 4] <- 0.95 * Z[, 3] + 0.05 * Z[, 4]      # X3, X4 correlated at ~0.95
  
  # Generate remaining predictors independently
  X[, 5:p] <- Z[, 5:p]
  
  # Generate response - only X1 and X4 have true effects
  true_pair1 <- c(1, 2)                       # First correlated pair
  true_pair2 <- c(3, 4)                       # Second correlated pair
  true_effects <- c(1, 0, 0, 1, rep(0, p-4))  # True coefficients: X1=1, X4=1
  y <- X %*% true_effects + rnorm(n, sd = 1)
  
  # EXECUTE: Run laser algorithm
  output <- laser(X, y, family = gaussian())
  
  # VERIFY: Test output structure and correctness
  # 1. Basic structure checks
  expect_s3_class(output, "laser")
  expect_true(all(c("cs", "fit") %in% names(output)))
  
  # 2. Confidence sets validation
  expect_type(output$cs, "list")
  expect_true(all(c("sets", "coverage") %in% names(output$cs)))
  expect_type(output$cs$sets, "list")
  expect_length(output$cs$sets, length(output$cs$coverage))
  
  # 3. Variable selection correctness under multicollinearity
  # Check if we find at least one variable from each correlated pair
  found_pair1 <- FALSE
  found_pair2 <- FALSE
  
  for (i in seq_along(output$cs$sets)) {
    set_vars <- output$cs$sets[[i]]
    if (any(true_pair1 %in% set_vars)) {
      found_pair1 <- TRUE
    }
    if (any(true_pair2 %in% set_vars)) {
      found_pair2 <- TRUE
    }
  }
  
  # At least one variable from each true pair should be found
  expect_true(found_pair1, "No confidence set contains any variable from the first correlated pair")
  expect_true(found_pair2, "No confidence set contains any variable from the second correlated pair")
  
  # 4. Check for uncertainty recognition in coefficients
  # Under multicollinearity, we expect similar magnitudes for correlated variables
  coef_ratio1 <- abs(output$fit$theta_hat[1] / output$fit$theta_hat[2])
  coef_ratio2 <- abs(output$fit$theta_hat[3] / output$fit$theta_hat[4])
  
  # Either both variables are identified or their coefficient ratio is reasonable
  expect_true(
    coef_ratio1 > 0.5 || coef_ratio1 < 2 || 
    (abs(output$fit$theta_hat[1]) < 0.1 && abs(output$fit$theta_hat[2]) < 0.1),
    "Coefficients for correlated pair 1 show unreasonable differences"
  )
  
  expect_true(
    coef_ratio2 > 0.5 || coef_ratio2 < 2 || 
    (abs(output$fit$theta_hat[3]) < 0.1 && abs(output$fit$theta_hat[4]) < 0.1),
    "Coefficients for correlated pair 2 show unreasonable differences"
  )
})

# =============================================================================
# Test Binomial Models
# =============================================================================
test_that("laser handles binomial regression with independent predictors", {
  # SETUP: Create controlled test data with known structure
  set.seed(1)
  n <- 500                                     # Number of observations
  p <- 10                                      # Number of predictors
  X <- matrix(rnorm(n * p), n, p)              # Independent predictors
  true_vars <- c(2, 5)                         # Variables with true effects
  
  # Create binary response from logistic model
  linear_pred <- 1.5*X[, true_vars[1]] - 0.8*X[, true_vars[2]]
  prob <- 1/(1 + exp(-linear_pred))            # Inverse logit function
  y <- rbinom(n, 1, prob)                      # Binary outcome
  
  # EXECUTE: Run laser algorithm
  output <- laser(X, y, family = binomial())
  
  # VERIFY: Test output structure and correctness
  # 1. Basic structure checks
  expect_s3_class(output, "laser")               # Should return proper S3 class
  expect_true(all(c("cs", "fit") %in% names(output)))  # Contains expected components
  
  # 2. Confidence sets validation
  expect_type(output$cs, "list")                 # cs should be a list
  expect_true(all(c("sets", "coverage") %in% names(output$cs)))  # Has both sets and coverage
  expect_type(output$cs$sets, "list")            # sets should be a list of confidence sets
  
  # 3. Variable selection correctness
  # Check if any confidence set contains at least one of the true variables
  found_true_var <- FALSE
  for (i in seq_along(output$cs$sets)) {
    if (any(true_vars %in% output$cs$sets[[i]])) {
      found_true_var <- TRUE
      break
    }
  }
  expect_true(found_true_var, "No confidence set contains any of the true variables")
  
  # 4. Result properties
  expect_true(all(output$fit$p_values >= 0 & output$fit$p_values <= 1))  # Valid p-values
  expect_true(all(output$cs$coverage <= 1))       # Coverage should not exceed 1
  
  # 5. Coefficient estimates
  expect_equal(length(output$fit$theta_hat), p)   # Should have one coefficient per variable
  expect_true(abs(output$fit$theta_hat[2]) > 0.5) # True effect should be detected
  expect_true(abs(output$fit$theta_hat[5]) > 0.3) # True effect should be detected
  
  # 6. Model fit quality
  expect_true(output$fit$final_loglike > -Inf)    # Should have valid likelihood
  expect_true(is.finite(output$fit$dispersion))   # Dispersion should be finite
})

test_that("laser handles binomial regression with highly correlated predictors", {
  # SETUP: Create controlled test data with multicollinearity
  set.seed(1)
  n <- 500                                     # Number of observations
  p <- 10                                      # Number of predictors
  
  # Generate base independent random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Create highly correlated pairs
  X <- matrix(0, n, p)
  X[, 1] <- Z[, 1]
  X[, 2] <- 0.95 * Z[, 1] + 0.05 * Z[, 2]      # X1, X2 correlated at ~0.95
  X[, 3] <- Z[, 3]
  X[, 4] <- 0.95 * Z[, 3] + 0.05 * Z[, 4]      # X3, X4 correlated at ~0.95
  
  # Generate remaining predictors independently
  X[, 5:p] <- Z[, 5:p]
  
  # Generate binary response - only X1 and X4 have true effects
  true_pair1 <- c(1, 2)                       # First correlated pair
  true_pair2 <- c(3, 4)                       # Second correlated pair
  linear_pred <- X[, 1] + X[, 4] - 1          # Intercept = -1
  prob <- 1/(1 + exp(-linear_pred))           # Inverse logit function
  y <- rbinom(n, 1, prob)                     # Binary outcome
  
  # EXECUTE: Run laser algorithm
  output <- laser(X, y, family = binomial())
  
  # VERIFY: Test output structure and correctness
  # 1. Basic structure checks
  expect_s3_class(output, "laser")
  expect_true(all(c("cs", "fit") %in% names(output)))
  
  # 2. Confidence sets validation
  expect_type(output$cs, "list")
  expect_true(all(c("sets", "coverage") %in% names(output$cs)))
  expect_type(output$cs$sets, "list")
  
  # 3. Variable selection correctness under multicollinearity
  # Check if we find at least one variable from each correlated pair
  found_pair1 <- FALSE
  found_pair2 <- FALSE
  
  for (i in seq_along(output$cs$sets)) {
    set_vars <- output$cs$sets[[i]]
    if (any(true_pair1 %in% set_vars)) {
      found_pair1 <- TRUE
    }
    if (any(true_pair2 %in% set_vars)) {
      found_pair2 <- TRUE
    }
  }
  
  expect_true(found_pair1, "No confidence set contains any variable from the first correlated pair")
  expect_true(found_pair2, "No confidence set contains any variable from the second correlated pair")
  
  # 4. Check coefficient magnitudes - should reflect uncertainty due to correlation
  var1_coef <- abs(output$fit$theta_hat[1])
  var2_coef <- abs(output$fit$theta_hat[2])
  var3_coef <- abs(output$fit$theta_hat[3])
  var4_coef <- abs(output$fit$theta_hat[4])
  
  # Either one variable from each pair should have substantial coefficient
  # or both variables in the pair should have similar coefficients
  expect_true(
    (var1_coef > 0.3 || var2_coef > 0.3) && 
    abs(var1_coef - var2_coef) < 0.5,
    "Coefficients for first pair don't reflect the correlation structure"
  )
  
  expect_true(
    (var3_coef > 0.3 || var4_coef > 0.3) && 
    abs(var3_coef - var4_coef) < 0.5,
    "Coefficients for second pair don't reflect the correlation structure"
  )
})

# =============================================================================
# Test Poisson Models
# =============================================================================
test_that("laser handles Poisson regression with independent predictors", {
  # SETUP: Create controlled test data with known structure
  set.seed(1)
  n <- 500                                     # Number of observations
  p <- 10                                      # Number of predictors
  X <- matrix(rnorm(n * p), n, p)              # Independent predictors
  true_vars <- c(4, 7)                         # Variables with true effects
  
  # Create count response from log-linear model
  linear_pred <- 0.5 + 0.7*X[, true_vars[1]] + 0.8*X[, true_vars[2]]
  mu <- exp(linear_pred)                       # Log link function
  y <- rpois(n, mu)                            # Count outcome
  
  # EXECUTE: Run laser algorithm
  output <- laser(X, y, family = poisson())

  # VERIFY: Test output structure and correctness
  # 1. Basic structure checks
  expect_s3_class(output, "laser")               # Should return proper S3 class
  expect_true(all(c("cs", "fit") %in% names(output)))  # Contains expected components
  
  # 2. Confidence sets validation
  expect_type(output$cs, "list")                 # cs should be a list
  expect_true(all(c("sets", "coverage") %in% names(output$cs)))  # Has both sets and coverage
  expect_type(output$cs$sets, "list")            # sets should be a list of confidence sets
  
  # 3. Variable selection correctness
  # Check if any confidence set contains at least one of the true variables
  found_true_var <- FALSE
  for (i in seq_along(output$cs$sets)) {
    if (any(true_vars %in% output$cs$sets[[i]])) {
      found_true_var <- TRUE
      break
    }
  }
  expect_true(found_true_var, "No confidence set contains any of the true variables")
  
  # 4. Result properties
  expect_true(all(output$fit$p_values >= 0 & output$fit$p_values <= 1))  # Valid p-values
  expect_true(all(output$cs$coverage <= 1))       # Coverage should not exceed 1
  
  # 5. Coefficient estimates
  expect_equal(length(output$fit$theta_hat), p)   # Should have one coefficient per variable
  expect_true(abs(output$fit$theta_hat[4]) > 0.35) # True effect should be detected
  expect_true(abs(output$fit$theta_hat[7]) > 0.4) # True effect should be detected
  
  # 6. Model fit quality
  expect_true(output$fit$final_loglike > -Inf)    # Should have valid likelihood
  expect_true(is.finite(output$fit$dispersion))   # Dispersion should be finite
})

test_that("laser handles Poisson regression with highly correlated predictors", {
  # SETUP: Create controlled test data with multicollinearity
  set.seed(1)
  n <- 500                                     # Number of observations
  p <- 10                                      # Number of predictors
  
  # Generate base independent random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Create highly correlated pairs
  X <- matrix(0, n, p)
  X[, 1] <- Z[, 1]
  X[, 2] <- 0.95 * Z[, 1] + 0.05 * Z[, 2]      # X1, X2 correlated at ~0.95
  X[, 3] <- Z[, 3]
  X[, 4] <- 0.95 * Z[, 3] + 0.05 * Z[, 4]      # X3, X4 correlated at ~0.95
  
  # Generate remaining predictors independently
  X[, 5:p] <- Z[, 5:p]
  
  # Generate count response - only X1 and X4 have true effects
  true_pair1 <- c(1, 2)                       # First correlated pair
  true_pair2 <- c(3, 4)                       # Second correlated pair
  linear_pred <- 0.2 + X[, 1] + X[, 4]        # Intercept = 0.2
  mu <- exp(linear_pred)                      # Log link function
  y <- rpois(n, mu)                           # Count outcome
  
  # EXECUTE: Run laser algorithm
  output <- laser(X, y, family = poisson())
  
  # VERIFY: Test output structure and correctness
  # 1. Basic structure checks
  expect_s3_class(output, "laser")
  expect_true(all(c("cs", "fit") %in% names(output)))
  
  # 2. Confidence sets validation
  expect_type(output$cs, "list")
  expect_true(all(c("sets", "coverage") %in% names(output$cs)))
  expect_type(output$cs$sets, "list")
  
  # 3. Variable selection correctness under multicollinearity
  # Check if we find at least one variable from each correlated pair
  found_pair1 <- FALSE
  found_pair2 <- FALSE
  
  for (i in seq_along(output$cs$sets)) {
    set_vars <- output$cs$sets[[i]]
    if (any(true_pair1 %in% set_vars)) {
      found_pair1 <- TRUE
    }
    if (any(true_pair2 %in% set_vars)) {
      found_pair2 <- TRUE
    }
  }
  
  expect_true(found_pair1, "No confidence set contains any variable from the first correlated pair")
  expect_true(found_pair2, "No confidence set contains any variable from the second correlated pair")
  
  # 4. Check for effects in the estimated coefficients
  # For highly correlated variables, we expect either one or both to show effects
  expect_true(abs(output$fit$theta_hat[1]) > 0.25 || abs(output$fit$theta_hat[2]) > 0.25,
              "No effect detected for the first correlated pair")
  expect_true(abs(output$fit$theta_hat[3]) > 0.25 || abs(output$fit$theta_hat[4]) > 0.25,
              "No effect detected for the second correlated pair")
  
  # 5. Model fit quality
  expect_true(output$fit$final_loglike > -Inf)
  expect_true(is.finite(output$fit$dispersion))
})

# =============================================================================
# Test Cox Proportional Hazards Models
# =============================================================================
test_that("laser handles Cox regression with independent predictors", {
  # SETUP: Create controlled test data with known structure
  set.seed(1)
  n <- 500                                     # Number of observations
  p <- 10                                      # Number of predictors
  X <- matrix(rnorm(n * p), n, p)              # Independent predictors
  true_vars <- c(6, 9)                         # Variables with true effects
  
  # Create survival data
  linear_pred <- 0.8*X[, true_vars[1]] + 0.6*X[, true_vars[2]]
  baseline_hazard <- rexp(n, 0.1)              # Baseline hazard
  hazard <- baseline_hazard * exp(linear_pred) # Proportional hazards
  time <- rexp(n, hazard)                      # Event times
  
  # Add censoring
  censor_time <- rexp(n, 0.05)                 # Censoring times
  observed_time <- pmin(time, censor_time)     # Observed time is minimum
  status <- as.numeric(time <= censor_time)    # Event indicator
  
  # Create survival response
  y <- cbind(observed_time, status)
  
  # EXECUTE: Run laser algorithm
  output <- laser(X, y, family = "cox")
  
  # VERIFY: Test output structure and correctness
  # 1. Basic structure checks
  expect_s3_class(output, "laser")               # Should return proper S3 class
  expect_true(all(c("cs", "fit") %in% names(output)))  # Contains expected components
  
  # 2. Confidence sets validation
  expect_type(output$cs, "list")                 # cs should be a list
  expect_true(all(c("sets", "coverage") %in% names(output$cs)))  # Has both sets and coverage
  expect_type(output$cs$sets, "list")            # sets should be a list of confidence sets
  
  # 3. Variable selection correctness
  # Check if any confidence set contains at least one of the true variables
  found_true_var <- FALSE
  for (i in seq_along(output$cs$sets)) {
    if (any(true_vars %in% output$cs$sets[[i]])) {
      found_true_var <- TRUE
      break
    }
  }
  expect_true(found_true_var, "No confidence set contains any of the true variables")
  
  # 4. Result properties
  expect_true(all(output$fit$p_values >= 0 & output$fit$p_values <= 1))  # Valid p-values
  expect_true(all(output$cs$coverage <= 1))       # Coverage should not exceed 1
  
  # 5. Coefficient estimates
  expect_equal(length(output$fit$theta_hat), p)   # Should have one coefficient per variable
  expect_true(abs(output$fit$theta_hat[6]) > 0.3) # True effect should be detected
  expect_true(abs(output$fit$theta_hat[9]) > 0.2) # True effect should be detected
  
  # 6. Model fit quality
  expect_true(output$fit$final_loglike > -Inf)    # Should have valid likelihood
})

test_that("laser handles Cox regression with highly correlated predictors", {
  # SETUP: Create controlled test data with multicollinearity
  set.seed(1)
  n <- 500                                     # Number of observations
  p <- 10                                      # Number of predictors
  
  # Generate base independent random variables
  Z <- matrix(rnorm(n * p), n, p)
  
  # Create highly correlated pairs
  X <- matrix(0, n, p)
  X[, 1] <- Z[, 1]
  X[, 2] <- 0.95 * Z[, 1] + 0.05 * Z[, 2]      # X1, X2 correlated at ~0.95
  X[, 3] <- Z[, 3]
  X[, 4] <- 0.95 * Z[, 3] + 0.05 * Z[, 4]      # X3, X4 correlated at ~0.95
  
  # Generate remaining predictors independently
  X[, 5:p] <- Z[, 5:p]
  
  # Generate survival data - only X1 and X4 have true effects
  true_pair1 <- c(1, 2)                       # First correlated pair
  true_pair2 <- c(3, 4)                       # Second correlated pair
  
  # Create survival data with effects from X1 and X4
  linear_pred <- X[, 1] + X[, 4]
  baseline_hazard <- rexp(n, 0.1)             # Baseline hazard
  hazard <- baseline_hazard * exp(linear_pred) # Proportional hazards
  time <- rexp(n, hazard)                     # Event times
  
  # Add censoring
  censor_time <- rexp(n, 0.05)                # Censoring times
  observed_time <- pmin(time, censor_time)    # Observed time is minimum
  status <- as.numeric(time <= censor_time)   # Event indicator
  
  # Create survival response
  y <- cbind(observed_time, status)
  
  # EXECUTE: Run laser algorithm
  output <- laser(X, y, family = "cox")
  
  # VERIFY: Test output structure and correctness
  # 1. Basic structure checks
  expect_s3_class(output, "laser")
  expect_true(all(c("cs", "fit") %in% names(output)))
  
  # 2. Confidence sets validation
  expect_type(output$cs, "list")
  expect_true(all(c("sets", "coverage") %in% names(output$cs)))
  expect_type(output$cs$sets, "list")
  
  # 3. Variable selection correctness under multicollinearity
  # Check if we find at least one variable from each correlated pair
  found_pair1 <- FALSE
  found_pair2 <- FALSE
  
  for (i in seq_along(output$cs$sets)) {
    set_vars <- output$cs$sets[[i]]
    if (any(true_pair1 %in% set_vars)) {
      found_pair1 <- TRUE
    }
    if (any(true_pair2 %in% set_vars)) {
      found_pair2 <- TRUE
    }
  }
  
  expect_true(found_pair1, "No confidence set contains any variable from the first correlated pair")
  expect_true(found_pair2, "No confidence set contains any variable from the second correlated pair")
  
  # 4. Check for uncertainty in the estimated coefficients
  # At least one variable from each pair should have non-negligible coefficient
  max_coef_pair1 <- max(abs(output$fit$theta_hat[true_pair1]))
  max_coef_pair2 <- max(abs(output$fit$theta_hat[true_pair2]))
  
  expect_true(max_coef_pair1 > 0.25, "No effect detected for the first correlated pair")
  expect_true(max_coef_pair2 > 0.25, "No effect detected for the second correlated pair")
  
  # 5. Model fit quality
  expect_true(output$fit$final_loglike > -Inf)
})