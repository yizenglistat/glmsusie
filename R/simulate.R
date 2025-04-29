#' Simulate Data for Generalized Linear, or Cox Models
#'
#' @description 
#' Generate synthetic data with controlled correlation structure for benchmarking 
#' variable selection methods.
#'
#' @param n Integer. Number of observations to generate (default: 500).
#' @param p Integer. Number of predictors to generate (default: 10).
#' @param intercept Numeric. Intercept term (default: -1, ignored for Cox models).
#' @param dispersion Numeric. Dispersion parameter for two-parameter GLMs (default: 9).
#' @param family Family object or "cox" string. The family object provides the link function.
#' @param settings Character. Correlation structure (default: "example-1").
#' @param cor Numeric. Correlation coefficient for "example-1" setting (default: 0.9).
#' @param censoring_rate Numeric. Censoring rate for survival data (default: 0.3).
#' @param seed Integer. Random seed for reproducibility (default: NULL).
#'
#' @details
#' Generates synthetic datasets with specific correlation structures:
#' 
#' - "example-1": Two highly correlated predictors, with only one truly associated
#'   with the response.
#'   
#' - "example-2": Block correlation structure with 5 predictors in two correlated groups.
#'
#' @return A list containing:
#' \itemize{
#'   \item X: Design matrix (n × p)
#'   \item y: Response variable (vector for GLMs, 2-column matrix for Cox)
#'   \item intercept: Intercept used in data generation
#'   \item dispersion: Dispersion parameter used in data generation
#'   \item family: Family object or "cox" string
#'   \item theta: True coefficient vector
#'   \item settings: Correlation setting used
#' }
#'
#' @examples
#' # Generate Gaussian data with correlated predictors
#' sim_data <- simulate_data(n=200, p=10, family=gaussian())
#' 
#' # Generate binary data for logistic regression
#' sim_data <- simulate_data(n=500, p=5, family=binomial(), settings="example-2")
#' 
#' # Generate survival data for Cox model with 40% censoring
#' sim_data <- simulate_data(n=300, p=10, family="cox", censoring_rate=0.4)
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
#' @importFrom stats rgamma rbinom rpois rnorm runif
#'
#' @export
simulate_data <- function(n = 500, 
                          p = 10, 
                          intercept = -1, 
                          dispersion = 9, 
                          family = gaussian(), 
                          settings = "example-1",
                          cor = 0.9,
                          censoring_rate = 0.3,
                          seed = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Check dimensional requirements for each setting
  if (settings == "example-1" && p < 2) {
    stop("For 'example-1' setting, p must be at least 2.")
  }
  if (settings == "example-2" && p < 5) {
    stop("For 'example-2' setting, p must be at least 5.")
  }

  # Extract family information
  if (is.character(family) && family == "cox") {
    family_name <- "cox"
    family_link <- NULL
  } else {
    family_name <- family$family
    family_link <- family$link
  }
  
  # Initialize X matrix and coefficient vector
  X <- matrix(0, nrow = n, ncol = p)
  theta <- rep(0, p)
  
  # Generate correlated predictors based on settings
  if (settings == "example-1") {
    # Two highly correlated predictors (X1, X2)
    Sigma <- matrix(c(1, cor, cor, 1), nrow = 2)
    X12 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
    
    # Add the two correlated predictors to X
    X[, 1:2] <- X12
    
    # Add remaining predictors as independent standard normal
    if (p > 2) {
      X[, 3:p] <- matrix(rnorm(n * (p - 2)), nrow = n)
    }
    
    # Set true coefficients (only X1 has effect)
    theta[1] <- 1
    
  } else if (settings == "example-2") {
    # Block correlation structure with 5 main predictors
    S <- rbind(
      c(1.0,  0.92, 0.7,  0.7,  0.9),
      c(0.92, 1.0,  0.7,  0.7,  0.7),
      c(0.7,  0.7,  1.0,  0.92, 0.8),
      c(0.7,  0.7,  0.92, 1.0,  0.8),
      c(0.9,  0.7,  0.8,  0.8,  1.0)
    )
    
    # Ensure positive definiteness
    S <- as.matrix(Matrix::nearPD(S)$mat)
    
    # Generate correlated predictors
    X_core <- MASS::mvrnorm(n, mu = rep(0, 5), Sigma = S)
    
    # Add to X matrix
    X[, 1:5] <- X_core
    
    # Add remaining predictors as independent
    if (p > 5) {
      X[, 6:p] <- matrix(rnorm(n * (p - 5)), nrow = n)
    }
    
    # Set true coefficients (X2 and X3 have effects)
    theta[2] <- 1
    theta[3] <- 1
    
  } else {
    stop("Unknown settings parameter. Use 'example-1' or 'example-2'.")
  }
  
  # Calculate linear predictor
  eta <- intercept + X %*% theta
  
  # Generate response variable based on family
  if (family_name == "cox") {
    # Cox model (survival data)
    intercept <- 0  # No intercept in Cox model
    dispersion <- 1  # Not used for Cox
    
    # Generate survival times from exponential distribution
    # with hazard proportional to exp(eta)
    lambda <- exp(eta)
    true_times <- rexp(n, rate = lambda)
    
    # Add censoring based on specified rate
    C <- runif(n, 0, quantile(true_times, 1 - censoring_rate))
    times <- pmin(true_times, C)
    status <- as.numeric(true_times <= C)
    
    # Create survival response (time, status)
    y <- cbind(times, status)
    colnames(y) <- c("time", "status")
    
  } else {
    # For GLM models, use the link function from the family object
    
    # Get inverse link function (for mean calculation)
    inv_link <- family$linkinv
    
    # Calculate mean using the inverse link function
    mu <- inv_link(eta)
    
    # Generate based on family type
    if (family_name == "gaussian") {
      # Linear model (Gaussian)
      y <- as.numeric(mu + rnorm(n, 0, sqrt(dispersion)))
      
    } else if (family_name == "binomial") {
      # Logistic/probit/etc. regression based on link
      dispersion <- 1  # Not used for binomial
      y <- rbinom(n, 1, mu)
      
    } else if (family_name == "poisson") {
      # Poisson regression
      dispersion <- 1  # Not used for Poisson
      y <- rpois(n, mu)
      
    } else if (family_name == "Gamma") {
      # Gamma regression
      shape <- 1/dispersion  # Shape parameter
      scale <- mu * dispersion  # Scale parameter
      y <- rgamma(n, shape = shape, scale = scale)
      
    } else if (family_name == "inverse.gaussian") {
      # Inverse Gaussian regression
      lambda <- 1/dispersion  # Precision parameter
      z <- rnorm(n)
      y <- mu + (mu^2 * z^2) / (2 * lambda) - 
        mu * sqrt(mu * z^2 / lambda) / 2
      
    } else {
      # Generic case for other families
      warning("Using approximation for unsupported family: ", family_name)
      
      # For other families, use a normal approximation with variance function
      var_func <- family$variance
      y <- rnorm(n, mean = mu, sd = sqrt(var_func(mu) * dispersion))
    }
  }
  
  # Return results
  return(list(
    X = X,
    y = y,
    intercept = intercept,
    dispersion = dispersion,
    family = family,
    theta = theta,
    settings = settings
  ))
}