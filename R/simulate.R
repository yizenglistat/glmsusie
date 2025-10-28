#' @title Generate Synthetic Data for Benchmark Simulations
#'
#' @description
#' Creates a design matrix \code{X} with controlled correlation structure and simulates
#' a response \code{y} according to the specified \code{family}.  Supports Gaussian,
#' Binomial, Poisson, Gamma GLMs, and Cox survival data with exponential baseline hazard.
#'
#' @param n Integer; number of observations (default: 600).
#' @param theta Numeric vector of true coefficients (length p).
#' @param intercept Numeric; true intercept term (default: 0).
#' @param settings Character; correlation pattern:
#'   \itemize{
#'     \item \code{"S1"}: all p variables equally correlated.
#'     \item \code{"S2"}: first 5×5 block correlated, others independent.
#'     \item \code{"S3"}: two p/2 × p/2 blocks correlated.
#'   }  Default: \code{c("S1","S2","S3")}.
#' @param rho Numeric in \code{[0,1]}; within-block correlation (default: 0.9).
#' @param dispersion Numeric; dispersion for Gaussian/Gamma (default: 1).
#' @param family A \code{\link[stats]{family}} object (e.g. \code{\link[stats]{gaussian}()}, \code{\link[stats]{binomial}()},
#'   \code{\link[stats]{poisson}()}, \code{\link[stats]{Gamma}()}) or the string \code{"cox"} for survival data
#'   (default: \code{\link[stats]{gaussian}()}).
#' @param censoring_rate Numeric in \code{[0,1]}; fraction censored for Cox models
#'   (default: 0.3).
#' @param baseline_hazard Numeric or function; constant hazard rate for Cox, or a
#'   hazard function if provided (default: 1).
#'
#' @return A list with components:
#' \describe{
#'   \item{X}{n × p design matrix with specified correlation.}
#'   \item{y}{Response:
#'     \itemize{
#'       \item Numeric vector of length n (GLMs).
#'       \item n × 2 matrix (time, status) (Cox).
#'     }}
#'   \item{eta}{Linear predictor intercept + X %*% theta.}
#' }
#'
#' @examples
#' \dontrun{
#' # Gaussian data, p = 2
#' dat1 <- generate(n = 100, theta = c(1, 0), settings = "S1",
#'                  family = gaussian(), dispersion = 2)
#'
#' # Binomial data with probit link
#' dat2 <- generate(n = 200, theta = c(0.5, -0.5),
#'                  family = binomial(link = "probit"),
#'                  settings = "S2", rho = 0.7)
#'
#' # Cox survival data
#' dat3 <- generate(n = 150, theta = c(1, 1, 0),
#'                  settings = "S3", family = "cox",
#'                  censoring_rate = 0.2, baseline_hazard = 0.05)
#' }
#'
#' @export
generate <- function(
  n = 600,
  theta = c(1, 0),
  intercept = 0,
  settings = c("S1", "S2", "S3"),
  rho = 0.9,
  dispersion = 1,
  family = gaussian(),
  censoring_rate = 0.3,
  baseline_hazard = 1
) {
  settings <- match.arg(settings)
  p        <- length(theta)

  ## 1) construct correlation matrix Σ
  Sigma <- switch(settings,
    S1 = { M <- matrix(rho, p, p); diag(M) <- 1; M },
    S2 = {
      if (p < 5) stop("S2 requires length(theta) >= 5")
      M5 <- matrix(rho, 5, 5); diag(M5) <- 1
      M  <- diag(1, p); M[1:5,1:5] <- M5; M
    },
    S3 = {
      if (p %% 2 != 0) stop("S3 requires length(theta) even")
      b <- p/2
      B <- matrix(rho, b, b); diag(B) <- 1
      M <- diag(1, p)
      M[1:b,1:b]       <- B
      M[(b+1):p,(b+1):p] <- B
      M
    }
  )

  ## 2) sample X ~ MVN(0, Σ)
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

  ## 3) compute linear predictor
  eta <- as.numeric(intercept + X %*% theta)

  fam_name <- if (is.character(family)) family else family$family

  ## 4) simulate response y
  if (tolower(fam_name) == "cox") {
    rate    <- baseline_hazard * exp(eta)
    T_event <- rexp(n, rate = rate)
    if      (censoring_rate <= 0) C_time <- rep(Inf, n)
    else if (censoring_rate >= 1) C_time <- rep(0, n)
    else {
      medT  <- median(T_event)
      crate <- -log(1 - censoring_rate) / medT
      C_time <- rexp(n, rate = crate)
    }
    time   <- pmin(T_event, C_time)
    status <- as.integer(T_event <= C_time)
    y <- cbind(time = time, status = status)

  } else {
    mu <- family$linkinv(eta)
    y <- switch(fam_name,
      gaussian = as.numeric(mu + rnorm(n, sd = sqrt(dispersion))),
      binomial = rbinom(n, size = 1, prob = mu),
      poisson  = rpois(n, lambda = mu),
      Gamma    = {
        shape <- 1 / dispersion
        scale <- mu * dispersion
        rgamma(n, shape = shape, scale = scale)
      },
      stop("Unsupported family: ", fam_name)
    )
  }

  list(X = X, y = y, eta = eta)
}