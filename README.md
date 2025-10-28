# glmsusie: A simple implementation of variable selection for general regression models with highly correlated predictors

<!-- [![GitHub stars](https://img.shields.io/github/stars/yizenglistat/glmsusie.svg)](https://github.com/yizenglistat/glmsusie/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/yizenglistat/glmsusie.svg)](https://github.com/yizenglistat/glmsusie/network)
 -->
[![R-CMD-check](https://github.com/yizenglistat/glmsusie/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yizenglistat/glmsusie/actions/workflows/R-CMD-check.yaml)
[![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue.svg)](https://www.r-project.org/)
[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![CRAN status](https://www.r-pkg.org/badges/version/glmsusie)](https://CRAN.R-project.org/package=glmsusie)
[![Downloads](https://cranlogs.r-pkg.org/badges/glmsusie)](https://cran.r-project.org/package=glmsusie)
[![Codecov test coverage](https://codecov.io/gh/yizenglistat/glmsusie/graph/badge.svg)](https://app.codecov.io/gh/yizenglistat/glmsusie)


## Overview

**glmsusie** implements the *generalized sum of single effects (gSuSiE)* framework via a blockwise coordinate ascent algorithm, performing variable selection with uncertainty quantification for generalized linear and Cox models. 

This R package accompanies the manuscript entitled *"A simple implementation of variable selection for general regression models with highly correlated predictors"* and provides reproducible code for the simulation studies presented in the paper.

Key features:
- **Variable selection** with uncertainty quantification through BIC-based posterior inclusion probabilities (PIPs)
- **Credible sets** that capture selection uncertainty in variable inclusion
- **Flexible modeling**: Generalized linear models (e.g., Gaussian, logistic, Poisson, and Gamma), Cox proportional hazards regression and general regression models
- **Scalable performance** via optimized C++ implementation using Rcpp and Armadillo
- **Handles multicollinearity** by modeling correlated predictors through sum of single-effect components

## Installation

Development version from GitHub:

```r
# remotes
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("yizenglistat/glmsusie")
```

## Demo

```r
library(glmsusie)

SEED <- 42

set.seed(SEED)

n <- 500; p <- 5

# Create correlation matrix with specific structure
cor_matrix <- matrix(0, nrow = p, ncol = p)
diag(cor_matrix) <- 1

# X1, X2, X3 highly correlated
cor_matrix[1:3, 1:3] <- 0.98
diag(cor_matrix[1:3, 1:3]) <- 1

# X4, X5 highly correlated
cor_matrix[4:5, 4:5] <- 0.98
diag(cor_matrix[4:5, 4:5]) <- 1

# Generate multivariate normal data
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)

# coefficents
theta <- rep(0, p); 
theta[c(2, 5)] <- 1
prob <- plogis(0.5 + X %*% theta)
y <- rbinom(n, 1, prob)

fit <- glmsusie(X           = X, 
                y           = y,
                family      = binomial("logit"),
                L           = 10L,
                coverage    = 0.95,
                seed        = SEED)
summary(fit)
## Call:
## glmsusie(X = X, y = y, L = 10L, family = binomial("logit"), coverage = 0.95,
##     seed = SEED)

## Family: binomial

## Coefficients: (sorted by magnitude)
##    Estimate    PIP
## X5   0.8867 0.9128
## X3   0.6677 0.7747
## X2   0.1026 0.1224
## X1   0.0861 0.1028
## X4   0.0825 0.0872

## 95% Confidence Sets:
##           Set Coverage
## cs1    {4, 5}        1
## cs2 {1, 2, 3}        1

## Model converged after 3 iterations.
## Computation time: 0.01 seconds.
```

## To reproduce the simulation results in the paper

We can reproduce the results by setting different sample size `n` and `family` accordingly below. 

#### S1 simulation results reproduce
```{r}
library(glmsusie)
seed <- 42
n <- 50
family <- gaussian()
# family <- binomial()
# family <- "cox"

# -----------------------------------------------------------------------------
# Experiment 1: SuSiE "S1" setting (two covariates highly correlated, p=2)
# True signal: (1, 0)
# -----------------------------------------------------------------------------
res_s1 <- benchmark(
  settings    = "S1",           # correlation structure
  n_sims      = 1000,           # number of Monte Carlo replicates
  n           = n,              # sample size
  rho         = 0.98,           # within-block correlation
  family      = family,         # family distribution
  true_theta  = c(1, 0),        # true coefficients
  intercept   = 0,              # true intercept
  dispersion  = 9,              # error variance
  parallel    = TRUE,           # parallel computing
  seed        = seed            # reproducibility
)

res_s1$glmsusie
res_s1$susie
res_s1$lasso
res_s1$enet
```

#### S2 simulation results reproduce
```{r}
library(glmsusie)
seed <- 42
n <- 50
family <- gaussian()
# family <- binomial()
# family <- "cox"

# -----------------------------------------------------------------------------
# Experiment 2: SuSiE "S2" setting (5 variables in a 5×5 block, p=5)
# True signal: (0, 1, 1, 0, 0)
# -----------------------------------------------------------------------------
res_s2 <- benchmark(
  settings    = "S2",
  n_sims      = 1000,
  n           = n,
  rho         = 0.9,
  family      = family,
  true_theta  = c(0, 1, 1, 0, 0),
  intercept   = 0,
  dispersion  = 9,
  parallel    = TRUE,
  seed        = seed
)

# Compare the performance of gSuSiE vs SuSiE in this block-correlation scenario
res_s2$glmsusie
res_s2$susie
res_s2$lasso
res_s2$enet
```

#### S3 simulation results reproduce
```{r}
library(glmsusie)
seed <- 42
n <- 50
family <- gaussian()
# family <- binomial()
# family <- "cox"

# -----------------------------------------------------------------------------
# Experiment 3: Additional "S3" setting (two blocks of size p/2 × p/2, p=4)
# True signal: (0, 1, 0, 1)
# -----------------------------------------------------------------------------
res_s3 <- benchmark(
  settings    = "S3",
  n_sims      = 1000,
  n           = n,
  rho         = 0.98,
  family      = family,
  true_theta  = c(0, 1, 0, 1),
  intercept   = 0,
  dispersion  = 9,
  parallel    = TRUE,
  seed        = seed
)

res_s3$glmsusie
res_s3$susie
res_s3$lasso
res_s3$enet
```

## Open Issues & Support

[![GitHub issues](https://img.shields.io/github/issues-raw/yizenglistat/glmsusie.svg)](https://github.com/yizenglistat/glmsusie/issues)

## License

GPL-3 | © 2025 Yizeng Li & Wei Pan

---

> _“Extensible and generalizable variable selection under strong multicollinear settings---glmsusie delivers both statistical rigor and computational efficiency.”_  