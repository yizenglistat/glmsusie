# glmsusie: Generalized Linear Models with Confidence Sets

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

**glmsusie** implements *likelihood-based additive single-effect regression (LASER)* via a block-coordinate ascent algorithm, and provides *confidence sets (CSs)* for variable selection with uncertainty quantification.  

Key features:
- **Robust variable selection** under high multicollinearity  
- **Confident sets** that capture posterior uncertainty in variable inclusion  
- **Flexible modeling**: Gaussian, generalized linear (e.g. logistic, Poisson), and Cox regression  
- **High/Scalable performance** via a C++ (Rcpp & Armadillo) backend  
- **Extensions**: applicable to general likelihood-based regression variable selection problems

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


## Open Issues & Support

[![GitHub issues](https://img.shields.io/github/issues-raw/yizenglistat/glmsusie.svg)](https://github.com/yizenglistat/glmsusie/issues)

## License

GPL-3 | © 2025 Yizeng Li & Wei Pan

---

> _“Extensible and generalizable variable selection with confidence in high‐dimensional, multicollinear settings---glmsusie delivers both statistical rigor and computational efficiency.”_  