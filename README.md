# glmcs: Generalized Linear Models with Confident Sets

[![GitHub stars](https://img.shields.io/github/stars/yizenglistat/glmcs.svg)](https://github.com/yizenglistat/glmcs/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/yizenglistat/glmcs.svg)](https://github.com/yizenglistat/glmcs/network)[![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue.svg)](https://www.r-project.org/)[![CRAN status](https://www.r-pkg.org/badges/version/glmcs)](https://CRAN.R-project.org/package=glmcs)[![Downloads](https://cranlogs.r-pkg.org/badges/glmcs)](https://cran.r-project.org/package=glmcs)
[![Codecov test coverage](https://codecov.io/gh/yizenglistat/glmcs/branch/main/graph/badge.svg)](https://codecov.io/gh/yizenglistat/glmcs)[![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

## Overview

**glmcs** implements *likelihood-based additive single-effect regression (LASER)* via a block-coordinate ascent algorithm, and provides **confidence sets (CSs)** for variable selection with uncertainty quantification.  

Key features:
- **Robust variable selection** under high multicollinearity  
- **Confident sets** that capture posterior uncertainty in variable inclusion  
- **Flexible modeling**: linear, generalized linear (e.g. logistic, Poisson), and Cox regression  
- **High/Scalable performance** via a C++ (Rcpp & Armadillo) backend  
- **Configurable**: tune number of effects (`L`), coverage, convergence tolerance, update strategy  

## Installation

Development version from GitHub:

```r
# devtools
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("yizenglistat/glmcs")
```

## Demo

```r
library(glmcs)

SEED <- 42

set.seed(SEED)

n <- 500; p <- 5

# Create correlation matrix with specific structure
cor_matrix <- matrix(0, nrow = p, ncol = p)
diag(cor_matrix) <- 1

# X1, X2, X3 highly correlated
cor_matrix[1:3, 1:3] <- 0.95
diag(cor_matrix[1:3, 1:3]) <- 1

# X4, X5 highly correlated
cor_matrix[4:5, 4:5] <- 0.95
diag(cor_matrix[4:5, 4:5]) <- 1

# Generate multivariate normal data
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)

# coefficents
theta <- rep(0, p); 
theta[c(2, 5)] <- 1
prob <- plogis(0.5 + X %*% theta)
y <- rbinom(n, 1, prob)

res <- glmcs(X           = X, 
             y           = y,
             family      = binomial("logit"),
             L           = 10L,
             coverage    = 0.95,
             min_abs_corr= 0.,
             seed        = SEED
)

print(res$cs)
# $sets
# $sets[[1]]
# [1] 5

# $sets[[2]]
# [1] 1 2


# $coverage
# [1] 0.9947149 0.9872938
```

![demo](https://github.com/user-attachments/assets/80e7d259-78ff-46ee-b820-c22796026bc2)


## Main Functions

- `glmcs(X, y, family, L, coverage, ...)`  
  End‐to‐end wrapper: fits the LASER model, variable selection, builds confidence sets.

- `get_laser_fit()`  
  Core fitting engine (C++): returns raw posterior model probabilities, point estimates, p-values, etc.

- `get_included()`, `get_cs()`  
  Select significant single‐effect components and construct credible sets with purity filtering (if applicable).

- `get_glm_fit()`, `get_cox_fit()`, `null_glm_fit()`, `null_cox_fit()`  
  Fast univariate fits (Gaussian, Binomial, Poisson, Gamma, inverse-Gaussian, Cox) via closed-form or IRLS.

- Low-level utilities:  
  - `get_loglike()`, `update_intercept()`, `update_dispersion()`, `get_purity()`, etc.

## Documentation & Support

- **Reference manual**: [CRAN](https://CRAN.R-project.org/package=glmcs)  
- **Issues**: [![GitHub issues](https://img.shields.io/github/issues-raw/yizenglistat/glmcs.svg)](https://github.com/yizenglistat/glmcs/issues)

## License

[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) | © 2025 Yizeng Li & Wei Pan

---

> _“Extensible and generalizable confidence‐driven variable selection in high‐dimensional, multicollinear settings---glmcs delivers both statistical rigor and computational efficiency.”_  