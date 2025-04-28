# glmcs: Scalable Additive‐Effect Regression with Credible Sets

[![CRAN status](https://www.r-pkg.org/badges/version/glmcs)](https://CRAN.R-project.org/package=glmcs)  
[![Build Status](https://github.com/yizenglistat/glmcs/actions/workflows/check.yaml/badge.svg)](https://github.com/yizenglistat/glmcs/actions)  
[![Codecov test coverage](https://codecov.io/gh/yizenglistat/glmcs/branch/main/graph/badge.svg)](https://codecov.io/gh/yizenglistat/glmcs)  

## Overview

`glmcs` implements **L-Effect Sparse Regression (LASER)** models:  
a family of additive‐effect regression methods for high-dimensional data with correlated predictors.  
Key features:

- **Univariate and multivariate** routines for linear, GLM (binomial, Poisson, Gamma, inverse-Gaussian) and Cox models.  
- **Block-coordinate ascent** algorithm with both “greedy” and “shuffle” updates.  
- **Credible sets** of variables with posterior mass and purity filtering to handle correlation.  
- **High performance** via C++ (Rcpp & Armadillo) backend.  
- **Easy tuning** of number of effects, coverage level, convergence tolerance, and more.  

## Installation

Stable CRAN release:

```r
install.packages("glmcs")
```

Development version from GitHub:

```r
# using pak
pak::pak("yizenglistat/glmcs")

# or with devtools
# devtools::install_github("yizenglistat/glmcs")
```

## Quick Start

```r
library(glmcs)

# Simulate correlated predictors and continuous response
set.seed(123)
n <- 200; p <- 50
X <- matrix(rnorm(n * p), n, p)
β_true <- rep(0, p); β_true[c(3, 7, 20)] <- c(2, -1.5, 1.2)
y <- X %*% β_true + rnorm(n)

# Fit LASER model
res <- laser(
  X, y,
  family      = gaussian("identity"),
  L           = 5L,
  coverage    = 0.90,
  standardize = TRUE,
  method      = "greedy",
  seed        = 42
)

# Inspect results
print(res$theta)        # Estimated effect vector
print(res$cs)           # Credible sets of variables
plot(res$pmp, type = "h")  # Posterior marginal inclusion probabilities
```

## Main Functions

- `laser()`  
  End‐to‐end wrapper: fits the LASER model, selects effects, builds credible sets, and returns a tidy result object.

- `get_laser_fit()`  
  Core fitting engine (C++): returns raw posterior matrices, timing, etc.

- `get_included()`, `get_cs()`  
  Select significant single‐effect components and construct credible sets with purity filtering.

- `get_glm_fit()`, `get_cox_fit()`, `null_glm_fit()`, `null_cox_fit()`  
  Fast univariate fits (Gaussian, Binomial, Poisson, Gamma, inverse-Gaussian, Cox) via closed-form or IRLS.

- Low-level utilities:  
  - `get_loglike()`, `update_intercept()`, `update_dispersion()`, `get_purity()`, `kl_divergence()`.

## Documentation & Support

- **Reference manual**: [CRAN](https://CRAN.R-project.org/package=glmcs)  
- **Vignettes**:  
  ```r
  browseVignettes("glmcs")
  ```  
- **Source code & issues**:  
  [GitHub: yizenglistat/glmcs](https://github.com/yizenglistat/glmcs)  

## License

GPL-3 | © 2025 Yizeng Li & Wei Pan

---

> _“Confident variable selection in high dimensions demands both statistical rigor and computational efficiency—`glmcs` delivers on both fronts.”_  