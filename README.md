# glmcs: Generalized Linear models with Confidence Sets

[![CRAN status](https://www.r-pkg.org/badges/version/glmcs)](https://CRAN.R-project.org/package=glmcs)[![Codecov test coverage](https://codecov.io/gh/yizenglistat/glmcs/branch/main/graph/badge.svg)](https://codecov.io/gh/yizenglistat/glmcs)  

## Overview

**glmcs** implements *likelihood-based additive single-effect regression (LASER)* via a block-coordinate ascent algorithm, and provides **confident sets (CSs)** for variable selection with uncertainty quantification.  

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

## Quick Start

```r
library(glmcs)

# Simulate correlated predictors and continuous response
set.seed(42)
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
```

## Main Functions

- `glmcs(X, y, family, L, coverage, ...)`  
  End‐to‐end wrapper: fits the LASER model, variable selection, builds confident sets.

- `get_laser_fit()  
  Core fitting engine (C++): returns raw posterior matrices, timing, etc.

- `get_included()`, `get_cs()`  
  Select significant single‐effect components and construct credible sets with purity filtering.

- `get_glm_fit()`, `get_cox_fit()`, `null_glm_fit()`, `null_cox_fit()`  
  Fast univariate fits (Gaussian, Binomial, Poisson, Gamma, inverse-Gaussian, Cox) via closed-form or IRLS.

- Low-level utilities:  
  - `get_loglike()`, `update_intercept()`, `update_dispersion()`, `get_purity()`, etc.

## Documentation & Support

- **Reference manual**: [CRAN](https://CRAN.R-project.org/package=glmcs)  
- **Source code & issues**:  
  [GitHub: yizenglistat/glmcs](https://github.com/yizenglistat/glmcs)  

## License

GPL-3 | © 2025 Yizeng Li & Wei Pan

---

> _“Confident variable selection in high dimensions demands both statistical rigor and computational efficiency—`glmcs` delivers on both fronts.”_  