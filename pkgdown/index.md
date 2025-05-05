# glmcs

## About

The `glmcs` package provides an innovative approach to variable selection in generalized linear models and other regression frameworks, including Cox regression for survival analysis. Our method excels when predictors are highly correlated and the true model is sparse.

At its core, `glmcs` implements the "Likelihood-based Additive Single-Effect Regression" (LASER) model – a novel sparse regression approach that essentially enhances traditional forward selection. The algorithm produces "Confidence Sets" (CSs), which are groups of correlated variables with high probability of containing at least one non-zero effect.

These confidence sets offer a statistical guarantee about variable importance while acknowledging the inherent uncertainty in highly correlated predictor spaces. Rather than forcing selection of a single variable that might be unstable, the CS approach identifies groups where you can be confident at least one variable is influential.

Developed by Yizeng Li and Wei Pan from the Division of Biostatistics at the University of Minnesota Twin Cities.


## Quick start

You can install the released version of glmcs from CRAN:

```r
install.packages("glmcs")
```

Alternatively, install the latest development version of `glmcs` from Github:

```r
# install.packages("remotes")
remotes::install_github("yizenglistat/glmcs")
```

See [here](articles/mwe.html) for a brief illustraton of `glmcs`. For more documentation and examples please visit [https://github.com/yizenglistat/glmcs](https://github.com/yizenglistat/glmcs).