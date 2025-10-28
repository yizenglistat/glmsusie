## About

The `glmsusie` package implements the **generalized Sum of Single Effects (gSuSiE)** model, a likelihood-based extension of the Bayesian SuSiE framework for variable selection in generalized linear and survival (Cox) regression models. gSuSiE is designed for settings with highly correlated predictors and sparse true signals, where traditional penalized regression methods often struggle.

At its core, gSuSiE decomposes regression effects into additive “single-effect” components, enabling stable variable selection and interpretable uncertainty quantification through **credible sets (CSs)** — groups of correlated variables that collectively have high probability of containing at least one true effect. 

Beyond credible sets and posterior inclusion probabilities (PIPs), gSuSiE also provides **standard errors, $p$-values, and 2BIC statistics** as additional diagnostic measures, offering both Bayesian and frequentist perspectives on variable importance.

Developed by **Yizeng Li** and **Wei Pan**, Division of Biostatistics and Health data Science, University of Minnesota Twin Cities.


## Quick start

You can install the released version of glmsusie from CRAN:

```r
install.packages("glmsusie")
```

Alternatively, install the latest development version of `glmsusie` from Github:

```r
# install.packages("remotes")
remotes::install_github("yizenglistat/glmsusie")
```

See [here](articles/mwe.html) for a brief illustraton of `glmsusie`. For more documentation and examples please visit [https://github.com/yizenglistat/glmsusie](https://github.com/yizenglistat/glmsusie).