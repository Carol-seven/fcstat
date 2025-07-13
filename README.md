# fcstat <img src="man/figure/logo.png" align="right" alt="" width="150"/>


[![GitHub last commit](https://img.shields.io/github/last-commit/Carol-seven/fcstat)](https://github.com/Carol-seven/fcstat/commits/master)
[![R-CMD-check](https://github.com/Carol-seven/fcstat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Carol-seven/fcstat/actions/workflows/R-CMD-check.yaml)
[![GitHub License](https://img.shields.io/github/license/Carol-seven/fcstat?color=blue)](https://github.com/Carol-seven/fcstat/blob/master/LICENSE.md)


The goal of **fcstat** is to provide classical statistical methods for estimating
functional connectivity analysis in brain networks, making them user-friendly
and useful for researchers and practitioners in the field of neuroimaging.


## Installation

You can install the development version of **fcstat** from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Carol-seven/fcstat")
```


## Example

``` r
library(fcstat)

X <- matrix(rnorm(200), 10, 20)

## Statistical methods for estimating the precision matrix,
## including the estimation and selection process
fcstat(X, method = "glasso", pkgopt = "glasso", crit = "CV", fold = 5)
```
