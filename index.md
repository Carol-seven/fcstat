# spice

## Sparse Precision (Inverse Covariance) Estimation

[![GitHub last
commit](https://img.shields.io/github/last-commit/Carol-seven/fcstat)](https://github.com/Carol-seven/fcstat/commits/master)
[![R-CMD-check](https://github.com/Carol-seven/fcstat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Carol-seven/fcstat/actions/workflows/R-CMD-check.yaml)
[![GitHub
License](https://img.shields.io/github/license/Carol-seven/fcstat?color=blue)](https://github.com/Carol-seven/fcstat/blob/master/LICENSE.md)

The goal of **spice** is to provide classical statistical methods for
estimating sparse precision (inverse covariance) matrix for functional
connectivity analysis in brain networks, making these methods accessible
and easy to use for researchers and practitioners in neuroimaging.

## Installation

You can install the development version of **spice** from
[GitHub](https://github.com/Carol-seven/spice) with:

``` r
# install.packages("devtools")
devtools::install_github("Carol-seven/spice")
```

## Example

``` r
library(spice)

set.seed(123)

X <- matrix(rnorm(200), 10, 20)

## Statistical methods for estimating the precision matrix,
## including the estimation and selection process
spice(X, method = "glasso", pkg = "glasso", crit = "CV", fold = 5)
```
