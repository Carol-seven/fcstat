# Ledoit-Wolf shrinkage estimator

Compute the Ledoit-Wolf shrinkage estimator for the
covariance/correlation matrix.

## Usage

``` r
ledoit_wolf_est(X, method = "linshrink", res = "cov")
```

## Arguments

- X:

  A data matrix.

- method:

  A character string (default = "linshrink") specifying the method used
  in shrinkage, includes:

  1.  "linshrink": Linear shrinkage (Ledoit and Wolf 2004) .

  2.  "nlshrink": Non-linear shrinkage (Ledoit and Wolf 2015; Ledoit and
      Wolf 2017) .

  See
  [`linshrink_cov`](https://rdrr.io/pkg/nlshrink/man/linshrink_cov.html)
  and
  [`nlshrink_cov`](https://rdrr.io/pkg/nlshrink/man/nlshrink_cov.html)
  for details.

- res:

  A character string (default = "cov") specifying the result matrix to
  be obtained, either the covariance matrix ("cov") or the correlation
  matrix ("cor").

## Value

A numeric matrix.

## References

Ledoit O, Wolf M (2004). “A Well-Conditioned Estimator for
Large-Dimensional Covariance Matrices.” *Journal of Multivariate
Analysis*, **88**(2), 365–411.
[doi:10.1016/S0047-259X(03)00096-4](https://doi.org/10.1016/S0047-259X%2803%2900096-4)
.  
  
Ledoit O, Wolf M (2015). “Spectrum Estimation: A Unified Framework for
Covariance Matrix Estimation and PCA in Large Dimensions.” *Journal of
Multivariate Analysis*, **139**, 360–384.
[doi:10.1016/j.jmva.2015.04.006](https://doi.org/10.1016/j.jmva.2015.04.006)
.  
  
Ledoit O, Wolf M (2017). “Numerical Implementation of the QuEST
Function.” *Computational Statistics & Data Analysis*, **115**, 199–223.
[doi:10.1016/j.csda.2017.06.004](https://doi.org/10.1016/j.csda.2017.06.004)
.
