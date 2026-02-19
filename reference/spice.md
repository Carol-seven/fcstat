# Sparse Precision (Inverse Covariance) Estimation

Provide a collection of statistical methods to estimate a precision
matrix.

## Usage

``` r
spice(
  X,
  method,
  base = "cov",
  n = NULL,
  lambda = NULL,
  nlambda = 20,
  lambda.min.ratio = 0.01,
  gamma = NA,
  initial = "glasso",
  pkg = "glasso",
  crit = "CV",
  fold = 5,
  ebic.tuning = 0.5,
  cores = 1
)
```

## Arguments

- X:

  1.  An \\n \times p\\ data matrix with sample size \\n\\ and dimension
      \\p\\.

  2.  A \\p \times p\\ sample covariance/correlation matrix with
      dimension \\p\\.

- method:

  A character string specifying the statistical method for estimating
  precision matrix. Available options include:

  1.  "glasso": Graphical lasso (Friedman et al. 2008) .

  2.  "ridge": Graphical ridge (van Wieringen and Peeters 2016) .

  3.  "elnet": Graphical elastic net (Zou and Hastie 2005) .

  4.  "clime": Constrained L1-minimization for inverse (covariance)
      matrix estimation (Cai et al. 2011) .

  5.  "tiger": Tuning-insensitive graph estimation and regression (Liu
      and Wang 2017) , which is only applicable when `X` is the \\n
      \times p\\ data matrix.

  6.  "adapt": Adaptive lasso (Zou 2006; Fan et al. 2009) .

  7.  "atan": Arctangent type penalty (Wang and Zhu 2016) .

  8.  "exp": Exponential type penalty (Wang et al. 2018) .

  9.  "mcp": Minimax concave penalty (Zou 2006) .

  10. "scad": Smoothly clipped absolute deviation (Fan and Li 2001; Fan
      et al. 2009) .

- base:

  A character string (default = "cov") specifying the calculation base:

  1.  "cov": The covariance matrix.

  2.  "cor": The correlation matrix.

  This is only applicable when `X` is the \\n \times p\\ data matrix.

- n:

  An integer (default = NULL) specifying the sample size. This is only
  required when the input matrix `X` is a \\p \times p\\ sample
  covariance/correlation matrix with dimension \\p\\.

- lambda:

  A non-negative numeric vector specifying the grid for the
  regularization parameter. The default is `NULL`, which generates its
  own `lambda` sequence based on `nlambda` and `lambda.min.ratio`. For
  `method = "clime"` combined with `pkg = "clime"`, the `lambda`
  sequence is based on `nlambda`, `lambda.min` and `lambda.max`.

- nlambda:

  An integer (default = 20) specifying the number of `lambda` values to
  generate when `lambda = NULL`.

- lambda.min.ratio:

  A numeric value \> 0 (default = 0.01) specifying the fraction of the
  maximum `lambda` value \\\lambda\_{max}\\ to generate the minimum
  `lambda` \\\lambda\_{min}\\. If `lambda = NULL`, a `lambda` grid of
  length `nlambda` is automatically generated on a log scale, ranging
  from \\\lambda\_{max}\\ down to \\\lambda\_{min}\\.

- gamma:

  A numeric value specifying the additional parameter for the chosen
  `method`. Default values:

  1.  "elnet": A sequence from 0.1 to 0.9 with increments of 0.1

  2.  "adapt": 0.5

  3.  "atan": 0.005

  4.  "exp": 0.01

  5.  "mcp": 3

  6.  "scad": 3.7

- initial:

  A \\p \times p\\ matrix or a \\p \times p \times \mathrm{npara}\\ (the
  number of all combinations of `lambda` and `gamma`) array specifying
  the initial estimate for `method` set to `"atan"`, `"exp"`, `"scad"`,
  and `"mcp"`; or specifying \\\tilde{\Omega}\\ of the adaptive weight
  for `method = "adapt"`, calculated as
  \\\lvert\tilde{\omega}\_{ij}\rvert^{-\gamma}\\, where \\\tilde{\Omega}
  := (\tilde{\omega}\_{ij})\\. Some options are also offered when a
  character string is provided (default = "glasso"), including:

  - "glasso": Use the precision matrix estimate derived from the
    graphical lasso.

  - "invS": Use the inverse calculation base matrix if the matrix is
    invertible.

  - "linshrink": Use the precision matrix estimate derived from
    Ledoit-Wolf linear shrinkage estimator of the population covariance
    matrix (Ledoit and Wolf 2004) .

  - "nlshrink": Use the precision matrix estimate derived from
    Ledoit-Wolf non-linear shrinkage estimator of the population
    covariance matrix (Ledoit and Wolf 2015; Ledoit and Wolf 2017) .

- pkg:

  A character string specifying the package option to use. The available
  options depend on the chosen `method`:

  1.  For `method = "glasso"`:

      - "ADMMsigma": The function from
        [`ADMMsigma`](https://rdrr.io/pkg/ADMMsigma/man/ADMMsigma.html).

      - "CovTools": The function from
        [`PreEst.glasso`](https://rdrr.io/pkg/CovTools/man/PreEst.glasso.html).

      - "CVglasso": The function from
        [`CVglasso`](https://rdrr.io/pkg/CVglasso/man/CVglasso.html).

      - "Glarmadillo": The function from
        [`glarma`](https://rdrr.io/pkg/Glarmadillo/man/glarma.html).

      - "glasso": The function from
        [`glasso`](https://rdrr.io/pkg/glasso/man/glasso.html).

      - "GLassoElnetFast": The function from
        [gelnet](https://github.com/TobiasRuckstuhl/GLassoElnetFast).

      - "glassoFast": The function from
        [`glassoFast`](https://rdrr.io/pkg/glassoFast/man/glassoFast.html).

      - "huge": The function from
        [`huge.glasso`](https://rdrr.io/pkg/huge/man/huge.glasso.html).

  2.  For `method = "ridge"`:

      - "ADMMsigma": The function from
        [`RIDGEsigma`](https://rdrr.io/pkg/ADMMsigma/man/RIDGEsigma.html).

      - "GLassoElnetFast": The function from
        [gelnet](https://github.com/TobiasRuckstuhl/GLassoElnetFast).

      - "porridge": The function from
        [`ridgePgen`](https://rdrr.io/pkg/porridge/man/ridgePgen.html).

      - "rags2ridges": The function from
        [`ridgeP`](https://cfwp.github.io/rags2ridges/reference/ridgeP.html).

  3.  For `method = "elnet"`:

      - "ADMMsigma": The function from
        [`ADMMsigma`](https://rdrr.io/pkg/ADMMsigma/man/ADMMsigma.html).

      - "GLassoElnetFast": The function from
        [gelnet](https://github.com/TobiasRuckstuhl/GLassoElnetFast).

  4.  For `method = "clime"`:

      - "clime": The function from
        [`clime`](https://rdrr.io/pkg/clime/man/clime.html).

      - "flare": The function from
        [`sugm`](https://rdrr.io/pkg/flare/man/sugm.html).

  5.  For `method = "tiger"`:

      - "flare": The function from
        [`sugm`](https://rdrr.io/pkg/flare/man/sugm.html).

      - "huge": The function from
        [`huge.tiger`](https://rdrr.io/pkg/huge/man/huge.tiger.html).

  6.  For `method` set to `"adapt"`, `"atan"`, `"exp"`, `"scad"`, and
      `"mcp"`:

      - "glasso": The function from
        [`glasso`](https://rdrr.io/pkg/glasso/man/glasso.html).

      - "GLassoElnetFast": The function from
        [gelnet](https://github.com/TobiasRuckstuhl/GLassoElnetFast).

      - "glassoFast": The function from
        [`glassoFast`](https://rdrr.io/pkg/glassoFast/man/glassoFast.html).

- crit:

  A character string (default = "CV") specifying the parameter selection
  method to use. Available options include:

  1.  "AIC": Akaike information criterion (Akaike 1973) .

  2.  "BIC": Bayesian information criterion (Schwarz 1978) .

  3.  "EBIC": extended Bayesian information criterion (Foygel and
      Drton 2010) .

  4.  "HBIC": high dimensional Bayesian information criterion (Wang et
      al. 2013; Fan et al. 2017) .

  5.  "CV": k-fold cross validation with negative log-likelihood loss.

- fold:

  An integer (default = 5) specifying the number of folds used for
  `crit = "CV"`.

- ebic.tuning:

  A numeric value in \[0, 1\] (default = 0.5) specifying the tuning
  parameter to calculate for `crit = "EBIC"`.

- cores:

  An integer (default = 1) specifying the number of cores to use for
  parallel execution.

## Value

An object with S3 class `"spice"` containing the following components:

- hatOmega_opt:

  The estimated precision matrix.

- lambda_opt:

  The optimal regularization parameter.

- gamma_opt:

  The optimal hyperparameter.

- hatOmega:

  A list of estimated precision matrices for `lambda` grid and `gamma`
  grid.

- lambda:

  The actual lambda grid used in the program, corresponding to
  `hatOmega`.

- gamma:

  The actual gamma grid used in the program, corresponding to
  `hatOmega`.

- CV.loss:

  Matrix of CV losses, with rows for CV folds and columns for parameter
  combinations, when `crit = "CV"`.

- IC.score:

  The information criterion score for each parameter combination when
  `crit` is set to `"AIC"`, `"BIC"`, `"EBIC"`, or `"HBIC"`.

## Note

For `method = "tiger"`, the estimation process solely relies on the raw
\\n \times p\\ data `X` and does not utilize the argument `base`. This
argument is not applicable for `method = "tiger"` and will have no
effect if provided.

## References

Akaike H (1973). “Information Theory and an Extension of the Maximum
Likelihood Principle.” In Petrov BN, Csáki F (eds.), *Second
International Symposium on Information Theory*, 267–281. Akad\\emiai
Kiad\\o, Budapest, Hungary.  
  
Cai TT, Liu W, Luo X (2011). “A Constrained \\\ell\\1 Minimization
Approach to Sparse Precision Matrix Estimation.” *Journal of the
American Statistical Association*, **106**(494), 594–607.
[doi:10.1198/jasa.2011.tm10155](https://doi.org/10.1198/jasa.2011.tm10155)
.  
  
Fan J, Feng Y, Wu Y (2009). “Network Exploration via the Adaptive LASSO
and SCAD Penalties.” *The Annals of Applied Statistics*, **3**(2),
521–541. [doi:10.1214/08-aoas215](https://doi.org/10.1214/08-aoas215)
.  
  
Fan J, Li R (2001). “Variable Selection via Nonconcave Penalized
Likelihood and its Oracle Properties.” *Journal of the American
Statistical Association*, **96**(456), 1348–1360.
[doi:10.1198/016214501753382273](https://doi.org/10.1198/016214501753382273)
.  
  
Fan J, Liu H, Ning Y, Zou H (2017). “High Dimensional Semiparametric
Latent Graphical Model for Mixed Data.” *Journal of the Royal
Statistical Society Series B: Statistical Methodology*, **79**(2),
405–421. [doi:10.1111/rssb.12168](https://doi.org/10.1111/rssb.12168)
.  
  
Foygel R, Drton M (2010). “Extended Bayesian Information Criteria for
Gaussian Graphical Models.” In Lafferty J, Williams C, Shawe-Taylor J,
Zemel R, Culotta A (eds.), *Advances in Neural Information Processing
Systems 23 (NIPS 2010)*, 604–612.  
  
Friedman J, Hastie T, Tibshirani R (2008). “Sparse Inverse Covariance
Estimation with the Graphical Lasso.” *Biostatistics*, **9**(3),
432–441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)
.  
  
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
  
Liu H, Wang L (2017). “TIGER: A Tuning-Insensitive Approach for
Optimally Estimating Gaussian Graphical Models.” *Electronic Journal of
Statistics*, **11**(1), 241–294.
[doi:10.1214/16-EJS1195](https://doi.org/10.1214/16-EJS1195) .  
  
Schwarz G (1978). “Estimating the Dimension of a Model.” *The Annals of
Statistics*, **6**(2), 461–464.
[doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136) .  
  
van Wieringen WN, Peeters CFW (2016). “Ridge Estimation of Inverse
Covariance Matrices from High-Dimensional Data.” *Computational
Statistics & Data Analysis*, **103**, 284–303.
[doi:10.1016/j.csda.2016.05.012](https://doi.org/10.1016/j.csda.2016.05.012)
.  
  
Wang L, Kim Y, Li R (2013). “Calibrating Nonconvex Penalized Regression
in Ultra-High Dimension.” *The Annals of Statistics*, **41**(5),
2505–2536. [doi:10.1214/13-AOS1159](https://doi.org/10.1214/13-AOS1159)
.  
  
Wang Y, Fan Q, Zhu L (2018). “Variable Selection and Estimation using a
Continuous Approximation to the \\L_0\\ Penalty.” *Annals of the
Institute of Statistical Mathematics*, **70**(1), 191–214.
[doi:10.1007/s10463-016-0588-3](https://doi.org/10.1007/s10463-016-0588-3)
.  
  
Wang Y, Zhu L (2016). “Variable Selection and Parameter Estimation with
the Atan Regularization Method.” *Journal of Probability and
Statistics*, **2016**, 6495417.
[doi:10.1155/2016/6495417](https://doi.org/10.1155/2016/6495417) .  
  
Zou H (2006). “The Adaptive Lasso and Its Oracle Properties.” *Journal
of the American Statistical Association*, **101**(476), 1418–1429.
[doi:10.1198/016214506000000735](https://doi.org/10.1198/016214506000000735)
.  
  
Zou H, Hastie T (2005). “Regularization and Variable Selection via the
Elastic Net.” *Journal of the Royal Statistical Society Series B:
Statistical Methodology*, **67**(2), 301–320.
[doi:10.1111/j.1467-9868.2005.00527.x](https://doi.org/10.1111/j.1467-9868.2005.00527.x)
.
