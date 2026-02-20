
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spice <img src="man/figures/logo.png" align="right" alt="" width="150"/>

## Sparse Precision (Inverse Covariance) Estimation

[![GitHub last commit](https://img.shields.io/github/last-commit/Carol-seven/spice)](https://github.com/Carol-seven/spice/commits/master)
[![R-CMD-check](https://github.com/Carol-seven/spice/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Carol-seven/spice/actions/workflows/R-CMD-check.yaml)
[![GitHub License](https://img.shields.io/github/license/Carol-seven/spice?color=blue)](https://github.com/Carol-seven/spice/blob/master/LICENSE.md)

The goal of **spice** is to provide classical statistical methods for estimating
sparse precision (inverse covariance) matrix for functional connectivity
analysis in brain networks, making these methods accessible and easy to use for
researchers and practitioners in neuroimaging.

## Methods

| Method | Reference |
|:---|:---|
| Graphical lasso (`method = "lasso"`) | Friedman et al. ([2008](#ref-friedman2008sparse)) |
| Graphical ridge (`method = "ridge"`) | Wieringen and Peeters ([2016](#ref-vanwieringen2016ridge)) |
| Graphical elastic net (`method = "elnet"`) | Zou and Hastie ([2005](#ref-zou2005regularization)) |
| CLIME (`method = "clime"`) | Cai et al. ([2011](#ref-cai2011aconstrained)) |
| TIGER (`method = "tiger"`) | Liu and Wang ([2017](#ref-liu2017tiger)) |
| Graphical adaptive lasso (`method = "adapt"`) | Zou ([2006](#ref-zou2006adaptive)); Fan et al. ([2009](#ref-fan2009network)) |
| Arctangent type penalty (`method = "atan"`) | Wang and Zhu ([2016](#ref-wang2016variable)) |
| Exponential type penalty (`method = "exp"`) | Wang et al. ([2018](#ref-wang2018variable)) |
| MCP (`method = "mcp"`) | Zhang ([2010](#ref-zhang2010nearly)) |
| SCAD (`method = "scad"`) | Fan and Li ([2001](#ref-fan2001variable)); Fan et al. ([2009](#ref-fan2009network)) |

## Installation

You can install the development version of **spice** from
[GitHub](https://github.com/Carol-seven/spice) with:

``` r
# install.packages("devtools")
devtools::install_github("Carol-seven/spice")
```

## Example

    library(spice)

    set.seed(123)

    X <- matrix(rnorm(200), 10, 20)

    ## Statistical methods for estimating the precision matrix,
    ## including the estimation and selection process
    spice(X, method = "glasso", pkg = "glasso", crit = "CV", fold = 5)

## Reference

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-cai2011aconstrained" class="csl-entry">

Cai, Tony T., Weidong Liu, and Xi Luo. 2011. “A Constrained $\ell$<!-- -->1 Minimization Approach to Sparse Precision Matrix Estimation.” *Journal of the American Statistical Association* 106 (494): 594–607. <https://doi.org/10.1198/jasa.2011.tm10155>.

</div>

<div id="ref-fan2009network" class="csl-entry">

Fan, Jianqing, Yang Feng, and Yichao Wu. 2009. “Network Exploration via the Adaptive LASSO and SCAD Penalties.” *The Annals of Applied Statistics* 3 (2): 521–41. <https://doi.org/10.1214/08-aoas215>.

</div>

<div id="ref-fan2001variable" class="csl-entry">

Fan, Jianqing, and Runze Li. 2001. “Variable Selection via Nonconcave Penalized Likelihood and Its Oracle Properties.” *Journal of the American Statistical Association* 96 (456): 1348–60. <https://doi.org/10.1198/016214501753382273>.

</div>

<div id="ref-friedman2008sparse" class="csl-entry">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. “Sparse Inverse Covariance Estimation with the Graphical Lasso.” *Biostatistics* 9 (3): 432–41. <https://doi.org/10.1093/biostatistics/kxm045>.

</div>

<div id="ref-liu2017tiger" class="csl-entry">

Liu, Han, and Lie Wang. 2017. “TIGER: A Tuning-Insensitive Approach for Optimally Estimating Gaussian Graphical Models.” *Electronic Journal of Statistics* 11 (1): 241–94. <https://doi.org/10.1214/16-EJS1195>.

</div>

<div id="ref-wang2018variable" class="csl-entry">

Wang, Yanxin, Qibin Fan, and Li Zhu. 2018. “Variable Selection and Estimation Using a Continuous Approximation to the $L_0$ Penalty.” *Annals of the Institute of Statistical Mathematics* 70 (1): 191–214. <https://doi.org/10.1007/s10463-016-0588-3>.

</div>

<div id="ref-wang2016variable" class="csl-entry">

Wang, Yanxin, and Li Zhu. 2016. “Variable Selection and Parameter Estimation with the Atan Regularization Method.” *Journal of Probability and Statistics* 2016: 6495417. <https://doi.org/10.1155/2016/6495417>.

</div>

<div id="ref-vanwieringen2016ridge" class="csl-entry">

Wieringen, Wessel N. van, and Carel F. W. Peeters. 2016. “Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data.” *Computational Statistics & Data Analysis* 103: 284–303. <https://doi.org/10.1016/j.csda.2016.05.012>.

</div>

<div id="ref-zhang2010nearly" class="csl-entry">

Zhang, Cun-Hui. 2010. “Nearly Unbiased Variable Selection Under Minimax Concave Penalty.” *The Annals of Statistics* 38 (2): 894–942. <https://doi.org/10.1214/09-AOS729>.

</div>

<div id="ref-zou2006adaptive" class="csl-entry">

Zou, Hui. 2006. “The Adaptive Lasso and Its Oracle Properties.” *Journal of the American Statistical Association* 101 (476): 1418–29. <https://doi.org/10.1198/016214506000000735>.

</div>

<div id="ref-zou2005regularization" class="csl-entry">

Zou, Hui, and Trevor Hastie. 2005. “Regularization and Variable Selection via the Elastic Net.” *Journal of the Royal Statistical Society Series B: Statistical Methodology* 67 (2): 301–20. <https://doi.org/10.1111/j.1467-9868.2005.00527.x>.

</div>

</div>
