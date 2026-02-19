# Derivative of a Non-Convex Penalty

Compute the derivative of a specified non-convex regularization penalty.

## Usage

``` r
ncv_deriv(penalty, Omega, lambda, gamma)
```

## Arguments

- penalty:

  A character string specifying the non-convex penalty to use. Available
  options include:

  1.  "adapt": Adaptive lasso (Zou 2006; Fan et al. 2009) .

  2.  "atan": Arctangent type penalty (Wang and Zhu 2016) .

  3.  "exp": Exponential type penalty (Wang et al. 2018) .

  4.  "mcp": Minimax concave penalty (Zou 2006) .

  5.  "scad": Smoothly clipped absolute deviation (Fan and Li 2001; Fan
      et al. 2009) .

- Omega:

  The precision matrix.

- lambda:

  A non-negative numeric value specifying the regularization parameter.

- gamma:

  A numeric value specifying the additional parameter for the penalty
  function. The defaults are:

  1.  "adapt": 0.5

  2.  "atan": 0.005

  3.  "exp": 0.01

  4.  "mcp": 3

  5.  "scad": 3.7

## Value

A numeric matrix.

## References

Fan J, Feng Y, Wu Y (2009). “Network Exploration via the Adaptive LASSO
and SCAD Penalties.” *The Annals of Applied Statistics*, **3**(2),
521–541. [doi:10.1214/08-aoas215](https://doi.org/10.1214/08-aoas215)
.  
  
Fan J, Li R (2001). “Variable Selection via Nonconcave Penalized
Likelihood and its Oracle Properties.” *Journal of the American
Statistical Association*, **96**(456), 1348–1360.
[doi:10.1198/016214501753382273](https://doi.org/10.1198/016214501753382273)
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
