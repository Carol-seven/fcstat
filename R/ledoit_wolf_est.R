#' Ledoit-Wolf shrinkage estimator
#'
#' @description
#' Compute the Ledoit-Wolf shrinkage estimator for the covariance or correlation matrix.
#'
#' @param X A data matrix.
#'
#' @param method A string (default = "lin") specifying the method used in shrinkage,
#' includes: \enumerate{
#' \item Linear shrinkage: "lin" (Ledoit and Wolf, 2004)
#' \item Non-linear shrinkage (Ledoit and Wolf, 2015, 2017): \itemize{
#' \item "nlminb": non-linear shrinakge using the optimization routine "nlminb".
#' \item "nloptr": non-linear shrainkage using the optimization routine "nloptr".
#' }
#' }
#' See \code{\link[stats]{nlminb}} and \code{\link[nloptr]{nloptr}} for details.
#'
#' @param res A string (default = "cov") specifying the result matrix to be obtained,
#' either the covariance matrix ("cov") or the correlation matrix ("cor").
#'
#' @import nlshrink
#' @importFrom stats cov2cor
#'
#' @return A numeric matrix.
#'
#' @references \itemize{
#' \item Ledoit, Olivier and Wolf, Michael. (2004).
#' A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices.
#' \emph{Journal of Multivariate Analysis}, 88(2), 365--411.
#' \item Ledoit, Olivier and Wolf, Michael. (2015).
#' Spectrum Estimation: A Unified Framework for Covariance Matrix Estimation and PCA in
#' Large dimensions.
#' \emph{Journal of Multivariate Analysis}, 139, 360--384.
#' \item Ledoit, Olivier and Wolf, Michael. (2017).
#' Numerical Implementation of the QuEST Function.
#' \emph{Computational Statistics & Data Analysis}, 115, 199--223.
#' }
#'
#' @export

ledoit_wolf_est <- function(X, method = "lin", res = "cov") {
  if (method == "lin") {
    Sigma <- nlshrink::linshrink_cov(X)
  } else if (method == "nlminb") {
    Sigma <- nlshrink::nlshrink_cov(X, method = "nlminb")
  } else if (method == "nloptr") {
    Sigma <- nlshrink::nlshrink_cov(X, method = "nloptr")
  }
  if (res == "cov") {
    result <- Sigma
  } else if (res == "cor") {
    result <- cov2cor(Sigma)
  }
  return(result)
}

