#' Ledoit-Wolf shrinkage estimator
#'
#' @description
#' Compute the Ledoit-Wolf shrinkage estimator for the covariance/correlation
#' matrix.
#'
#' @param X A data matrix.
#'
#' @param method A character string (default = "linshrink") specifying
#' the method used in shrinkage, includes: \enumerate{
#' \item "linshrink": Linear shrinkage \insertCite{ledoit2004well}{spice}.
#' \item "nlshrink": Non-linear shrinkage
#' \insertCite{ledoit2015spectrum,ledoit2017numerical}{spice}.
#' }
#' See \code{\link[nlshrink]{linshrink_cov}} and
#' \code{\link[nlshrink]{nlshrink_cov}} for details.
#'
#' @param res A character string (default = "cov") specifying the result matrix
#' to be obtained, either the covariance matrix ("cov") or the correlation
#' matrix ("cor").
#'
#' @return
#' A numeric matrix.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom nlshrink linshrink_cov nlshrink_cov
#' @importFrom stats cov2cor
#' @importFrom Rdpack reprompt
#'
#' @export

ledoit_wolf_est <- function(X, method = "linshrink", res = "cov") {
  if (method == "linshrink") {
    Sigma <- nlshrink::linshrink_cov(X)
  } else if (method == "nlshrink") {
    Sigma <- nlshrink::nlshrink_cov(X, method = "nlminb")
  }
  if (res == "cov") {
    result <- Sigma
  } else if (res == "cor") {
    result <- cov2cor(Sigma)
  }
  return(result)
}

