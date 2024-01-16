#' Ledoit-Wolf shrinkage estimator
#'
#' @description
#' Compute the Ledoit-Wolf shrinkage estimator for the covariance or correlation matrix.
#'
#' @param X A data matrix.
#'
#' @param method A character string (default = "lin") specifying the method used in
#' shrinkage, includes: \enumerate{
#' \item Linear shrinkage: "lin" \insertCite{ledoit2004well}{fcstat}.
#' \item Non-linear shrinkage \insertCite{ledoit2015spectrum,ledoit2017numerical}{fcstat}:
#' \itemize{
#' \item "nlminb": non-linear shrinkage using the optimization routine "nlminb".
#' \item "nloptr": non-linear shrinkage using the optimization routine "nloptr".
#' }
#' }
#' See \code{\link[stats]{nlminb}} and \code{\link[nloptr]{nloptr}} for details.
#'
#' @param res A character string (default = "cov") specifying the result matrix to be
#' obtained, either the covariance matrix ("cov") or the correlation matrix ("cor").
#'
#' @import nlshrink
#' @importFrom stats cov2cor
#' @importFrom Rdpack reprompt
#'
#' @return A numeric matrix.
#'
#' @references
#' \insertAllCited{}
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

