#' Initial
#'
#' @description
#' The initial estimate for \code{method} set to \code{"atan"}, \code{"exp"},
#' \code{"scad"}, or \code{"mcp"}; or specifies \eqn{\tilde{\Omega}} of the adaptive
#' weight for \code{method = "adapt"}, calculated as \eqn{|\tilde{\omega}_{ij}|^{-\gamma}},
#' where \eqn{\tilde{\Omega} := (\tilde{\omega}_{ij})}.
#'
#' @param X An n-by-p data matrix with sample size n and dimension p.
#'
#' @param S A p-by-p sample covariance/correlation matrix with dimension p.
#'
#' @param base A character string specifying the calculation base, either the covariance
#' matrix ("cov") or the correlation matrix ("cor").
#'
#' @param initial A p-by-p matrix or a p-by-p-by-npara (the number of all combinations of
#' \code{lambda} and \code{gamma}) array specifying the initial estimate for \code{method}
#' or \eqn{\tilde{\Omega}} of the adaptive weight. Some options are also offered when a
#' character string is provided (default "linshrink"), including:
#' \itemize{
#' \item "invS": use the inverse calculation base matrix.
#' \item "linshrink": use the precision matrix estimate derived from Ledoit-Wolf linear
#' shrinakge estimator of the population covariance matrix
#' \insertCite{ledoit2004well}{fcstat}.
#' \item "nlshrink": use the precision matrix estimate derived from Ledoit-Wolf non-linear
#' shrinakge estimator of the population covariance matrix
#' \insertCite{ledoit2015spectrum,ledoit2017numerical}{fcstat}.
#' \item "glasso": use the precision matrix estimate derived from the graphical lasso.
#' }
#'
#' @param parameter The parameter grid.
#'
#' @importFrom glassoFast glassoFast
#' @importFrom methods is
#'
#' @return The initial estimate or \eqn{\tilde{\Omega}} of the adaptive weight.
#'
#' @noRd

gen_initial <- function(X, S, base, initial, parameter) {

  if (is(initial, "matrix")) {
    Omega <- replicate(nrow(parameter), initial, simplify = FALSE)

  } else if (is(initial, "array")) {
    Omega <- lapply(1:dim(initial)[3], function(z) initial[,,z])

  } else if (initial == "invS") {
    Omega <- replicate(nrow(parameter), solve(S), simplify = FALSE)

  } else if (initial == "linshrink") {
    Omega <- replicate(nrow(parameter),
                       solve(ledoit_wolf_est(X, method = "linshrink", res = base)),
                       simplify = FALSE)

  } else if (initial == "nlshrink") {
    Omega <- replicate(nrow(parameter),
                       solve(ledoit_wolf_est(X, method = "nlshrink", res = base)),
                       simplify = FALSE)

  } else if (initial == "glasso") {
    Omega <- lapply(parameter$lambda, function(z) glassoFast::glassoFast(S, rho = z)$wi)
  }

  return(Omega)
}
