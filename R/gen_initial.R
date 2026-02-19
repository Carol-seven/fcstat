#' Initial
#'
#' @description
#' The initial estimate for \code{method} set to \code{"atan"}, \code{"exp"},
#' \code{"scad"}, and \code{"mcp"}; or specifies \eqn{\tilde{\Omega}} of
#' the adaptive weight for \code{method = "adapt"}, calculated as
#' \eqn{|\tilde{\omega}_{ij}|^{-\gamma}},
#' where \eqn{\tilde{\Omega} := (\tilde{\omega}_{ij})}.
#'
#' @param X An \eqn{n \times p} data matrix with sample size \eqn{n} and
#' dimension \eqn{p}.
#'
#' @param S A \eqn{p \times p} sample covariance/correlation matrix with
#' dimension \eqn{p}.
#'
#' @param base A character string specifying the calculation base,
#' either the covariance matrix ("cov") or the correlation matrix ("cor").
#'
#' @param initial A \eqn{p \times p} matrix or
#' a \eqn{p \times p \times \mathrm{npara}}
#' (the number of all combinations of \code{lambda} and \code{gamma}) array
#' specifying the initial estimate for \code{method} set to \code{"atan"},
#' \code{"exp"}, \code{"scad"}, and \code{"mcp"};
#' or specifying \eqn{\tilde{\Omega}} of the adaptive weight for
#' \code{method = "adapt"}, calculated as
#' \eqn{\lvert\tilde{\omega}_{ij}\rvert^{-\gamma}},
#' where \eqn{\tilde{\Omega} := (\tilde{\omega}_{ij})}.
#' Some options are also offered when a character string is provided
#' (default = "glasso"), including:
#' \itemize{
#' \item "glasso": Use the precision matrix estimate derived from the graphical
#' lasso.
#' \item "invS": Use the inverse calculation base matrix if the matrix is
#' invertible.
#' \item "linshrink": Use the precision matrix estimate derived from
#' Ledoit-Wolf linear shrinkage estimator of the population covariance matrix
#' \insertCite{ledoit2004well}{spice}.
#' \item "nlshrink": Use the precision matrix estimate derived from
#' Ledoit-Wolf non-linear shrinkage estimator of the population covariance
#' matrix \insertCite{ledoit2015spectrum,ledoit2017numerical}{spice}.
#' }
#'
#' @param pkg A character string specifying the package option to use.
#' This argument only works when \code{initial = "glasso"}. \itemize{
#' \item "glasso": The function from \code{\link[glasso]{glasso}}.
#' \item "GLassoElnetFast": The function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "glassoFast": The function from \code{\link[glassoFast]{glassoFast}}.
#' }
#'
#' @param lambda A non-negative numeric vector specifying the grid for
#' the regularization parameter.
#'
#' @importFrom glassoFast glassoFast
#' @importFrom methods is
#'
#' @return
#' The initial estimate or \eqn{\tilde{\Omega}} of the adaptive weight.
#'
#' @noRd

gen_initial <- function(X, S, base, initial, lambda, pkg) {

  if (is(initial, "matrix")) {
    Omega <- replicate(length(lambda), initial, simplify = FALSE)

  } else if (is(initial, "array")) {
    Omega <- lapply(1:dim(initial)[3], function(z) initial[,,z])

  } else if (initial == "invS") {
    Omega <- replicate(length(lambda), solve(S), simplify = FALSE)

  } else if (initial == "linshrink") {
    Omega <- replicate(length(lambda),
                       solve(ledoit_wolf_est(X, method = "linshrink", res = base)),
                       simplify = FALSE)

  } else if (initial == "nlshrink") {
    Omega <- replicate(length(lambda),
                       solve(ledoit_wolf_est(X, method = "nlshrink", res = base)),
                       simplify = FALSE)

  } else if (initial == "glasso") {
    if (pkg == "glasso") {
      Omega <- lapply(lambda, function(z) {
        glasso::glasso(s = S, rho = z, penalize.diagonal = TRUE, start = "cold")$wi
      })
    } else if (pkg == "GLassoElnetFast") {
      Omega <- lapply(lambda, function(z) {
        GLassoElnetFast::gelnet(S = S, lambda = z, alpha = 1, penalize.diagonal = TRUE)$Theta
      })
    } else if (pkg == "glassoFast") {
      Omega <- lapply(lambda, function(z) {
        glassoFast::glassoFast(S = S, rho = z, start = "cold")$wi
      })
    }
  }

  return(Omega)
}

