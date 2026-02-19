#' Estimation for Precision Matrix
#'
#' @description
#' Estimate a sparse precision matrix using a collection of statistical methods.
#'
#' @param method A character string specifying the statistical method for
#' estimating precision matrix. Available options include: \enumerate{
#' \item "glasso": Graphical lasso \insertCite{friedman2008sparse}{spice}.
#' \item "ridge": Graphical ridge \insertCite{vanwieringen2016ridge}{spice}.
#' \item "elnet": Graphical elastic net \insertCite{zou2005regularization}{spice}.
#' \item "clime": Constrained L1-minimization for inverse (covariance) matrix
#' estimation \insertCite{cai2011aconstrained}{spice}.
#' \item "tiger": Tuning-insensitive graph estimation and regression
#' \insertCite{liu2017tiger}{spice}.
#' }
#'
#' @param X An \eqn{n \times p} data matrix with sample size \eqn{n} and
#' dimension \eqn{p}.
#'
#' @param S A \eqn{p \times p} sample covariance/correlation matrix with
#' dimension \eqn{p}.
#'
#' @param lambda A non-negative numeric value specifying the regularization
#' parameter.
#'
#' @param gamma A numeric value specifying the additional parameter for
#' \code{method} set to \code{"elnet"}.
#'
#' @param pkg A character string specifying the package option to use.
#' The available options depend on the selected method: \enumerate{
#' \item For \code{method = "glasso"}: \itemize{
#' \item "ADMMsigma": The function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "CovTools": The function from \code{\link[CovTools]{PreEst.glasso}}.
#' \item "CVglasso": The function from \code{\link[CVglasso]{CVglasso}}.
#' \item "Glarmadillo": The function from \code{\link[Glarmadillo]{glarma}}.
#' \item "glasso": The function from \code{\link[glasso]{glasso}}.
#' \item "GLassoElnetFast": The function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "glassoFast": The function from \code{\link[glassoFast]{glassoFast}}.
#' \item "huge": The function from \code{\link[huge]{huge.glasso}}.
#' }
#' \item For \code{method = "ridge"}: \itemize{
#' \item "ADMMsigma": The function from \code{\link[ADMMsigma]{RIDGEsigma}}.
#' \item "GLassoElnetFast": The function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "porridge": The function from \code{\link[porridge]{ridgePgen}}.
#' \item "rags2ridges": The function from \code{\link[rags2ridges]{ridgeP}}.
#' }
#' \item For \code{method = "elnet"}: \itemize{
#' \item "ADMMsigma": The function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "GLassoElnetFast": The function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' }
#' \item For \code{method = "clime"}: \itemize{
#' \item "clime": The function from \code{\link[clime]{clime}}.
#' \item "flare": The function from \code{\link[flare]{sugm}}.
#' }
#' \item For \code{method = "tiger"}: \itemize{
#' \item "flare": The function from \code{\link[flare]{sugm}}.
#' \item "huge": The function from \code{\link[huge]{huge.tiger}}.
#' }
#' }
#'
#' @import RBGL
#' @import graph
#' @importFrom ADMMsigma ADMMsigma
#' @importFrom CovTools PreEst.glasso
#' @importFrom CVglasso CVglasso
#' @importFrom Glarmadillo glarma
#' @importFrom glasso glasso
#' @importFrom GLassoElnetFast gelnet
#' @importFrom glassoFast glassoFast
#' @importFrom huge huge.glasso
#' @importFrom ADMMsigma RIDGEsigma
#' @importFrom porridge ridgePgen
#' @importFrom rags2ridges ridgeP
#' @importFrom clime clime
#' @importFrom flare sugm
#' @importFrom huge huge.tiger
#' @importFrom Rdpack reprompt
#'
#' @return Estimated precision matrix.
#'
#' @noRd

spice_method <- function(method, X = NULL, S = NULL,
                         lambda = NULL, gamma = NULL,
                         pkg = NULL) {
  if (method == "glasso") {
    if (pkg == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = 1, diagonal = TRUE)$Z
    } else if (pkg == "CovTools") {
      hatOmega <- CovTools::PreEst.glasso(X = X, method = list(type = "fixed", param = lambda))$C
    } else if (pkg == "CVglasso") {
      hatOmega <- CVglasso::CVglasso(S = S, lam = lambda, diagonal = TRUE)$Omega
    } else if (pkg == "Glarmadillo") {
      hatOmega <- Glarmadillo::glarma(s = S, rho = lambda)$Theta
    } else if (pkg == "glasso") {
      hatOmega <- glasso::glasso(s = S, rho = lambda, penalize.diagonal = TRUE, start = "cold")$wi
    } else if (pkg == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 1, penalize.diagonal = TRUE)$Theta
    } else if (pkg == "glassoFast") {
      hatOmega <- glassoFast::glassoFast(S = S, rho = lambda, start = "cold")$wi
    } else if (pkg == "huge") {
      hatOmega <- huge::huge.glasso(x = S, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  } else if (method == "ridge") {
    if (pkg == "ADMMsigma") {
      hatOmega <- ADMMsigma::RIDGEsigma(S = S, lam = lambda)$Omega
    } else if (pkg == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 0)$Theta
    } else if (pkg == "porridge") {
      hatOmega <- porridge::ridgePgen(S = S, lambda = matrix(lambda, ncol(S), ncol(S)), target = matrix(0, ncol(S), ncol(S)))
    } else if (pkg == "rags2ridges") {
      hatOmega <- rags2ridges::ridgeP(S = S, lambda = lambda, target = matrix(0, ncol(S), ncol(S)))
    }
  } else if (method == "elnet") {
    if (pkg == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = gamma)$Z
    } else if (pkg == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = gamma)$Theta
    }
  } else if (method == "clime") {
    if (pkg == "clime") {
      hatOmega <- clime::clime(x = S, lambda = lambda, sigma = TRUE, standardize = FALSE, linsolver = "simplex")$Omegalist[[1]]
    } else if (pkg == "flare") {
      hatOmega <- flare::sugm(data = S, lambda = lambda, method = "clime", verbose = FALSE)$icov[[1]]
    }
  } else if (method == "tiger") {
    if (pkg == "flare") {
      hatOmega <- flare::sugm(data = X, lambda = lambda, method = "tiger", verbose = FALSE)$icov[[1]]
    } else if (pkg == "huge") {
      hatOmega <- huge::huge.tiger(x = X, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  }
  return(hatOmega)
}

