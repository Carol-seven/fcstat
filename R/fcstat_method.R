#' Estimation for precision matrix
#'
#' @description
#' Estimate a sparse precision matrix using a collection of statistical methods.
#'
#' @param method A character string specifying the statistical method for estimating
#' precision matrix. Available options include: \enumerate{
#' \item "glasso": graphical lasso \insertCite{friedman2008sparse}{fcstat}.
#' \item "ridge": graphical ridge \insertCite{vanwieringen2016ridge}{fcstat}.
#' \item "elnet": graphical elastic net \insertCite{zou2005regularization}{fcstat}.
#' \item "clime": constrained L1-minimization for inverse (covariance) matrix estimation
#' \insertCite{cai2011aconstrained}{fcstat}.
#' \item "tiger": tuning-insensitive graph estimation and regression
#' \insertCite{liu2017tiger}{fcstat}.
#' }
#'
#' @param X An n-by-p data matrix with sample size n and dimension p.
#'
#' @param S A p-by-p sample covariance/correlation matrix with dimension p.
#'
#' @param lambda A non-negative scalar specifying the regularization parameter.
#'
#' @param gamma Grid of scalars specifying the hyperparameter for the chosen \code{method}.
#' Default values: \enumerate{
#' \item "elnet": a sequence from 0.1 to 0.9 with increments of 0.1
#' \item "adapt": 0.5
#' \item "atan": 0.005
#' \item "exp": 0.01
#' \item "scad": 3.7
#' \item "mcp": 3
#' }
#'
#' @param gamma A scalar specifying the hyperparameter for \code{method} set to
#' \code{"elnet"}.
#'
#' @param utilopt A character string specifying the utility option to use. The available
#' options depend on the chosen method: \enumerate{
#' \item For \code{method = "glasso"}: \itemize{
#' \item "ADMMsigma": the utility function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "CovTools": the utility function from \code{\link[CovTools]{PreEst.glasso}}.
#' \item "CVglasso": the utility function from \code{\link[CVglasso]{CVglasso}}.
#' \item "glasso": the utility function from \code{\link[glasso]{glasso}}.
#' \item "GLassoElnetFast": the utility function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "glassoFast": the utility function from \code{\link[glassoFast]{glassoFast}}.
#' \item "huge": the utility function from \code{\link[huge]{huge.glasso}}.
#' }
#' \item For \code{method = "ridge"}: \itemize{
#' \item "ADMMsigma": the utility function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "GLassoElnetFast": the utility function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "porridge": the utility function from \code{\link[porridge]{ridgePgen}}.
#' \item "rags2ridges": the utility function from \code{\link[rags2ridges]{ridgeP}}.
#' }
#' \item For \code{method = "elnet"}: \itemize{
#' \item "ADMMsigma": the utility function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "GLassoElnetFast": the utility function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' }
#' \item For \code{method = "clime"}: \itemize{
#' \item "clime": the utility function from \code{\link[clime]{clime}}.
#' \item "flare": the utility function from \code{\link[flare]{sugm}}.
#' }
#' \item For \code{method = "tiger"}: \itemize{
#' \item "flare": the utility function from \code{\link[flare]{sugm}}.
#' \item "huge": the utility function from \code{\link[huge]{huge.tiger}}.
#' }
#' }
#'
#' @import RBGL
#' @import graph
#' @importFrom ADMMsigma ADMMsigma
#' @importFrom CovTools PreEst.glasso
#' @importFrom CVglasso CVglasso
#' @importFrom glasso glasso
#' @importFrom GLassoElnetFast gelnet
#' @importFrom glassoFast glassoFast
#' @importFrom huge huge.glasso
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

fcstat_method <- function(method, X = NULL, S = NULL,
                          lambda = NULL, gamma = NULL,
                          utilopt = NULL) {
  if (method == "glasso") {
    if (utilopt == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = 1)$Omega
    } else if (utilopt == "CovTools") {
      hatOmega <- CovTools::PreEst.glasso(X = X, method = list(type = "fixed", param = lambda))$C
    } else if (utilopt == "CVglasso") {
      hatOmega <- CVglasso::CVglasso(S = S, lam = lambda)$Omega
    } else if (utilopt == "glasso") {
      hatOmega <- glasso::glasso(s = S, rho = lambda)$wi
    } else if (utilopt == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 1)$Theta
    } else if (utilopt == "glassoFast") {
      hatOmega <- glassoFast::glassoFast(S = S, rho = lambda)$wi
    } else if (utilopt == "huge") {
      hatOmega <- huge::huge.glasso(x = S, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  } else if (method == "ridge") {
    if (utilopt == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = 0)$Omega
    } else if (utilopt == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 0)$Theta
    } else if (utilopt == "porridge") {
      hatOmega <- porridge::ridgePgen(S = S, lambda = matrix(lambda, ncol(S), ncol(S)), target = matrix(0, ncol(S), ncol(S)))
    } else if (utilopt == "rags2ridges") {
      hatOmega <- rags2ridges::ridgeP(S = S, lambda = lambda, target = matrix(0, ncol(S), ncol(S)))
    }
  } else if (method == "elnet") {
    if (utilopt == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = gamma)$Omega
    } else if (utilopt == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = gamma)$Theta
    }
  } else if (method == "clime") {
    if (utilopt == "clime") {
      hatOmega <- clime::clime(x = S, lambda = lambda, sigma = TRUE, standardize = FALSE, linsolver = "simplex")$Omegalist[[1]]
    } else if (utilopt == "flare") {
      hatOmega <- flare::sugm(data = S, lambda = lambda, method = "clime", verbose = FALSE)$icov[[1]]
    }
  } else if (method == "tiger") {
    if (utilopt == "flare") {
      hatOmega <- flare::sugm(data = X, lambda = lambda, method = "tiger", verbose = FALSE)$icov[[1]]
    } else if (utilopt == "huge") {
      hatOmega <- huge::huge.tiger(x = X, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  }
  return(hatOmega)
}

