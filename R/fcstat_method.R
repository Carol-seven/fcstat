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
#' @param gamma A scalar specifying the hyperparameter for \code{method} set to
#' \code{"elnet"}.
#'
#' @param pkgopt A character string specifying the package option to use. The available
#' options depend on the selected method: \enumerate{
#' \item For \code{method = "glasso"}: \itemize{
#' \item "ADMMsigma_diagtrue": the function from \code{\link[ADMMsigma]{ADMMsigma}} with
#' warm-starts, and the diagonal elements of the estimated precision matrix are penalized.
#' \item "ADMMsigma_diagfalse": the function from \code{\link[ADMMsigma]{ADMMsigma}} with
#' warm-starts, and the diagonal elements of the estimated precision matrix are not
#' penalized.
#' \item "CovTools": the function from \code{\link[CovTools]{PreEst.glasso}}.
#' \item "CVglasso_diagtrue": the function from \code{\link[CVglasso]{CVglasso}} with
#' warm-starts, and the diagonal elements of the estimated precision matrix are penalized.
#' \item "CVglasso_diagfalse": the function from \code{\link[CVglasso]{CVglasso}} with
#' warm-starts, and the diagonal elements of the estimated precision matrix are not
#' penalized.
#' \item "glasso_cold_diagtrue": the function from \code{\link[glasso]{glasso}}, and the
#' diagonal elements of the estimated precision matrix are penalized.
#' \item "glasso_cold_diagfalse": the function from \code{\link[glasso]{glasso}}, and the
#' diagonal elements of the estimated precision matrix are not penalized.
#' \item "glasso_warm_diagtrue": the function from \code{\link[glasso]{glasso}} with
#' warm-starts, and the diagonal elements of the estimated precision matrix are penalized.
#' \item "glasso_warm_diagfalse": the function from \code{\link[glasso]{glasso}} with
#' warm-starts, and the diagonal elements of the estimated precision matrix are not
#' penalized.
#' \item "GLassoElnetFast_cold_diagtrue": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}, and the diagonal
#' elements of the estimated precision matrix are penalized.
#' \item "GLassoElnetFast_cold_diagfalse": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}, and the diagonal
#' elements of the estimated precision matrix are not penalized.
#' \item "GLassoElnetFast_warm_diagtrue": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet} with warm-starts, and
#' the diagonal elements of the estimated precision matrix are penalized.
#' \item "GLassoElnetFast_warm_diagfalse": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet} with warm-starts, and
#' the diagonal elements of the estimated precision matrix are not penalized.
#' \item "glassoFast_cold": the function from \code{\link[glassoFast]{glassoFast}}, and
#' the diagonal elements of the estimated precision matrix are penalized.
#' \item "glassoFast_warm": the function from \code{\link[glassoFast]{glassoFast}} with
#' warm-starts, and the diagonal elements of the estimated precision matrix are penalized.
#' \item "huge": the function from \code{\link[huge]{huge.glasso}}, and the diagonal
#' elements of the estimated precision matrix are penalized.
#' }
#' \item For \code{method = "ridge"}: \itemize{
#' \item "ADMMsigma": the function from \code{\link[ADMMsigma]{RIDGEsigma}}.
#' \item "GLassoElnetFast": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "porridge": the function from \code{\link[porridge]{ridgePgen}}.
#' \item "rags2ridges": the function from \code{\link[rags2ridges]{ridgeP}}.
#' }
#' \item For \code{method = "elnet"}: \itemize{
#' \item "ADMMsigma": the function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "GLassoElnetFast": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' }
#' \item For \code{method = "clime"}: \itemize{
#' \item "clime": the function from \code{\link[clime]{clime}}.
#' \item "flare": the function from \code{\link[flare]{sugm}}.
#' }
#' \item For \code{method = "tiger"}: \itemize{
#' \item "flare": the function from \code{\link[flare]{sugm}}.
#' \item "huge": the function from \code{\link[huge]{huge.tiger}}.
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

fcstat_method <- function(method, X = NULL, S = NULL,
                          lambda = NULL, gamma = NULL,
                          pkgopt = NULL) {
  if (method == "glasso") {
    if (pkgopt == "ADMMsigma_diagtrue") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = 1, diagonal = TRUE)$Z
    } else if (pkgopt == "ADMMsigma_diagfalse") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = 1, diagonal = FALSE)$Z
    } else if (pkgopt == "CovTools") {
      hatOmega <- CovTools::PreEst.glasso(X = X, method = list(type = "fixed", param = lambda))$C
    } else if (pkgopt == "CVglasso_diagtrue") {
      hatOmega <- CVglasso::CVglasso(S = S, lam = lambda, diagonal = TRUE)$Omega
    } else if (pkgopt == "CVglasso_diagfalse") {
      hatOmega <- CVglasso::CVglasso(S = S, lam = lambda, diagonal = FALSE)$Omega
    } else if (pkgopt == "glasso_cold_diagtrue") {
      hatOmega <- glasso::glasso(s = S, rho = lambda, penalize.diagonal = TRUE, start = "cold")$wi
    } else if (pkgopt == "glasso_cold_diagfalse") {
      hatOmega <- glasso::glasso(s = S, rho = lambda, penalize.diagonal = FALSE, start = "cold")$wi
    } else if (pkgopt == "glasso_warm_diagtrue") {
      hatOmega <- glasso::glasso(s = S, rho = lambda, penalize.diagonal = TRUE, start = "warm",
                                 w.init = S + lambda, wi.init = diag(ncol(S)))$wi
    } else if (pkgopt == "glasso_warm_diagfalse") {
      hatOmega <- glasso::glasso(s = S, rho = lambda, penalize.diagonal = FALSE, start = "warm",
                                 w.init = S + lambda, wi.init = diag(ncol(S)))$wi
    } else if (pkgopt == "GLassoElnetFast_cold_diagtrue") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 1, penalize.diagonal = TRUE)$Theta
    } else if (pkgopt == "GLassoElnetFast_cold_diagfalse") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 1, penalize.diagonal = FALSE)$Theta
    } else if (pkgopt == "GLassoElnetFast_warm_diagtrue") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 1, penalize.diagonal = TRUE,
                                          W = S + lambda, Theta = diag(ncol(S)))$Theta
    } else if (pkgopt == "GLassoElnetFast_warm_diagfalse") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 1, penalize.diagonal = FALSE,
                                          W = S + lambda, Theta = diag(ncol(S)))$Theta
    } else if (pkgopt == "glassoFast_cold") {
      hatOmega <- glassoFast::glassoFast(S = S, rho = lambda, start = "cold")$wi
    } else if (pkgopt == "glassoFast_warm") {
      hatOmega <- glassoFast::glassoFast(S = S, rho = lambda, start = "warm",
                                         w.init = S + lambda, wi.init = diag(ncol(S)))$wi
    } else if (pkgopt == "huge") {
      hatOmega <- huge::huge.glasso(x = S, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  } else if (method == "ridge") {
    if (pkgopt == "ADMMsigma") {
      hatOmega <- ADMMsigma::RIDGEsigma(S = S, lam = lambda)$Omega
    } else if (pkgopt == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 0)$Theta
    } else if (pkgopt == "porridge") {
      hatOmega <- porridge::ridgePgen(S = S, lambda = matrix(lambda, ncol(S), ncol(S)), target = matrix(0, ncol(S), ncol(S)))
    } else if (pkgopt == "rags2ridges") {
      hatOmega <- rags2ridges::ridgeP(S = S, lambda = lambda, target = matrix(0, ncol(S), ncol(S)))
    }
  } else if (method == "elnet") {
    if (pkgopt == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = gamma)$Z
    } else if (pkgopt == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = gamma)$Theta
    }
  } else if (method == "clime") {
    if (pkgopt == "clime") {
      hatOmega <- clime::clime(x = S, lambda = lambda, sigma = TRUE, standardize = FALSE, linsolver = "simplex")$Omegalist[[1]]
    } else if (pkgopt == "flare") {
      hatOmega <- flare::sugm(data = S, lambda = lambda, method = "clime", verbose = FALSE)$icov[[1]]
    }
  } else if (method == "tiger") {
    if (pkgopt == "flare") {
      hatOmega <- flare::sugm(data = X, lambda = lambda, method = "tiger", verbose = FALSE)$icov[[1]]
    } else if (pkgopt == "huge") {
      hatOmega <- huge::huge.tiger(x = X, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  }
  return(hatOmega)
}

