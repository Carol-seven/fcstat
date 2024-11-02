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
#' \item "ADMMsigma": the function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "CovTools": the function from \code{\link[CovTools]{PreEst.glasso}}.
#' \item "CVglasso": the function from \code{\link[CVglasso]{CVglasso}}.
#' \item "Glarmadillo": the function from \code{\link[Glarmadillo]{glarma}}.
#' \item "glasso": the function from \code{\link[glasso]{glasso}}.
#' \item "GLassoElnetFast": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "glassoFast": the function from \code{\link[glassoFast]{glassoFast}}.
#' \item "huge": the function from \code{\link[huge]{huge.glasso}}.
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
#' @return A list with the following: \describe{
#' \item{hatOmega}{The estimated precision matrix.}
#' \item{niter}{Number of iterations.}
#' }
#'
#' @noRd

fcstat_method <- function(method, X = NULL, S = NULL,
                          lambda = NULL, gamma = NULL,
                          pkgopt = NULL) {
  if (method == "glasso") {
    if (pkgopt == "ADMMsigma") {
      res <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = 1, diagonal = TRUE)
      result <- list(hatOmega = res$Z, niter = res$Iterations)
    } else if (pkgopt == "CovTools") {
      res <- CovTools::PreEst.glasso(X = X, method = list(type = "fixed", param = lambda))
      result <- list(hatOmega = res$C, niter = NULL)
    } else if (pkgopt == "CVglasso") {
      res <- CVglasso::CVglasso(S = S, lam = lambda, diagonal = TRUE)
      result <- list(hatOmega = res$Omega, niter = res$Iterations)
    } else if (pkgopt == "Glarmadillo") {
      res <- Glarmadillo::glarma(s = S, rho = lambda)
      result <- list(hatOmega = res$Theta, niter = res$iter)
    } else if (pkgopt == "glasso") {
      res <- glasso::glasso(s = S, rho = lambda, penalize.diagonal = TRUE, start = "cold")
      result <- list(hatOmega = res$wi, niter = res$niter)
    } else if (pkgopt == "GLassoElnetFast") {
      res <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 1, penalize.diagonal = TRUE)
      result <- list(hatOmega = res$Theta, niter = res$niter)
    } else if (pkgopt == "glassoFast") {
      res <- glassoFast::glassoFast(S = S, rho = lambda, start = "cold")
      result <- list(hatOmega = res$wi, niter = res$niter)
    } else if (pkgopt == "huge") {
      res <- huge::huge.glasso(x = S, lambda = lambda, verbose = FALSE)
      result <- list(hatOmega = res$icov[[1]], niter = NULL)
    }
  } else if (method == "ridge") {
    if (pkgopt == "ADMMsigma") {
      res <- ADMMsigma::RIDGEsigma(S = S, lam = lambda)
      result <- list(hatOmega = res$Omega, niter = NULL)
    } else if (pkgopt == "GLassoElnetFast") {
      res <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 0)
      result <- list(hatOmega = res$Theta, niter = res$niter)
    } else if (pkgopt == "porridge") {
      res <- porridge::ridgePgen(S = S, lambda = matrix(lambda, ncol(S), ncol(S)), target = matrix(0, ncol(S), ncol(S)))
      result <- list(hatOmega = res, niter = NULL)
    } else if (pkgopt == "rags2ridges") {
      res <- rags2ridges::ridgeP(S = S, lambda = lambda, target = matrix(0, ncol(S), ncol(S)))
      result <- list(hatOmega = res, niter = NULL)
    }
  } else if (method == "elnet") {
    if (pkgopt == "ADMMsigma") {
      res <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = gamma)
      result <- list(hatOmega = res$Z, niter = res$Iterations)
    } else if (pkgopt == "GLassoElnetFast") {
      res <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = gamma)
      result <- list(hatOmega = res$Theta, niter = res$niter)
    }
  } else if (method == "clime") {
    if (pkgopt == "clime") {
      res <- clime::clime(x = S, lambda = lambda, sigma = TRUE, standardize = FALSE, linsolver = "simplex")
      result <- list(hatOmega = res$Omegalist[[1]], niter = NULL)
    } else if (pkgopt == "flare") {
      res <- flare::sugm(data = S, lambda = lambda, method = "clime", verbose = FALSE)
      result <- list(hatOmega = res$icov[[1]], niter = res$ite)
    }
  } else if (method == "tiger") {
    if (pkgopt == "flare") {
      res <- flare::sugm(data = X, lambda = lambda, method = "tiger", verbose = FALSE)
      result <- list(hatOmega = res$icov[[1]], niter = res$ite)
    } else if (pkgopt == "huge") {
      res <- huge::huge.tiger(x = X, lambda = lambda, verbose = FALSE)
      result <- list(hatOmega = res$icov[[1]], niter = NULL)
    }
  }
  return(result)
}

