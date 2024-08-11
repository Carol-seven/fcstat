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
#' @param target A p-by-p symmetric matrix or a scalar (default = 0) serves as the value
#' for all elements, specifying the target matrix for \code{method = "ridge"} and
#' \code{method = "elnet"}.
#'
#' @param utilopt A character string specifying the utility option to use. The available
#' options depend on the chosen method: \enumerate{
#' \item For \code{method = "glasso"}: \itemize{
#' \item "CVglasso": the utility function from \code{\link[CVglasso]{CVglasso}}.
#' \item "CovTools": the utility function from \code{\link[CovTools]{PreEst.glasso}}.
#' \item "glasso": the utility function from \code{\link[glasso]{glasso}}.
#' \item "glassoFast": the utility function from \code{\link[glassoFast]{glassoFast}}.
#' \item "huge": the utility function from \code{\link[huge]{huge.glasso}}.
#' }
#' \item For \code{method = "ridge"}: \itemize{
#' \item "porridge": the utility function from \code{\link[porridge]{ridgePgen}}.
#' \item "rags2ridges": the utility function from \code{\link[rags2ridges]{ridgeP}}.
#' }
#' \item For \code{method = "elnet"}: \itemize{
#' \item "ADMMsigma": the utility function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "GLassoElnetFast": the utility function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' }
#' \item For \code{method = "clime"}: \itemize{
#' \item "clime_primaldual": the utility function from \code{\link[clime]{clime}}
#' with the linsolver \code{primaldual}.
#' \item "clime_simplex": the utility function from \code{\link[clime]{clime}}
#' with the linsolver \code{simplex}.
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
#' @importFrom CVglasso CVglasso
#' @importFrom CovTools PreEst.glasso
#' @importFrom glasso glasso
#' @importFrom glassoFast glassoFast
#' @importFrom huge huge.glasso
#' @importFrom porridge ridgePgen
#' @importFrom rags2ridges ridgeP
#' @importFrom ADMMsigma ADMMsigma
#' @importFrom GLassoElnetFast gelnet
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
                          target = 0, utilopt = NULL) {
  if (method == "glasso") {
    if (utilopt == "CVglasso") {
      hatOmega <- CVglasso::CVglasso(S = S, lam = lambda)$Omega
    } else if (utilopt == "CovTools") {
      hatOmega <- CovTools::PreEst.glasso(X, method = list(type = "fixed", param = lambda))$C
    } else if (utilopt == "glasso") {
      hatOmega <- glasso::glasso(s = S, rho = lambda)$wi
    } else if (utilopt == "glassoFast") {
      hatOmega <- glassoFast::glassoFast(S = S, rho = lambda)$wi
    } else if (utilopt == "huge") {
      hatOmega <- huge::huge.glasso(x = S, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  } else if (method == "ridge") {
    if (is.numeric(target) & length(target) == 1) {
      target <- matrix(target, ncol(S), ncol(S))
    }
    if (utilopt == "porridge") {
      hatOmega <- porridge::ridgePgen(S = S, lambda = matrix(lambda, ncol(S), ncol(S)), target = target)
    } else if (utilopt == "rags2ridges") {
      hatOmega <- rags2ridges::ridgeP(S = S, lambda = lambda, target = target)
    }
  } else if (method == "elnet") {
    if (is.numeric(target) & length(target) == 1) {
      target <- matrix(target, ncol(S), ncol(S))
    }
    if (utilopt == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = gamma)$Omega
    } else if (utilopt == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = gamma, Target = target)$Theta
    }
  } else if (method == "clime") {
    if (utilopt == "clime_primaldual") {
      hatOmega <- clime::clime(S, lambda = lambda, sigma = TRUE, standardize = FALSE,
                               linsolver = "primaldual")$Omegalist[[1]]
    } else if (utilopt == "clime_simplex") {
      hatOmega <- clime::clime(S, lambda = lambda, sigma = TRUE, standardize = FALSE,
                               linsolver = "simplex")$Omegalist[[1]]
    } else if (utilopt == "flare") {
      hatOmega <- flare::sugm(S, lambda = lambda, method = "clime", verbose = FALSE)$icov[[1]]
    }
  } else if (method == "tiger") {
    if (utilopt == "flare") {
      hatOmega <- flare::sugm(X, lambda = lambda, method = "tiger", verbose = FALSE)$icov[[1]]
    } else if (utilopt == "huge") {
      hatOmega <- huge::huge.tiger(X, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  }
  return(hatOmega)
}


##----------------------------------------------------------------------------------------
#' Graphical lasso
#'
#' @description
#' Estimate a sparse precision matrix using a lasso (L1) penalty.
#'
#' @param S The sampmle covariance or correlation matrix, a p-by-p symmetric matrix.
#'
#' @param lambda A non-negative scalar, p-by-p matrix, or vector of length p specifying
#' the regularization parameter for penalty. In the vector case, the parameter matrix has
#' ij-th element lambda[i].
#'
#' @param ... Additional arguments passed to \code{fcstat}.
#'
#' @importFrom glassoFast glassoFast
#'
#' @return Estimated precision matrix.
#'
#' @noRd

fcstat_glasso <- function(S, lambda, ...) {
  hatOmega <- glassoFast::glassoFast(S, rho = lambda)$wi
  return(hatOmega)
}


##----------------------------------------------------------------------------------------
#' Graphical ridge
#'
#' @description
#' Estimate a precision matrix using a ridge (L2) penalty.
#'
#' @param S The sampmle covariance or correlation matrix, a p-by-p symmetric matrix.
#'
#' @param lambda A non-negative scalar specifying the regularization parameter for penalty.
#'
#' @param target A p-by-p symmetric matrix or a scalar (default = 0) serves as the value
#' for all elements, specifying the target matrix.
#'
#' @param ... Additional arguments passed to \code{fcstat}.
#'
#' @import RBGL
#' @import graph
#' @importFrom rags2ridges ridgeP
#'
#' @return Estimated precision matrix.
#'
#' @noRd

fcstat_ridge <- function(S, lambda, target = 0, ...) {
  if (is.numeric(target) & length(target) == 1) {
    target <- matrix(target, ncol(S), ncol(S))
  }
  hatOmega <- rags2ridges::ridgeP(S, lambda = lambda, target = target)
  return(hatOmega)
}


##----------------------------------------------------------------------------------------
#' Graphical elastic net
#'
#' @description
#' Estimate a sparse precision matrix using an elastic net type penalty.
#'
#' @param S The sampmle covariance or correlation matrix, a p-by-p symmetric matrix.
#'
#' @param lambda A non-negative scalar, p-by-p matrix, or vector of length p specifying
#' the regularization parameter for penalty. In the vector case, the parameter matrix has
#' ij-th element \eqn{\sqrt{lambda[i]*lambda[j]}}.
#'
#' @param gamma A scalar between 0 and 1 specifying the tuning parameter to balance L1 and
#' L2 penalties.
#'
#' @param target A p-by-p symmetric matrix or a scalar (default = 0) serves as the value
#' for all elements, specifying the target matrix.
#'
#' @param ... Additional arguments passed to \code{fcstat}.
#'
#' @importFrom GLassoElnetFast gelnet
#'
#' @return Estimated precision matrix.
#'
#' @noRd

fcstat_elnet <- function(S, lambda, gamma, target = 0, ...) {
  if (is.numeric(target) & length(target) == 1) {
    target <- matrix(target, ncol(S), ncol(S))
  }
  hatOmega <- GLassoElnetFast::gelnet(S, lambda = lambda, alpha = gamma, Target = target)$Theta
  return(hatOmega)
}


##----------------------------------------------------------------------------------------
#' Constrained L1-minimization for inverse (covariance) matrix estimation
#'
#' @description
#' Estimate a sparse precision matrix using a constrained L1 minimization approach.
#'
#' @param S The sampmle covariance or correlation matrix, a p-by-p symmetric matrix.
#'
#' @param lambda A non-negative scalar specifying the regularization parameter.
#'
#' @param utilopt A character string specifying the utility option to use for
#' \code{method = "clime"}: \itemize{
#' \item "clime_primaldual": the utility function \code{\link[clime]{clime}} from
#' the package \code{clime} with the linsolver \code{primaldual}.
#' \item "clime_simplex": the utility function \code{\link[clime]{clime}} from
#' the package \code{clime} with the linsolver \code{simplex}.
#' \item "flare": the utility function \code{\link[flare]{sugm}} from
#' the package \code{flare}.
#' }
#'
#' @param ... Additional arguments passed to \code{fcstat}.
#'
#' @importFrom clime clime
#' @importFrom flare sugm
#'
#' @return Estimated precision matrix.
#'
#' @noRd

fcstat_clime <- function(S, lambda, utilopt,...) {
  if (utilopt == "clime_primaldual") {
    hatOmega <- clime::clime(S, lambda = lambda, sigma = TRUE, standardize = FALSE,
                             linsolver = "primaldual")$Omegalist[[1]]
  } else if (utilopt == "clime_simplex") {
    hatOmega <- clime::clime(S, lambda = lambda, sigma = TRUE, standardize = FALSE,
                             linsolver = "simplex")$Omegalist[[1]]
  } else if (utilopt == "flare") {
    hatOmega <- flare::sugm(S, lambda = lambda, method = "clime", verbose = FALSE)$icov[[1]]
  }
  return(hatOmega)
}


##----------------------------------------------------------------------------------------
#' Tuning-insensitive graph estimation and regression
#'
#' @description
#' Estimate a sparse precision matrix using the tuning-insensitive graph estimation and
#' regression.
#'
#' @param X A data matrix.
#'
#' @param lambda A non-negative scalar specifying the regularization parameter.
#'
#' @importFrom huge huge.tiger
#'
#' @return Estimated precision matrix.
#'
#' @noRd

fcstat_tiger <- function(X, lambda, ...) {
  hatOmega <- huge::huge.tiger(X, lambda = lambda, verbose = FALSE)$icov[[1]]
  return(hatOmega)
}

