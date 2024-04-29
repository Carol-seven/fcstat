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
#' Estimate a precision matrix using a ridge penalty.
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
#' @param ... Additional arguments passed to \code{fcstat}.
#'
#' @importFrom flare sugm
#'
#' @return Estimated precision matrix.
#'
#' @noRd

fcstat_clime <- function(S, lambda, ...) {
  hatOmega <- flare::sugm(S, lambda = lambda, method = "clime", verbose = FALSE)$icov[[1]]
  # clime::clime(S, lambda = lambda, sigma = TRUE, standardize = FALSE)$Omegalist[[1]]
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
#' @importFrom flare sugm
#'
#' @return Estimated precision matrix.
#'
#' @noRd

fcstat_tiger <- function(X, lambda, ...) {
  hatOmega <- flare::sugm(X, lambda = lambda, method = "tiger", verbose = FALSE)$icov[[1]]
  return(hatOmega)
}

