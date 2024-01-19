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
#' @import glassoFast
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
#' @import clime
#'
#' @return Estimated precision matrix.
#'
#' @noRd

fcstat_clime <- function(S, lambda, ...) {
  hatOmega <- clime::clime(S, lambda = lambda, sigma = TRUE, standardize = FALSE)$Omegalist[[1]]
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


##----------------------------------------------------------------------------------------
#' Statistical methods for estimating precision matrix
#'
#' @description
#' Provide a collection of statistical methods to estimate a precision matrix.
#'
#' @param X An n-by-p data matrix with sample size n and dimension p.
#'
#' @param method A character string specifying the statistical methods for estimating
#' precision matrix. Available options include: \enumerate{
#' \item "glasso": graphical lasso \insertCite{friedman2008sparse}{fcstat}.
#' \item "ridge": graphical ridge \insertCite{vanwieringen2016ridge}{fcstat}.
#' \item "elnet": graphical elastic net \insertCite{zou2005regularization}{fcstat}.
#' \item "clime": constrained L1-minimization for inverse (covariance) matrix estimation
#' \insertCite{cai2011aconstrained}{fcstat}.
#' \item "tiger": tuning-insensitive graph estimation and regression
#' \insertCite{liu2017tiger}{fcstat}.
#' \item "adapt": adaptive lasso \insertCite{zou2006adaptive,fan2009network}{fcstat}.
#' \item "atan": arctangent type penalty \insertCite{wang2016variable}{fcstat}.
#' \item "exp": exponential type penalty \insertCite{wang2018variable}{fcstat}.
#' \item "mcp": minimax concave penalty \insertCite{zou2006adaptive}{fcstat}.
#' \item "scad": smoothly clipped absolute deviation \insertCite{fan2001variable,fan2009network}{fcstat}.
#' }
#'
#' @param base A character string (default = "cov") specifying the calculation base,
#' either the covariance matrix ("cov") or the correlation matrix ("cor").
#'
#' @param approach A character string (default = "smp") specifying the approach used to
#' compute the calculation base. Available options include: \enumerate{
#' \item "smp": sample (conventional approach).
#' \item "lin": linear shrinkage \insertCite{ledoit2004well}{fcstat}.
#' \item "nlminb": non-linear shrinkage using the optimization routine "nlminb"
#' \insertCite{ledoit2015spectrum,ledoit2017numerical}{fcstat}.
#' \item "nloptr": non-linear shrinkage using the optimization routine "nloptr"
#' \insertCite{ledoit2015spectrum,ledoit2017numerical}{fcstat}.
#' }
#'
#' @param lambda Grid of non-negative scalars for the regularization parameter.
#' The default is \code{NULL}, which generates its own \code{lambda} sequence based on
#' \code{nlambda} and \code{lambda.min.ratio}. For \code{method = "clime"}, the
#' \code{lambda} sequence is based on \code{nlambda}, \code{lambda.min} and \code{lambda.max}.
#'
#' @param nlambda An integer (default = 50) specifying the number of \code{lambda} values
#' to be generated when \code{lambda = NULL}.
#'
#' @param lambda.min.ratio A scalar specifying the fraction of the maximum \code{lambda}
#' value \eqn{\lambda_{max}} to generate the minimum \code{lambda} \eqn{\lambda_{min}}.
#' If \code{lambda = NULL}, the program automatically generates a \code{lambda} grid as a
#' sequence of length \code{nlambda} in log scale, starting from \eqn{\lambda_{min}} to
#' \eqn{\lambda_{max}}, expect for \code{method = "clime"}. The default value is 0.4 for
#' \code{method = "tiger"} and 0.01 for others.
#'
#' @param lambda.min A scalar specifying the minimum value of program generated
#' \code{lambda} grid for \code{method = "clime"}. Default is 1e-4 (\eqn{n>p}) or 1e-2
#' (\eqn{n<p}).
#'
#' @param lambda.max A scalar (default = 0.8) specifying the maximum value of program
#' generated \code{lambda} grid for \code{method = "clime"}.
#'
#' @param gamma A scalar specifying the hyperparameter for the chosen \code{method}.
#' Default values: \enumerate{
#' \item "elnet": a sequence from 0.1 to 0.9 with increments of 0.1
#' \item "adapt": 0.5
#' \item "atan": 0.005
#' \item "exp": 0.01
#' \item "scad": 3.7
#' \item "mcp": 3
#' }
#'
#' @param target A p-by-p symmetric matrix or a scalar (default = 0) serves as the value
#' for all elements, specifying the target matrix for \code{method = "ridge"} or
#' \code{method = "elnet"}.
#'
#' @param initial A character string (default = "glasso") specifying \itemize{
#' \item The initial estimate for \code{method} set to \code{"atan"}, \code{"exp"},
#' \code{"scad"}, or \code{"mcp"}. \itemize{
#' \item "glasso": use the precision matrix estimate derived from the graphical lasso.
#' \item "Fan": use the inverse calculation base matrix for \eqn{p < n}; use the precision
#' matrix estimate derived from the graphical lasso for \eqn{p \geq n}.
#' }
#' \item The adaptive weight for \code{method = "adapt"}, calculated as
#' \eqn{|\tilde{\omega}_{ij}|^{-\gamma}}, where \eqn{\tilde{\Omega} := (\tilde{\omega}_{ij})}.
#' \itemize{
#' \item "glasso": use the precision matrix estimate derived from the graphical lasso as
#' \eqn{\tilde{\Omega}}.
#' \item "Fan": use the inverse calculation base matrix as \eqn{\tilde{\Omega}} for
#' \eqn{p < n};  use the precision matrix estimate derived from the graphical lasso as
#' \eqn{\tilde{\Omega}} for \eqn{p \geq n}.
#' }
#' }
#'
#' @param crit A string specifying the tuning parameter selection criterion to use.
#' Available options include: \enumerate{
#' \item "AIC": Akaike information criterion \insertCite{akaike1973information}{fcstat}.
#' \item "BIC": Bayesian information criterion \insertCite{schwarz1978estimating}{fcstat}.
#' \item "EBIC": extended Bayesian information criterion \insertCite{foygel2010extended}{fcstat}.
#' \item "HBIC": high dimensional Bayesian information criterion \insertCite{wang2013calibrating,fan2017high}{fcstat}.
#' }
#'
#' @param ebic.tuning A scalar (default = 0.5) specifying the tuning parameter to
#' calculate when \code{crit = "EBIC"}.
#'
#' @param scaling.ratio A scalar (default = 0.5) specifying the scaling ratio to adjust
#' the \code{lambda} grid by multiplying the minimum and maximum with the scaling ratio
#' for the case where the optimal regularization parameter \code{lambda} is found at the
#' first index of \code{lambda} grid, make the process run iteratively with the expanded
#' grid until the optimal parameter is no longer at the first index.
#'
#' @import foreach
#' @import glassoFast
#' @importFrom Rdpack reprompt
#'
#' @return A list containing the following components: \describe{
#' \item{hatOmega_opt}{Estimated precision matrix.}
#' \item{lambda_opt}{Optimal regularization parameter.}
#' \item{gamma_opt}{Optimal hyperparameter.}
#' \item{score_opt}{Information criterion score.}
#' \item{lambda}{Actual lambda grid used in the program.}
#' \item{gamma}{Actual gamma grid used in the program.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export

fcstat <- function(X, method, base = "cov", approach = "smp",
                   lambda = NULL, nlambda = 50, lambda.min.ratio = NULL,
                   lambda.min = NULL, lambda.max = NULL, ## for clime
                   gamma = NA, ## for elnet, adapt (with adapt.weight = "Fan"), atan, exp, scad, mcp
                   target = 0, ## for ridge, elnet
                   initial = "glasso", ## initial estimator for atan, exp, scad, mcp; adaptive weight for adapt
                   crit = "BIC", ebic.tuning = 0.5,
                   scaling.ratio = 0.5) {

  if (!method %in% c("glasso", "ridge", "elnet", "clime", "tiger",
                     "adapt", "atan", "exp", "mcp", "scad")) {
    stop("Error in method.
         Available options: glasso, ridge, elnet, clime, tiger, adapt, atan, exp, mcp, scad")
  }

  ## sample size
  n <- nrow(X)
  ## dimensionality
  p <- ncol(X)

  ## the calculation base S is a customized combination of 'base' and 'approach'
  if (approach == "smp") {
    S <- eval(parse(text =  paste0(base, "(X)")))
  } else if (approach %in% c("lin", "nlminb", "nloptr")) {
    S <- ledoit_wolf_est(X, method = approach, res = base)
  }

  ## lambda grid
  if(is.null(lambda)) {
    if (method == "clime") { ## pkg:clime
      if (is.null(lambda.min)) {
        lambda.min <- ifelse(n > p, 1e-4, 1e-2)
      }
      if (is.null(lambda.max)) {
        lambda.max <- 0.8
      }
    } else if (method == "tiger") { ## pkg:flare
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio <- 0.4
      }
      lambda.max <- pi*sqrt(log(p)/n)
      lambda.min <- lambda.min.ratio*lambda.max
    } else {
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio <- 0.01
      }
      lambda.max <- max(abs(S - diag(p))) ## max(max(S-diag(p)), -min(S-diag(p)))
      lambda.min <- lambda.min.ratio*lambda.max
    }
    lambda <- exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  }
  nlambda <- length(lambda)

  ## gamma grid
  if(is.null(gamma)) {
    if (method == "elnet") {
      gamma <- seq(0.1, 0.9, 0.1)
    } else if (method == "adapt") {
      gamma <- 0.5
    } else if (method == "atan") {
      gamma <- 0.005
    } else if (method == "exp") {
      gamma <- 0.01
    } else if (method == "scad") {
      gamma <- 3.7
    } else if (method == "mcp") {
      gamma <- 3
    }
  }

  repeat{
    ## run the analysis with the current parameter grid

    ## parameter grid combination
    parameter <- expand.grid(lambda = lambda, gamma = gamma)

    ## compute the precision matrix estimator hatOmega along the parameter grid
    if (method %in% c("glasso", "ridge", "elnet", "clime", "tiger")) {
      hatOmega <- foreach(k = 1:nrow(parameter)) %do% {
        eval(parse(text = paste0(
          "fcstat_", method,
          "(X = X, S = S, lambda = parameter$lambda[k], gamma = parameter$gamma[k], target = target)"
        )))
      }
    } else if (method %in% c("adapt", "atan", "exp", "mcp", "scad")) {
      if (initial == "glasso") {
        hatOmega <- foreach(k = 1:nrow(parameter)) %do% {
          Omega <- glassoFast::glassoFast(S, rho = parameter$lambda[k])$wi
          lambda_mat <- eval(parse(text = paste0(
            "deriv(penalty = ", method, "Omega = Omega, lambda = parameter$lambda[k], gamma = parameter$gamma[k])"
          )))
          glassoFast::glassoFast(S, rho = lambda_mat)$wi
        }
      } else if (initial == "Fan") {
        ## Fan, J., Y. Feng, and Y. Wu (2009). Network exploration via the adaptive LASSO
        ## and SCAD penalties. The Annals of Applied Statistics, 3 (2), 521–541.
        hatOmega <- foreach(k = 1:nrow(parameter)) %do% {
          if (p < n) {
            Omega <- solve(S)
          } else {
            Omega <- glassoFast::glassoFast(S, rho = parameter$lambda[k])$wi
          }
          lambda_mat <- eval(parse(text = paste0(
            "deriv(penalty = ", method, "Omega = Omega, lambda = parameter$lambda[k], gamma = parameter$gamma[k])"
          )))
          glassoFast::glassoFast(S, rho = lambda_mat)$wi
        }
      }
    }

    ## select the optimal parameters among a set of possible values
    eval(parse(text = paste0(
      "score <- foreach(k = 1:nrow(parameter), .combine = 'c') %do% {",
      "criterion(hatOmega = hatOmega[[k]], S = S, n = n, crit = '", crit,
      "', ebic.tuning = ebic.tuning)}"
    )))

    ## find the optimal parameter
    index <- which.min(score)

    ## if the optimal lambda is not at the beginning, exit the loop
    if (nlambda == 1 | parameter$lambda[index] != min(lambda)) {
      break
    }

    message("Optimal lambda was at the beginning of the lambda grid. Rerunning with an expanded lambda range.")
    ## if the lambda index is at the beginning, rerun with an expanded range
    ## adjust lambda.max and lambda.min by multiplying them with the scaling ratio
    lambda.max <- lambda.max*scaling.ratio
    lambda.min <- lambda.min*scaling.ratio
    lambda <- exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  }

  result <- list(hatOmega_opt = hatOmega[[index]],
                 lambda_opt = parameter$lambda[index],
                 gamma_opt = parameter$gamma[index],
                 score_opt = score[index],
                 lambda = unique(parameter$lambda),
                 gamma = unique(parameter$gamma),
                 score = score)
  return(result)

}

utils::globalVariables(c("k", "score"))

