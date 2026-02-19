#' Estimation for Precision Matrix
#'
#' @description
#' Provide a collection of statistical methods to obtain a series of precision
#' matrix estimates across the grid of tuning parameters.
#'
#' @param X \enumerate{
#' \item An \eqn{n \times p} data matrix with sample size \eqn{n} and
#' dimension \eqn{p}.
#' \item A \eqn{p \times p} sample covariance/correlation matrix with
#' dimension \eqn{p}.
#' }
#'
#' @param method A character string specifying the statistical method for
#' estimating precision matrix. Available options include: \enumerate{
#' \item "glasso": Graphical lasso \insertCite{friedman2008sparse}{spice}.
#' \item "ridge": Graphical ridge \insertCite{vanwieringen2016ridge}{spice}.
#' \item "elnet": Graphical elastic net \insertCite{zou2005regularization}{spice}.
#' \item "clime": Constrained L1-minimization for inverse (covariance) matrix
#' estimation \insertCite{cai2011aconstrained}{spice}.
#' \item "tiger": Tuning-insensitive graph estimation and regression
#' \insertCite{liu2017tiger}{spice}, which is only applicable when \code{X} is
#' the \eqn{n \times p} data matrix.
#' \item "adapt": Adaptive lasso
#' \insertCite{zou2006adaptive,fan2009network}{spice}.
#' \item "atan": Arctangent type penalty \insertCite{wang2016variable}{spice}.
#' \item "exp": Exponential type penalty \insertCite{wang2018variable}{spice}.
#' \item "mcp": Minimax concave penalty \insertCite{zou2006adaptive}{spice}.
#' \item "scad": Smoothly clipped absolute deviation
#' \insertCite{fan2001variable,fan2009network}{spice}.
#' }
#'
#' @param base A character string (default = "cov") specifying the calculation
#' base: \enumerate{
#' \item "cov": The covariance matrix.
#' \item "cor": The correlation matrix.
#' }
#' This is only applicable when \code{X} is the \eqn{n \times p} data matrix.
#'
#' @param lambda A non-negative numeric vector specifying the grid for
#' the regularization parameter.
#' The default is \code{NULL}, which generates its own \code{lambda} sequence
#' based on \code{nlambda} and \code{lambda.min.ratio}.
#' For \code{method = "clime"} combined with \code{pkg = "clime"},
#' the \code{lambda} sequence is based on \code{nlambda}, \code{lambda.min} and
#' \code{lambda.max}.
#'
#' @param nlambda An integer (default = 20) specifying the number of
#' \code{lambda} values to generate when \code{lambda = NULL}.
#'
#' @param lambda.min.ratio A numeric value > 0 (default = 0.01) specifying
#' the fraction of the maximum \code{lambda} value \eqn{\lambda_{max}} to
#' generate the minimum \code{lambda} \eqn{\lambda_{min}}.
#' If \code{lambda = NULL}, a \code{lambda} grid of length \code{nlambda} is
#' automatically generated on a log scale, ranging from \eqn{\lambda_{max}}
#' down to \eqn{\lambda_{min}}.
#'
#' @param gamma A numeric value specifying the additional parameter for
#' the chosen \code{method}. Default values:
#' \enumerate{
#' \item "elnet": A sequence from 0.1 to 0.9 with increments of 0.1
#' \item "adapt": 0.5
#' \item "atan": 0.005
#' \item "exp": 0.01
#' \item "mcp": 3
#' \item "scad": 3.7
#' }
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
#' The available options depend on the selected method:
#' \enumerate{
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
#' \item For \code{method} set to \code{"adapt"}, \code{"atan"}, \code{"exp"},
#' \code{"scad"}, and \code{"mcp"}: \itemize{
#' \item "glasso": The function from \code{\link[glasso]{glasso}}.
#' \item "GLassoElnetFast": The function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "glassoFast": The function from \code{\link[glassoFast]{glassoFast}}.
#' }
#' }
#'
#' @param cores An integer (default = 1) specifying the number of cores to use
#' for parallel execution.
#'
#' @return
#' An object with S3 class "spice_internal" containing the following components:
#' \describe{
#' \item{hatOmega}{A list of estimated precision matrices for \code{lambda} grid
#' and \code{gamma} grid.}
#' \item{lambda}{The actual lambda grid used in the program, corresponding to
#' \code{hatOmega}.}
#' \item{gamma}{The actual gamma grid used in the program, corresponding to
#' \code{hatOmega}.}
#' \item{X}{The \eqn{n \times p} data matrix used in the program.}
#' \item{S}{The \eqn{p \times p} calculation base matrix used in the program.}
#' }
#'
#' @note
#' For the method \code{tiger}, the estimation process solely relies on the raw
#' \eqn{n \times p} data \code{X} and does not utilize the argument \code{base}.
#' This argument is not applicable for \code{tiger} and will have no effect
#' if provided.
#'
#' @references
#' \insertAllCited{}
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom glasso glasso
#' @importFrom GLassoElnetFast gelnet
#' @importFrom glassoFast glassoFast
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom Rdpack reprompt
#'
#' @autoglobal
#'
#' @noRd

spice_internal <- function(
    X, method, base = "cov",
    lambda = NULL, nlambda = 20, lambda.min.ratio = 0.01,
    gamma = NULL, ## for elnet, adapt, atan, exp, mcp, scad
    initial = "glasso", ## initial estimator for atan, exp, mcp, scad; adaptive weight for adapt
    pkg = "glasso", ## package option
    cores = 1) {

  if (!method %in% c("glasso", "ridge", "elnet", "clime", "tiger",
                     "adapt", "atan", "exp", "mcp", "scad")) {
    stop('Error in `method`!
         Available options: "glasso", "ridge", "elnet", "clime", "tiger", "adapt", "atan", "exp", "mcp", "scad".')
  }

  ## dimensionality
  p <- ncol(X)
  n <- nrow(X)

  ## sample covariance/correlation matrix
  if (isSymmetric(X)) {
    S <- X
    X <- NULL
  } else {
    if (base == "cov") {
      S <- (n-1)/n*cov(X)
    } else if (base == "cor") {
      S <- cor(X)
    }
  }

  ## lambda grid
  if(is.null(lambda)) {
    S_adjusted <- S - diag(p)
    lambda.max <- max(max(S_adjusted), -min(S_adjusted))
    lambda.min <- lambda.min.ratio*lambda.max
    lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }

  ## gamma grid
  if (all(is.na(gamma))) {
    gamma <- switch(method, "elnet" = seq(0.1, 0.9, 0.1), "adapt" = 0.5,
                    "atan" = 0.005, "exp"  = 0.01,
                    "mcp"  = 3, "scad" = 3.7, NA)
  }

  ## parameter grid combination
  parameter <- expand.grid(lambda = unique(lambda), gamma = unique(gamma))

  npara <- nrow(parameter)

  if (npara > 1 & cores > 1) {

    ## CPU cores
    num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
    if (cores > num.cores) {
      cat("The number of available CPU cores is ", num.cores, "!\n", sep = "")
    }
    if (cores > npara) {
      cores <- npara
    }
    cluster <- makeCluster(cores)
    registerDoParallel(cluster)

    ## compute the precision matrix estimator hatOmega along the parameter grid
    if (method %in% c("glasso", "ridge", "elnet", "clime", "tiger")) {
      hatOmega <- foreach(k = 1:npara, .packages = "spice",
                          .export = c("spice_method")) %dopar% {
        spice_method(method = method, X = X, S = S,
                     lambda = parameter$lambda[k], gamma = parameter$gamma[k],
                     pkg = pkg)
      }
    } else { ## method %in% c("adapt", "atan", "exp", "mcp", "scad")
      Omega <- gen_initial(X = X, S = S, base = base, initial = initial,
                           lambda = parameter$lambda, pkg = pkg)
      lambda_mat <- lapply(1:npara, function(k) {
        spice::deriv(penalty = method, Omega = Omega[[k]],
                     lambda = parameter$lambda[k], gamma = parameter$gamma[k])
      })
      if (pkg == "glasso") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          glasso::glasso(s = S, rho = lambda_mat[[k]], penalize.diagonal = TRUE, start = "cold")$wi
        }
      } else if (pkg == "GLassoElnetFast") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          GLassoElnetFast::gelnet(S = S, lambda = lambda_mat[[k]], alpha = 1, penalize.diagonal = TRUE)$Theta
        }
      } else if (pkg == "glassoFast") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          glassoFast::glassoFast(S = S, rho = lambda_mat[[k]], start = "cold")$wi
        }
      }
    }

    stopCluster(cluster)

  } else {

    ## compute the precision matrix estimator hatOmega along the parameter grid
    if (method %in% c("glasso", "ridge", "elnet", "clime", "tiger")) {
      hatOmega <- lapply(1:npara, function(k) {
        spice_method(method = method, X = X, S = S,
                     lambda = parameter$lambda[k], gamma = parameter$gamma[k],
                     pkg = pkg)
      })
    } else { ## method %in% c("adapt", "atan", "exp", "mcp", "scad")
      Omega <- gen_initial(X = X, S = S, base = base, initial = initial,
                           lambda = parameter$lambda, pkg = pkg)
      lambda_mat <- lapply(1:npara, function(k) {
        spice::deriv(penalty = method, Omega = Omega[[k]],
                     lambda = parameter$lambda[k], gamma = parameter$gamma[k])
      })
      if (pkg == "glasso") {
        hatOmega <- lapply(1:npara, function(k) {
          glasso::glasso(s = S, rho = lambda_mat[[k]], penalize.diagonal = TRUE, start = "cold")$wi
        })
      } else if (pkg == "GLassoElnetFast") {
        hatOmega <- lapply(1:npara, function(k) {
          GLassoElnetFast::gelnet(S = S, lambda = lambda_mat[[k]], alpha = 1, penalize.diagonal = TRUE)$Theta
        })
      } else if (pkg == "glassoFast") {
        hatOmega <- lapply(1:npara, function(k) {
          glassoFast::glassoFast(S = S, rho = lambda_mat[[k]], start = "cold")$wi
        })
      }
    }
  }

  result <- list(hatOmega = hatOmega,
                 lambda = parameter$lambda,
                 gamma = parameter$gamma,
                 X = X,
                 S = S)
  class(result) <- c("spice.internal")
  return(result)

}

