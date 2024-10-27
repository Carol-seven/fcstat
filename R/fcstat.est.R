#' Estimation for precision matrix
#'
#' @description
#' Provide a collection of statistical methods to obtain a series of precision matrix
#' estimates across the grid of tuning parameters.
#'
#' @param X \enumerate{
#' \item An n-by-p data matrix with sample size n and dimension p.
#' \item A p-by-p sample covariance/correlation matrix with dimension p.
#' }
#'
#' @param method A character string specifying the statistical method for estimating
#' precision matrix. Available options include: \enumerate{
#' \item "glasso": graphical lasso \insertCite{friedman2008sparse}{fcstat}.
#' \item "ridge": graphical ridge \insertCite{vanwieringen2016ridge}{fcstat}.
#' \item "elnet": graphical elastic net \insertCite{zou2005regularization}{fcstat}.
#' \item "clime": constrained L1-minimization for inverse (covariance) matrix estimation
#' \insertCite{cai2011aconstrained}{fcstat}.
#' \item "tiger": tuning-insensitive graph estimation and regression
#' \insertCite{liu2017tiger}{fcstat}, which is only applicable when \code{X} is the n-by-p
#' data matrix.
#' \item "adapt": adaptive lasso \insertCite{zou2006adaptive,fan2009network}{fcstat}.
#' \item "atan": arctangent type penalty \insertCite{wang2016variable}{fcstat}.
#' \item "exp": exponential type penalty \insertCite{wang2018variable}{fcstat}.
#' \item "mcp": minimax concave penalty \insertCite{zou2006adaptive}{fcstat}.
#' \item "scad": smoothly clipped absolute deviation \insertCite{fan2001variable,fan2009network}{fcstat}.
#' }
#'
#' @param base A character string (default = "cov") specifying the calculation base,
#' either the covariance matrix ("cov") or the correlation matrix ("cor"). This is only
#' applicable when \code{X} is the n-by-p data matrix.
#'
#' @param lambda Grid of non-negative scalars for the regularization parameter.
#' The default is \code{NULL}, which generates its own \code{lambda} sequence based on
#' \code{nlambda} and \code{lambda.min.ratio}. For \code{method = "clime"} combined with
#' \code{pkgopt = "clime"}, the \code{lambda} sequence is based on \code{nlambda},
#' \code{lambda.min} and \code{lambda.max}.
#'
#' @param nlambda An integer (default = 20) specifying the number of \code{lambda} values
#' to be generated when \code{lambda = NULL}.
#'
#' @param lambda.min.ratio A scalar (default = 0.01) specifying the fraction of
#' the maximum \code{lambda} value \eqn{\lambda_{max}} to generate the minimum
#' \code{lambda} \eqn{\lambda_{min}}. If \code{lambda = NULL}, the program automatically
#' generates a \code{lambda} grid as a sequence of length \code{nlambda} in log scale,
#' starting from \eqn{\lambda_{min}} to \eqn{\lambda_{max}}.
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
#' @param initial A p-by-p matrix or a p-by-p-by-npara (the number of all combinations of
#' \code{lambda} and \code{gamma}) array specifying the initial estimate for \code{method}
#' set to \code{"atan"}, \code{"exp"}, \code{"scad"}, and \code{"mcp"}; or specifying
#' \eqn{\tilde{\Omega}} of the adaptive weight for \code{method = "adapt"}, calculated as
#' \eqn{|\tilde{\omega}_{ij}|^{-\gamma}}, where \eqn{\tilde{\Omega} := (\tilde{\omega}_{ij})}.
#' Some options are also offered when a character string is provided (default "glasso"),
#' including:
#' \itemize{
#' \item "glasso": use the precision matrix estimate derived from the graphical lasso.
#' \item "invS": use the inverse calculation base matrix if the matrix is invertible.
#' \item "linshrink": use the precision matrix estimate derived from Ledoit-Wolf linear
#' shrinakge estimator of the population covariance matrix
#' \insertCite{ledoit2004well}{fcstat}.
#' \item "nlshrink": use the precision matrix estimate derived from Ledoit-Wolf non-linear
#' shrinakge estimator of the population covariance matrix
#' \insertCite{ledoit2015spectrum,ledoit2017numerical}{fcstat}.
#' }
#'
#' @param pkgopt A character string specifying the package option to use. The available
#' options depend on the selected method: \enumerate{
#' \item For \code{method = "glasso"}: \itemize{
#' \item "ADMMsigma": the function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "CovTools": the function from \code{\link[CovTools]{PreEst.glasso}}.
#' \item "CVglasso": the function from \code{\link[CVglasso]{CVglasso}}.
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
#' \item For \code{method} set to \code{"adapt"}, \code{"atan"}, \code{"exp"},
#' \code{"scad"}, and \code{"mcp"}: \itemize{
#' \item "glasso": the function from \code{\link[glasso]{glasso}}.
#' \item "GLassoElnetFast": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "glassoFast": the function from \code{\link[glassoFast]{glassoFast}}.
#' }
#' }
#'
#' @param cores An integer (default = 1) specifying the number of cores to use for
#' parallel execution.
#'
#' @note
#' For the method \code{tiger}, the estimation process solely relies on the raw n-by-p
#' data \code{X} and does not utilize the argument \code{base}. This argument is not
#' applicable for \code{tiger} and will have no effect if provided.
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom glasso glasso
#' @importFrom glassoFast glassoFast
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom Rdpack reprompt
#'
#' @return An object with S3 class "fcstat.est" containing the following components:
#' \describe{
#' \item{hatOmega}{A list of estimated precision matrices for \code{lambda} grid and
#' \code{gamma} grid.}
#' \item{lambda}{The actual lambda grid used in the program, corresponding to \code{hatOmega}.}
#' \item{gamma}{The actual gamma grid used in the program, corresponding to \code{hatOmega}.}
#' \item{X}{The n-by-p data matrix used in the program.}
#' \item{S}{The p-by-p calculation base matrix used in the program.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @autoglobal
#'
#' @noRd

fcstat.est <- function(
    X, method, base = "cov",
    lambda = NULL, nlambda = 20, lambda.min.ratio = 0.01,
    gamma = NULL, ## for elnet, adapt, atan, exp, mcp, scad
    initial = "glasso", ## initial estimator for atan, exp, mcp, scad; adaptive weight for adapt
    pkgopt = "glasso", ## package option
    cores = 1) {

  if (!method %in% c("glasso", "ridge", "elnet", "clime", "tiger",
                     "adapt", "atan", "exp", "mcp", "scad")) {
    stop("Error in method.
         Available options: 'glasso', 'ridge', 'elnet', 'clime', 'tiger', 'adapt', 'atan', 'exp', 'mcp', 'scad'.\n")
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
      S <- cov(X)
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
  if(all(is.na(gamma))) {
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
    } else {
      gamma <- NA
    }
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
    if (method %in% c("glasso", "ridge", "elnet", "clime")) {
      hatOmega <- foreach(k = 1:npara, .packages = "fcstat",
                          .export = c("fcstat_method")) %dopar% {
        fcstat_method(method = method, X = X, S = S,
                      lambda = parameter$lambda[k], gamma = parameter$gamma[k],
                      pkgopt = pkgopt)
      }
    } else if (method == "tiger") {
      hatOmega <- foreach(k = 1:npara, .packages = "fcstat",
                          .export = c("fcstat_method")) %dopar% {
        tryCatch(
          fcstat_method(method = method, X = X, S = S,
                        lambda = parameter$lambda[k], gamma = parameter$gamma[k],
                        pkgopt = pkgopt),
          error = function(e) {
            message("lambda[", k, "] ", parameter$lambda[k], " is too samll. Ignored!")
          }
        )
      }
    } else { ## method %in% c("adapt", "atan", "exp", "mcp", "scad")
      Omega <- gen_initial(X = X, S = S, base = base, initial = initial, lambda = parameter$lambda, pkgopt = pkgopt)
      lambda_mat <- lapply(1:npara, function(k) {
        fcstat::deriv(penalty = method, Omega = Omega[[k]],
                      lambda = parameter$lambda[k], gamma = parameter$gamma[k])
      })
      if (pkgopt == "glasso") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          glasso::glasso(s = S, rho = lambda_mat[[k]], penalize.diagonal = TRUE, start = "cold")$wi
        }
      } else if (pkgopt == "GLassoElnetFast") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          GLassoElnetFast::gelnet(S = S, lambda = lambda_mat[[k]], alpha = 1, penalize.diagonal = TRUE)$Theta
        }
      } else if (pkgopt == "glassoFast") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          glassoFast::glassoFast(S = S, rho = lambda_mat[[k]], start = "cold")$wi
        }
      }
    }

    stopCluster(cluster)

  } else {

    ## compute the precision matrix estimator hatOmega along the parameter grid
    if (method %in% c("glasso", "ridge", "elnet", "clime")) {
      hatOmega <- lapply(1:npara, function(k) {
        fcstat_method(method = method, X = X, S = S,
                      lambda = parameter$lambda[k], gamma = parameter$gamma[k],
                      pkgopt = pkgopt)
      })
    } else if (method == "tiger") {
      hatOmega <- lapply(1:npara, function(k) {
        tryCatch(
          fcstat_method(method = method, X = X, S = S,
                        lambda = parameter$lambda[k], gamma = parameter$gamma[k],
                        pkgopt = pkgopt),
          error = function(e) {
            message("lambda[", k, "] ", parameter$lambda[k], " is too samll. Ignored!")
          }
        )
      })
    } else { ## method %in% c("adapt", "atan", "exp", "mcp", "scad")
      Omega <- gen_initial(X = X, S = S, base = base, initial = initial, lambda = parameter$lambda, pkgopt = pkgopt)
      lambda_mat <- lapply(1:npara, function(k) {
        fcstat::deriv(penalty = method, Omega = Omega[[k]],
                      lambda = parameter$lambda[k], gamma = parameter$gamma[k])
      })
      if (pkgopt == "glasso") {
        hatOmega <- lapply(1:npara, function(k) {
          glasso::glasso(s = S, rho = lambda_mat[[k]], penalize.diagonal = TRUE, start = "cold")$wi
        })
      } else if (pkgopt == "GLassoElnetFast") {
        hatOmega <- lapply(1:npara, function(k) {
          GLassoElnetFast::gelnet(S = S, lambda = lambda_mat[[k]], alpha = 1, penalize.diagonal = TRUE)$Theta
        })
      } else if (pkgopt == "glassoFast") {
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
  class(result) <- c("fcstat.est")
  return(result)

}

