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
#' \code{nlambda} and \code{lambda.min.ratio}.
#'
#' @param nlambda An integer (default = 50) specifying the number of \code{lambda} values
#' to be generated when \code{lambda = NULL}.
#'
#' @param lambda.min.ratio A scalar specifying the fraction of the maximum \code{lambda}
#' value \eqn{\lambda_{max}} to generate the minimum \code{lambda} \eqn{\lambda_{min}}.
#' If \code{lambda = NULL}, the program automatically generates a \code{lambda} grid as a
#' sequence of length \code{nlambda} in log scale, starting from \eqn{\lambda_{min}} to
#' \eqn{\lambda_{max}}. The default value is 0.4 for \code{method = "clime"}, 0.1 for
#' \code{method = "tiger"} and 0.01 for other methods.
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
#' @param target A p-by-p symmetric matrix or a scalar (default = 0) serves as the value
#' for all elements, specifying the target matrix for \code{method = "ridge"} or
#' \code{method = "elnet"}.
#'
#' @param initial A p-by-p matrix or a p-by-p-by-npara (the number of all combinations of
#' \code{lambda} and \code{gamma}) array specifying the initial estimate for \code{method}
#' set to \code{"atan"}, \code{"exp"}, \code{"scad"}, or \code{"mcp"}; or specifying
#' \eqn{\tilde{\Omega}} of the adaptive weight for \code{method = "adapt"}, calculated as
#' \eqn{|\tilde{\omega}_{ij}|^{-\gamma}}, where \eqn{\tilde{\Omega} := (\tilde{\omega}_{ij})}.
#' Some options are also offered when a character string is provided (default "linshrink"),
#' including:
#' \itemize{
#' \item "glasso": use the precision matrix estimate derived from the graphical lasso.
#' \item "linshrink": use the precision matrix estimate derived from Ledoit-Wolf linear
#' shrinakge estimator of the population covariance matrix
#' \insertCite{ledoit2004well}{fcstat}.
#' \item "nlshrink": use the precision matrix estimate derived from Ledoit-Wolf non-linear
#' shrinakge estimator of the population covariance matrix
#' \insertCite{ledoit2015spectrum,ledoit2017numerical}{fcstat}.
#' \item "invS-glasso": use the inverse calculation base matrix if the matrix is
#' invertible; otherwise, use the precision matrix estimate derived from the graphical
#' lasso.
#' \item "invS-linshrink": use the inverse calculation base matrix if the matrix is
#' invertible; otherwise, use the precision matrix estimate derived from Ledoit-Wolf
#' linear shrinakge estimator of the population covariance matrix
#' \insertCite{ledoit2004well}{fcstat}.
#' \item "invS-nlshrink": use the inverse calculation base matrix if the matrix is
#' invertible; otherwise, use the precision matrix estimate derived from Ledoit-Wolf
#' non-linear shrinkage estimator of the population covariance matrix
#' \insertCite{ledoit2015spectrum,ledoit2017numerical}{fcstat}.
#' }
#'
#' @param utilopt A character string specifying the utility option to use for
#' \code{method = "clime"}: \itemize{
#' \item "clime_primaldual": utility function from the package \code{\link[clime]{clime}}
#' with the linsolver primaldual.
#' \item "clime_simplex": utility function from the package \code{\link[clime]{clime}}
#' with the linsolver simplex.
#' \item "flare": use the utility function from the package \code{\link[flare]{sugm}}.
#' }
#'
#' @note
#' For the method \code{tiger}, the estimation process solely relies on the raw n-by-p
#' data \code{X} and does not utilize the argument \code{base}. This argument is not
#' applicable for \code{tiger} and will have no effect if provided.
#'
#' @import foreach
#' @importFrom glassoFast glassoFast
#' @importFrom Rdpack reprompt
#'
#' @return An object with S3 class "fcstat.est" containing the following components: \describe{
#' \item{hatOmega}{A list of estimated precision matrices for \code{lambda} grid and \code{gamma} grid.}
#' \item{X}{The actual n-by-p data matrix used in the program.}
#' \item{S}{The p-by-p calculation base matrix used in the program.}
#' \item{method}{The statistical method used for estimating precision matrix.}
#' \item{base}{The calculation base.}
#' \item{lambda}{The actual lambda grid used in the program, corresponding to \code{hatOmega}.}
#' \item{gamma}{The actual gamma grid used in the program, corresponding to \code{hatOmega}.}
#' \item{target}{The target matrix.}
#' \item{intial}{The initial estimate or \eqn{\tilde{\Omega}} of the adaptive weight.}
#' \item{utilopt}{The utility option.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export

fcstat.est <- function(
    X, method, base = "cov",
    lambda = NULL, nlambda = 50, lambda.min.ratio = NULL,
    gamma = NULL, ## for elnet, adapt, atan, exp, mcp, scad
    target = 0, ## for ridge, elnet
    initial = "linshrink", ## initial estimator for atan, exp, mcp, scad; adaptive weight for adapt
    utilopt = "flare") { ## utility option for clime

  if (!method %in% c("glasso", "ridge", "elnet", "clime", "tiger",
                     "adapt", "atan", "exp", "mcp", "scad")) {
    stop("Error in method.
         Available options: glasso, ridge, elnet, clime, tiger, adapt, atan, exp, mcp, scad")
  }

  ## dimensionality
  p <- ncol(X)

  ## sample covariance/correlation matrix
  if (isSymmetric(X)) {
    S <- X
    X <- NULL
  } else {
    X <- scale(X)
    S <- eval(parse(text =  paste0(base, "(X)")))
  }

  ## lambda grid
  if(is.null(lambda)) {
    if (method == "clime") { ## pkg:flare
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio <- 0.4
      }
      lambda.max.tmp <- abs(range(S - diag(diag(S))))
      lambda.max <- ifelse(min(lambda.max.tmp) == 0, max(lambda.max.tmp), min(lambda.max.tmp))
      # ## pkg:clime
      # if (is.null(lambda.min)) {
      #   lambda.min <- ifelse(n > p, 1e-4, 1e-2)
      # }
      # if (is.null(lambda.max)) {
      #   lambda.max <- 0.8
      # }
    } else if (method == "tiger") { ## pkg:huge
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio <- 0.1
      }
      lambda.max <- max(abs(range(S - diag(p))))
      # ## pkg:flare
      # if (is.null(lambda.min.ratio)) {
      #   lambda.min.ratio <- 0.4
      # }
      # lambda.max <- pi*sqrt(log(p)/n)
      # lambda.min <- lambda.min.ratio*lambda.max
    } else {
      if (is.null(lambda.min.ratio)) {
        lambda.min.ratio <- 0.01
      }
      lambda.max <- max(abs(range(S - diag(p)))) ## max(max(S-diag(p)), -min(S-diag(p)))
    }
    lambda.min <- lambda.min.ratio*lambda.max
    lambda <- exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  }
  nlambda <- length(lambda)

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

  ## compute the precision matrix estimator hatOmega along the parameter grid
  if (method %in% c("glasso", "ridge", "elnet", "clime", "tiger")) {
    hatOmega <- foreach(k = 1:nrow(parameter)) %do% {
      eval(parse(text = paste0(
        "fcstat_", method,
        "(X = X, S = S, lambda = parameter$lambda[k], gamma = parameter$gamma[k], target = target, utilopt = utilopt)"
      )))
    }
  } else if (method %in% c("adapt", "atan", "exp", "mcp", "scad")) {
    if (all(grepl("^invS-", initial))) {
      Omega <- tryCatch({
        gen_initial(X, S, base, initial = "invS", parameter)
      }, error = function(e) {
         gen_initial(X, S, base, initial = sub("^invS-(.*)", "\\1", initial), parameter)
      })
    } else {
      Omega <- gen_initial(X, S, base, initial = initial, parameter)
    }

    hatOmega <- mapply(function(x, y, z) {
      lambda_mat <- deriv(penalty = method, Omega = z, lambda = x, gamma = y)
      glassoFast::glassoFast(S, rho = lambda_mat)$wi
    }, x = parameter$lambda, y = parameter$gamma, z = Omega, SIMPLIFY = FALSE)
  }

  result <- list(hatOmega = hatOmega,
                 X = X, S = S,
                 method = method, base = base,
                 lambda = parameter$lambda, gamma = parameter$gamma,
                 target = if (method %in% c("ridge", "elnet")) target else NULL,
                 initial = if (method %in% c("adapt", "atan", "exp", "mcp", "scad")) initial else NULL,
                 utilopt = if (method == "clime") utilopt else NULL)
  class(result) <- c("fcstat.est")
  return(result)

}

