#' Sparse Precision (Inverse Covariance) Estimation
#'
#' @description
#' Provide a collection of statistical methods to estimate a precision matrix.
#'
#' @param X \enumerate{
#' \item An \eqn{n \times p} data matrix with sample size \eqn{n} and
#' dimension \eqn{p}.
#' \item A \eqn{p \times p} sample covariance/correlation matrix with
#' dimension \eqn{p}.
#' }
#'
#' @param method A character string specifying the statistical method for
#' estimating precision matrix. Available options include:
#' \enumerate{
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
#' base:
#' \enumerate{
#' \item "cov": The covariance matrix.
#' \item "cor": The correlation matrix.
#' }
#' This is only applicable when \code{X} is the \eqn{n \times p} data matrix.
#'
#' @param n An integer (default = NULL) specifying the sample size.
#' This is only required when the input matrix \code{X} is a \eqn{p \times p}
#' sample covariance/correlation matrix with dimension \eqn{p}.
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
#' The available options depend on the chosen \code{method}:
#' \enumerate{
#' \item For \code{method = "glasso"}:
#' \itemize{
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
#' \item For \code{method = "ridge"}:
#' \itemize{
#' \item "ADMMsigma": The function from \code{\link[ADMMsigma]{RIDGEsigma}}.
#' \item "GLassoElnetFast": The function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "porridge": The function from \code{\link[porridge]{ridgePgen}}.
#' \item "rags2ridges": The function from \code{\link[rags2ridges]{ridgeP}}.
#' }
#' \item For \code{method = "elnet"}:
#' \itemize{
#' \item "ADMMsigma": The function from \code{\link[ADMMsigma]{ADMMsigma}}.
#' \item "GLassoElnetFast": The function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' }
#' \item For \code{method = "clime"}:
#' \itemize{
#' \item "clime": The function from \code{\link[clime]{clime}}.
#' \item "flare": The function from \code{\link[flare]{sugm}}.
#' }
#' \item For \code{method = "tiger"}:
#' \itemize{
#' \item "flare": The function from \code{\link[flare]{sugm}}.
#' \item "huge": The function from \code{\link[huge]{huge.tiger}}.
#' }
#' \item For \code{method} set to \code{"adapt"}, \code{"atan"}, \code{"exp"},
#' \code{"scad"}, and \code{"mcp"}:
#' \itemize{
#' \item "glasso": The function from \code{\link[glasso]{glasso}}.
#' \item "GLassoElnetFast": The function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "glassoFast": The function from \code{\link[glassoFast]{glassoFast}}.
#' }
#' }
#'
#' @param crit A character string (default = "CV") specifying the parameter
#' selection method to use. Available options include:
#' \enumerate{
#' \item "AIC": Akaike information criterion
#' \insertCite{akaike1973information}{spice}.
#' \item "BIC": Bayesian information criterion
#' \insertCite{schwarz1978estimating}{spice}.
#' \item "EBIC": extended Bayesian information criterion
#' \insertCite{foygel2010extended}{spice}.
#' \item "HBIC": high dimensional Bayesian information criterion
#' \insertCite{wang2013calibrating,fan2017high}{spice}.
#' \item "CV": k-fold cross validation with negative log-likelihood loss.
#' }
#'
#' @param fold An integer (default = 5) specifying the number of folds used for
#' \code{crit = "CV"}.
#'
#' @param ebic.tuning A numeric value in [0, 1] (default = 0.5) specifying
#' the tuning parameter to calculate for \code{crit = "EBIC"}.
#'
#' @param cores An integer (default = 1) specifying the number of cores to use
#' for parallel execution.
#'
#' @return
#' An object with S3 class \code{"spice"} containing the following components:
#' \describe{
#' \item{hatOmega_opt}{The estimated precision matrix.}
#' \item{lambda_opt}{The optimal regularization parameter.}
#' \item{gamma_opt}{The optimal hyperparameter.}
#' \item{hatOmega}{A list of estimated precision matrices for \code{lambda} grid
#' and \code{gamma} grid.}
#' \item{lambda}{The actual lambda grid used in the program, corresponding to
#' \code{hatOmega}.}
#' \item{gamma}{The actual gamma grid used in the program, corresponding to
#' \code{hatOmega}.}
#' \item{CV.loss}{Matrix of CV losses, with rows for CV folds and columns for
#' parameter combinations, when \code{crit = "CV"}.}
#' \item{IC.score}{The information criterion score for each parameter
#' combination when \code{crit} is set to \code{"AIC"}, \code{"BIC"},
#' \code{"EBIC"}, or \code{"HBIC"}.}
#' }
#'
#' @note
#' For \code{method = "tiger"}, the estimation process solely relies on the raw
#' \eqn{n \times p} data \code{X} and does not utilize the argument \code{base}.
#' This argument is not applicable for \code{method = "tiger"} and will have
#' no effect if provided.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom stats cor cov sd
#' @importFrom Rdpack reprompt
#'
#' @autoglobal
#'
#' @export

spice <- function(
    X, method, base = "cov", n = NULL,
    lambda = NULL, nlambda = 20, lambda.min.ratio = 0.01,
    gamma = NA, ## for elnet, adapt, atan, exp, mcp, scad
    initial = "glasso", ## initial estimator for atan, exp, mcp, scad; adaptive weight for adapt
    pkg = "glasso", ## package option
    crit = "CV", fold = 5, ebic.tuning = 0.5,
    cores = 1) {

  est.obj <- spice_internal(
    X = X, method = method, base = base,
    lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
    gamma = gamma, initial = initial, pkg = pkg,
    cores = cores)

  npara <- length(est.obj$hatOmega)

  if (npara > 1) {

    lambda <- est.obj$lambda
    gamma <- est.obj$gamma
    hatOmega <- est.obj$hatOmega

    if (crit == "CV") {

      if(is.null(X)) {
        stop("CV requires the n-by-p data matrix!")
      }

      n <- nrow(X)
      index <- sample(n)

      ## initialize the loss
      CV.loss <- matrix(NA, fold, npara)

      for (j in 1:fold) {

        ## indices for test and training sets; training and test sets
        ind.test <- index[((j-1)*floor(n/fold)+1):(j*floor(n/fold))]
        X.test <- X[ind.test, , drop = FALSE]
        ind.train <- index[index != ind.test]
        X.train <- X[ind.train, , drop = FALSE]

        ## sample covariance/correlation matrix
        if (base == "cov") {
          n.test <- length(ind.test)
          S.test <- (n.test-1)/n.test*cov(X.test)
        } else if (base == "cor") {
          S.test <- cor(X.test)
        }

        ## compute the precision matrix estimate
        cvlist <- spice_internal(X = X.train, method = method, base = base,
                                 lambda = lambda, gamma = gamma,
                                 initial = initial, pkg = pkg,
                                 cores = cores)

        ## loss: negative log-likelihood
        for (k in 1:npara) {
          CV.loss[j,k] <- - determinant(cvlist$hatOmega[[k]], logarithm = TRUE)$modulus[1] + sum(diag(S.test%*%cvlist$hatOmega[[k]]))
          # log(det(cvlist$hatOmega[[k]])) determinant(cvlist$hatOmega[[k]], logarithm = TRUE)$modulus[1]
        }
      }

      CV.loss[!is.finite(CV.loss)] <- Inf

      ## the mean and sd of the k-fold loss for each parameter grid value
      loss.mean <- apply(CV.loss, 2, mean)

      ## find the optimal parameter
      index <- which.min(loss.mean)

      result <- list(hatOmega_opt = hatOmega[[index]],
                     lambda_opt = lambda[index],
                     gamma_opt = gamma[index],
                     hatOmega = hatOmega,
                     lambda = lambda,
                     gamma = gamma,
                     CV.loss = CV.loss)

    } else {

      S <- est.obj$S

      if (is.null(X)) {
        if (is.null(n)) {
          stop("The input 'X' is the p-by-p sample covariance/correlation matrix.
               The selection requires the sample size 'n'.\n")
        }
      } else {
        n <- nrow(X)
      }

      ## select the optimal parameters among a set of possible values
      IC.score <- sapply(1:npara, function (k) {
        criterion(hatOmega = hatOmega[[k]], S = S, n = n, crit = crit, ebic.tuning = ebic.tuning)
      })

      ## find the optimal parameter
      index <- which.min(IC.score)

      result <- list(hatOmega_opt = hatOmega[[index]],
                     lambda_opt = lambda[index],
                     gamma_opt = gamma[index],
                     hatOmega = hatOmega,
                     lambda = lambda,
                     gamma = gamma,
                     IC.score = IC.score)
    }

  } else {
    result <- est.obj[!(names(est.obj) %in% c("X", "S"))]
  }

  class(result) <- c("spice")
  return(result)
}

