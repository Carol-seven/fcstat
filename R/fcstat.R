#' Statistical methods for estimating precision matrix
#'
#' @description
#' Provide a collection of statistical methods to estimate a precision matrix.
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
#' @param n An integer (default = NULL) specifying the sample size. This is only required
#' when the input matrix \code{X} is a p-by-p sample covariance/correlation matrix with
#' dimension p.
#'
#' @param lambda Grid of non-negative scalars for the regularization parameter.
#' The default is \code{NULL}, which generates its own \code{lambda} sequence based on
#' \code{nlambda} and \code{lambda.min.ratio}.
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
#' Some options are also offered when a character string is provided (default "linshrink"),
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
#' \item For \code{method} set to \code{"adapt"}, \code{"atan"}, \code{"exp"},
#' \code{"scad"}, and \code{"mcp"}: \itemize{
#' \item "glasso_cold": the function from \code{\link[glasso]{glasso}}.
#' \item "glasso_warm": the function from \code{\link[glasso]{glasso}} with warm-starts.
#' \item "GLassoElnetFast_cold": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet}.
#' \item "GLassoElnetFast_warm": the function from
#' \href{https://github.com/TobiasRuckstuhl/GLassoElnetFast}{gelnet} with warm-starts.
#' \item "glassoFast_cold": the function from \code{\link[glassoFast]{glassoFast}}.
#' \item "glassoFast_warm": the function from \code{\link[glassoFast]{glassoFast}} with
#' warm-starts.
#' }
#' }
#'
#' @param crit A string (default = "CV") specifying the parameter selection method to use.
#' Available options include: \enumerate{
#' \item "AIC": Akaike information criterion \insertCite{akaike1973information}{fcstat}.
#' \item "BIC": Bayesian information criterion \insertCite{schwarz1978estimating}{fcstat}.
#' \item "EBIC": extended Bayesian information criterion \insertCite{foygel2010extended}{fcstat}.
#' \item "HBIC": high dimensional Bayesian information criterion \insertCite{wang2013calibrating,fan2017high}{fcstat}.
#' \item "CV": k-fold cross validation with negative log-likelihood loss.
#' }
#'
#' @param fold An integer (default = 5) specifying the number of folds used for \code{crit = "CV"}.
#'
#' @param ebic.tuning A scalar (default = 0.5) specifying the tuning parameter to
#' calculate for \code{crit = "EBIC"}.
#'
#' @param cores An integer (default = 1) specifying the number of cores to use for
#' parallel execution.
#'
#' @note
#' For the method \code{tiger}, the estimation process solely relies on the raw n-by-p
#' data \code{X} and does not utilize the argument \code{base}. This argument is not
#' applicable for \code{tiger} and will have no effect if provided.
#'
#' @importFrom stats cor cov sd
#' @importFrom Rdpack reprompt
#'
#' @return
#' \itemize{
#' \item For \code{crit = "CV"}, an object with S3 class "fcstat" containing the following
#' components: \describe{
#' \item{hatOmega_opt}{The estimated precision matrix.}
#' \item{lambda_opt}{The optimal regularization parameter.}
#' \item{gamma_opt}{The optimal hyperparameter.}
#' \item{loss_opt}{The optimal k-fold loss.}
#' \item{hatOmega}{A list of estimated precision matrices for \code{lambda} grid and
#' \code{gamma} grid.}
#' \item{lambda}{The actual lambda grid used in the program.}
#' \item{gamma}{The actual gamma grid used in the program.}
#' \item{loss.mean}{The mean of k-fold loss for each parameter grid value.}
#' \item{loss.sd}{The standard deviation of k-fold loss for each parameter grid value.}
#' }
#' \item For other criteria, an object with S3 class "fcstat" containing the following
#' components: \describe{
#' \item{hatOmega_opt}{The estimated precision matrix.}
#' \item{lambda_opt}{The optimal regularization parameter.}
#' \item{gamma_opt}{The optimal hyperparameter.}
#' \item{score_opt}{The optimal information criterion score.}
#' \item{hatOmega}{A list of estimated precision matrices for \code{lambda} grid and
#' \code{gamma} grid.}
#' \item{lambda}{The actual lambda grid used in the program.}
#' \item{gamma}{The actual gamma grid used in the program.}
#' \item{score}{The information criterion score for each parameter grid value.}
#' }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @autoglobal
#'
#' @export

fcstat <- function(
    X, method, base = "cov", n = NULL,
    lambda = NULL, nlambda = 20, lambda.min.ratio = 0.01,
    gamma = NA, ## for elnet, adapt, atan, exp, mcp, scad
    initial = "glasso", ## initial estimator for atan, exp, mcp, scad; adaptive weight for adapt
    pkgopt = "glassoFast_cold_diagture", ## package option
    crit = "CV", fold = 5, ebic.tuning = 0.5,
    cores = 1) {

  est.obj <- fcstat.est(
    X = X, method = method, base = base,
    lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
    gamma = gamma, initial = initial, pkgopt = pkgopt,
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
      loss <- matrix(NA, fold, npara)

      for (j in 1:fold) {

        ## indices for test and training sets; training and test sets
        ind.test <- index[((j-1)*floor(n/fold)+1):(j*floor(n/fold))]
        X.test <- X[ind.test, , drop = FALSE]
        ind.train <- index[index != ind.test]
        X.train <- X[ind.train, , drop = FALSE]

        ## sample covariance/correlation matrix
        if (base == "cov") {
          S.test <- cov(X.test)
        } else if (base == "cor") {
          S.test <- cor(X.test)
        }

        ## compute the precision matrix estimate
        cvlist <- fcstat.est(X = X.train,
                             method = method, base = base,
                             lambda = lambda, gamma = gamma,
                             initial = initial, pkgopt = pkgopt,
                             cores = cores)

        ## loss: negative log-likelihood
        for (k in 1:npara) {
          loss[j,k] <- - determinant(cvlist$hatOmega[[k]], logarithm = TRUE)$modulus[1] + sum(diag(S.test%*%cvlist$hatOmega[[k]]))
          # log(det(cvlist$hatOmega[[k]])) determinant(cvlist$hatOmega[[k]], logarithm = TRUE)$modulus[1]
        }
      }

      loss[!is.finite(loss)] <- Inf

      ## the mean and sd of the k-fold loss for each parameter grid value
      loss.mean <- apply(loss, 2, mean)
      loss.sd <- apply(loss, 2, sd)

      ## find the optimal parameter
      index <- which.min(loss.mean)

      result <- list(hatOmega_opt = hatOmega[[index]],
                     lambda_opt = lambda[index],
                     gamma_opt = gamma[index],
                     loss_opt = loss.mean[index],
                     hatOmega = hatOmega,
                     lambda = lambda,
                     gamma = gamma,
                     loss.mean = loss.mean,
                     loss.sd = loss.sd)

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
      score <- sapply(1:npara, function (k) {
        criterion(hatOmega = hatOmega[[k]], S = S, n = n, crit = crit, ebic.tuning = ebic.tuning)
      })

      ## find the optimal parameter
      index <- which.min(score)

      result <- list(hatOmega_opt = hatOmega[[index]],
                     lambda_opt = lambda[index],
                     gamma_opt = gamma[index],
                     score_opt = score[index],
                     hatOmega = hatOmega,
                     lambda = lambda,
                     gamma = gamma,
                     score = score)
    }

  } else {
    result <- est.obj[!(names(est.obj) %in% c("X", "S"))]
  }

  class(result) <- c("fcstat")
  return(result)
}

