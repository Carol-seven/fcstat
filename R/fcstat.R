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
#' @param method A character string specifying the statistical methods for estimating
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
#' \item "clime_primaldual": the utility function \code{\link[clime]{clime}} from
#' the package \code{clime} with the linsolver \code{primaldual}.
#' \item "clime_simplex": the utility function \code{\link[clime]{clime}} from
#' the package \code{clime} with the linsolver \code{simplex}.
#' \item "flare": the utility function \code{\link[flare]{sugm}} from
#' the package \code{flare}.
#' }
#'
#' @param n An integer (default - NULL) specifying the sample size. This is only required
#' when the input matrix \code{X} is a p-by-p sample covariance/correlation matrix with
#' dimension p.
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
#' @note
#' For the method \code{tiger}, the estimation process solely relies on the raw n-by-p
#' data \code{X} and does not utilize the argument \code{base}. This argument is not
#' applicable for \code{tiger} and will have no effect if provided.
#'
#' @importFrom Rdpack reprompt
#'
#' @return
#' \itemize{
#' \item For \code{crit = "CV"}, an object with S3 class "fcstat.sel" containing the
#' following components: \describe{
#' \item{hatOmega_opt}{The estimated precision matrix.}
#' \item{lambda_opt}{The optimal regularization parameter.}
#' \item{gamma_opt}{The optimal hyperparameter.}
#' \item{loss_opt}{The optimal k-fold loss.}
#' \item{lambda}{The actual lambda grid used in the program.}
#' \item{gamma}{The actual gamma grid used in the program.}
#' \item{loss.mean}{The mean of k-fold loss for each parameter grid value.}
#' \item{loss.sd}{The standard deviation of k-fold loss for each parameter grid value.}
#' }
#' \item For other criteria, an object with S3 class "fcstat.sel" containing the following
#' components: \describe{
#' \item{hatOmega_opt}{The estimated precision matrix.}
#' \item{lambda_opt}{The optimal regularization parameter.}
#' \item{gamma_opt}{The optimal hyperparameter.}
#' \item{score_opt}{The optimal information criterion score.}
#' \item{lambda}{The actual lambda grid used in the program.}
#' \item{gamma}{The actual gamma grid used in the program.}
#' \item{score}{The information criterion score for each parameter grid value.}
#' }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export

fcstat <- function(
    X, method, base = "cov",
    lambda = NULL, nlambda = 50, lambda.min.ratio = NULL,
    gamma = NA, ## for elnet, adapt, atan, exp, mcp, scad
    target = 0, ## for ridge, elnet
    initial = "linshrink", ## initial estimator for atan, exp, mcp, scad; adaptive weight for adapt
    utilopt = "flare", ## utility option for clime
    n = NULL, crit = "CV", fold = 5, ebic.tuning = 0.5) {

  est.obj <- fcstat.est(
    X = X, method = method, base = base,
    lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
    gamma = gamma, target = target, initial = initial, utilopt = utilopt)

  result <- fcstat.sel(est.obj = est.obj, n = n, crit = crit, fold = fold, ebic.tuning = ebic.tuning)
  class(result) <- "fcstat"

  return(result)
}

