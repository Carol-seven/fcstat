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
#' @param crit A string (default = "CV") specifying the parameter selection method to use.
#' Available options include: \enumerate{
#' \item "AIC": Akaike information criterion \insertCite{akaike1973information}{fcstat}.
#' \item "BIC": Bayesian information criterion \insertCite{schwarz1978estimating}{fcstat}.
#' \item "EBIC": extended Bayesian information criterion \insertCite{foygel2010extended}{fcstat}.
#' \item "HBIC": high dimensional Bayesian information criterion \insertCite{wang2013calibrating,fan2017high}{fcstat}.
#' \item "CV": k-fold cross validation.
#' }
#'
#' @param fold An integer (default = 5) specifying the number of folds used for \code{crit = "CV"}.
#'
#' @param ebic.tuning A scalar (default = 0.5) specifying the tuning parameter to
#' calculate when \code{crit = "EBIC"}.
#'
#' @note
#' For the method \code{tiger}, the estimation process solely relies on the raw n-by-p
#' data \code{X} and does not utilize the arguments \code{base} and \code{approach}.
#' These arguments are not applicable for \code{tiger} and will have no effect if provided.
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
    X, method,
    base = "cov", approach = "smp",
    lambda = NULL, nlambda = 50, lambda.min.ratio = NULL,
    lambda.min = NULL, lambda.max = NULL, ## for clime
    gamma = NA, ## for elnet, adapt (with adapt.weight = "Fan"), atan, exp, scad, mcp
    target = 0, ## for ridge, elnet
    initial = "glasso", ## initial estimator for atan, exp, scad, mcp; adaptive weight for adapt
    crit = "CV", fold = 5, ebic.tuning = 0.5) {

  est.obj <- fcstat.est(
    X = X, method = method, base = base, approach = approach,
    lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
    lambda.min = lambda.min, lambda.max = lambda.max,
    gamma = gamma, target = target, initial = initial)

  result <- fcstat.sel(est.obj = est.obj, crit = crit, fold = fold, ebic.tuning = ebic.tuning)
  class(result) <- "fcstat"

  return(result)
}

