#' Derivative of the non-convex regularization penalty
#'
#' @description
#' Compute the derivative of the specified non-convex regularization penalty.
#'
#' @param penalty A character string specifying the non-convex penalty to use.
#' Available options include: \enumerate{
#' \item "adapt": adaptive lasso \insertCite{zou2006adaptive,fan2009network}{fcstat}.
#' \item "atan": arctangent type penalty \insertCite{wang2016variable}{fcstat}.
#' \item "exp": exponential type penalty \insertCite{wang2018variable}{fcstat}.
#' \item "mcp": minimax concave penalty \insertCite{zou2006adaptive}{fcstat}.
#' \item "scad": smoothly clipped absolute deviation \insertCite{fan2001variable,fan2009network}{fcstat}.
#' }
#'
#' @param Omega The precision matrix.
#'
#' @param lambda A scalar specifying the regularization parameter.
#'
#' @param gamma A scalar specifying the hyperparameter for the penalty function.
#' The defaults are: \enumerate{
#' \item "adapt": 0.5
#' \item "atan": 0.005
#' \item "exp": 0.01
#' \item "mcp": 3
#' \item "scad": 3.7
#' }
#'
#' @importFrom Rdpack reprompt
#'
#' @return A numeric matrix.
#'
#' @references
#' \insertAllCited{}
#'
#' @export

deriv <- function(penalty, Omega, lambda, gamma) {
  if (penalty == "adapt") {
    if (missing(gamma)) {
      gamma <- 0.5
    }
    weight <- abs(Omega)^(-gamma)
    weight[is.infinite(weight)] <- 9999999999
    lambda_mat <- lambda*weight
  } else if (penalty == "atan") {
    if (missing(gamma)) {
      gamma <- 0.005
    }
    lambda_mat <- lambda * gamma * (gamma + 2/pi) / (gamma^2 + Omega^2)
  } else if (penalty == "exp") {
    if (missing(gamma)) {
      gamma <- 0.01
    }
    lambda_mat <- (lambda/gamma) * exp(-abs(Omega)/gamma)
  } else if (penalty == "mcp") {
    if (missing(gamma)) {
      gamma <- 3
    }
    lambda_mat <- (lambda - abs(Omega)/gamma) * (abs(Omega) <= gamma*lambda)
  } else if (penalty == "scad") {
    if (missing(gamma)) {
      gamma <- 3.7
    }
    lambda_mat <- lambda * (abs(Omega) <= lambda) +
      pmax(gamma*lambda - abs(Omega), 0) / (gamma-1) * (abs(Omega) > lambda)
  }
  return(lambda_mat)
}

