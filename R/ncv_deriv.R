#' Derivative of a Non-Convex Penalty
#'
#' @description
#' Compute the derivative of a specified non-convex regularization penalty.
#'
#' @param penalty A character string specifying the non-convex penalty to use.
#' Available options include: \enumerate{
#' \item "adapt": Adaptive lasso
#' \insertCite{zou2006adaptive,fan2009network}{spice}.
#' \item "atan": Arctangent type penalty \insertCite{wang2016variable}{spice}.
#' \item "exp": Exponential type penalty \insertCite{wang2018variable}{spice}.
#' \item "mcp": Minimax concave penalty \insertCite{zou2006adaptive}{spice}.
#' \item "scad": Smoothly clipped absolute deviation
#' \insertCite{fan2001variable,fan2009network}{spice}.
#' }
#'
#' @param Omega The precision matrix.
#'
#' @param lambda A non-negative numeric value specifying the regularization
#' parameter.
#'
#' @param gamma A numeric value specifying the additional parameter for
#' the penalty function. The defaults are:
#' \enumerate{
#' \item "adapt": 0.5
#' \item "atan": 0.005
#' \item "exp": 0.01
#' \item "mcp": 3
#' \item "scad": 3.7
#' }
#'
#' @return
#' A numeric matrix.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom Rdpack reprompt
#'
#' @export

ncv_deriv <- function(penalty, Omega, lambda, gamma) {
  if (missing(gamma)) {
    gamma <- switch(penalty, "adapt" = 0.5, "atan" = 0.005,
                    "exp"  = 0.01, "mcp"  = 3, "scad" = 3.7, NA)
  }
  if (penalty == "adapt") {
    weight <- abs(Omega)^(-gamma)
    weight[is.infinite(weight)] <- 9999999999
    lambda_mat <- lambda*weight
  } else if (penalty == "atan") {
    lambda_mat <- lambda * gamma * (gamma + 2/pi) / (gamma^2 + Omega^2)
  } else if (penalty == "exp") {
    lambda_mat <- (lambda/gamma) * exp(-abs(Omega)/gamma)
  } else if (penalty == "mcp") {
    lambda_mat <- (lambda - abs(Omega)/gamma) * (abs(Omega) <= gamma*lambda)
  } else if (penalty == "scad") {
    lambda_mat <- lambda * (abs(Omega) <= lambda) +
      pmax(gamma*lambda - abs(Omega), 0) / (gamma-1) * (abs(Omega) > lambda)
  }
  return(lambda_mat)
}

