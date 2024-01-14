#' Non-convex regularization penalty
#'
#' @description
#' Compute the specified non-convex regularization penalty.
#'
#' @param penalty A string specifying the non-convex penalty to use. Options include:
#' \enumerate{
#' \item "atan": arctangent type penalty \insertCite{wang2016variable}{fcstat}.
#' \item "exp": exponential type penalty \insertCite{wang2018variable}{fcstat}.
#' \item "mcp": minimax concave penalty \insertCite{zou2006adaptive}{fcstat}.
#' \item "scad": smoothly clipped absolute deviation \insertCite{fan2001variable,fan2009network}{fcstat}.
#' }
#'
#' @param Omega The precision matrix.
#'
#' @param lambda A numeric value representing the regularization parameter.
#'
#' @param gamma A numeric value representing the hyperparameter for the penalty function.
#' The defaults are: \enumerate{
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

pen <- function(penalty, Omega, lambda, gamma) {
  if (penalty == "atan") {
    if (missing(gamma)) {
      gamma <- 0.005
    }
    pen_mat <- lambda * (gamma + 2/pi) * atan(abs(Omega)/gamma)
  } else if (penalty == "exp") {
    if (missing(gamma)) {
      gamma <- 0.01
    }
    pen_mat <- lambda * (1 - exp(-abs(Omega)/gamma))
  } else if (penalty == "mcp") {
    if (missing(gamma)) {
      gamma <- 3
    }
    pen_mat <- (lambda*abs(Omega) - Omega^2/(2*gamma)) * (abs(Omega) <= gamma*lambda) +
      (gamma*lambda^2/2) * (abs(Omega) > gamma*lambda)
  } else if (penalty == "scad") {
    if (missing(gamma)) {
      gamma <- 3.7
    }
    pen_mat <- lambda*abs(Omega) * (abs(Omega) <= lambda) +
      (2*gamma*lambda*abs(Omega)-Omega^2-lambda^2)/(2*(gamma-1)) * (lambda < abs(Omega) & abs(Omega) <= gamma*lambda) +
      lambda^2*(gamma+1)/2 * (abs(Omega) > gamma*lambda)
  }
  return(pen_mat)
}

