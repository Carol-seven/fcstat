#' Tuning parameter selection information criteria
#'
#' @description
#' Provide information criteria for selecting tuning parameters.
#'
#' @param hatOmega The estimated precision matrix.
#'
#' @param S The sample covariance matrix.
#'
#' @param n An integer specifying the sample size.
#'
#' @param crit A string specifying the tuning parameter selection criterion to use.
#' Available options include: \enumerate{
#' \item "AIC": Akaike information criterion \insertCite{akaike1973information}{fcstat}.
#' \item "BIC": Bayesian information criterion \insertCite{schwarz1978estimating}{fcstat}.
#' \item "EBIC": extended Bayesian information criterion \insertCite{foygel2010extended}{fcstat}.
#' \item "HBIC": high dimensional Bayesian information criterion \insertCite{wang2013calibrating,fan2017high}{fcstat}.
#' }
#'
#' @param ebic.tuning A numeric value (default = 0.5) specifying the tuning parameter to
#' calculate when \code{crit = "EBIC"}.
#'
#' @importFrom Rdpack reprompt
#'
#' @return A numeric value.
#'
#' @references
#' \insertAllCited{}
#'
#' @export

criterion <- function(hatOmega, S, n, crit, ebic.tuning = 0.5) {
  ## dimensionality
  p <- ncol(S)
  ## Gaussian log-likelihood
  loglik <- (n/2) * (log(det(hatOmega)) - sum(diag(S%*%hatOmega)))
  ## cardinality of the edge set
  edges <- sum(hatOmega[upper.tri(hatOmega)] != 0)
  if (crit == "AIC") {
    result <- -2*loglik + 2*edges
  } else if (crit == "BIC") {
    result <- -2*loglik + log(n)*edges
  } else if (crit == "EBIC") {
    result <- -2*loglik + log(n)*edges + 4*ebic.tuning*log(p)*edges
  } else if (crit == "HBIC") {
    result <- -2*loglik + log(log(n))*log(p)*edges
  }
  return(result)
}

