#' Tuning Parameter Selection Information Criteria
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
#' @param crit A character string specifying the tuning parameter selection
#' criterion to use.
#'
#' @param ebic.tuning A numeric value (default = 0.5) specifying the tuning
#' parameter to calculate for \code{crit = "EBIC"}.
#'
#' @return
#' A numeric value.
#'
#' @noRd

criterion <- function(hatOmega, S, n, crit, ebic.tuning = 0.5) {
  ## dimensionality
  p <- ncol(S)
  ## Gaussian log-likelihood
  loglik <- (n/2) * (determinant(hatOmega, logarithm = TRUE)$modulus[1] - sum(diag(S%*%hatOmega)))
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

