#' Performance Measures
#'
#' @description
#' Calculate various measures to evaluate the performance of the estimated
#' precision matrix.
#'
#' @param hatOmega The estimated precision matrix.
#'
#' @param Omega The reference precision matrix.
#'
#' @return A list containing the following components: \describe{
#' \item{Fnorm}{Frobenius (Hilbert-Schmidt) norm between the true and estimated
#' precision matrices.}
#' \item{KL}{Kullback-Leibler divergence between the true and estimated
#' precision matrices.}
#' \item{Snorm}{Spectral (operator) norm of the difference between the true and
#' estimated precision matrices.}
#' \item{precision}{Precision measure, the ratio of true positives to the total
#' predicted positives.}
#' \item{recall}{Recall measure, also known as Sensitivity, the ratio of true
#' positives to the total actual positives.}
#' \item{specificity}{Specificity measure, the ratio of true negatives to
#' the total actual negatives.}
#' \item{F1}{F1 score, the harmonic mean of Precision and Recall.}
#' \item{MCC}{Matthews correlation coefficient, a measure of the quality of
#' binary classifications.}
#' \item{sparsity}{The proportion of zeros among edges in the estimated
#' precision matrix.}
#' }
#'
#' @export

performance <- function(hatOmega, Omega) {
  Sigma <- solve(Omega)
  p <- ncol(hatOmega)
  Fnorm <- norm(Omega - hatOmega, "F")
  KL <- sum(diag(Sigma%*%hatOmega)) -
    determinant(Sigma%*%hatOmega, logarithm = TRUE)$modulus[1] - p
  Snorm <- svd(Omega - hatOmega)$d[1]
  Omega_edge <- Omega[upper.tri(Omega)]
  hatOmega_edge <- hatOmega[upper.tri(hatOmega)]
  FN <- sum(Omega_edge != 0 & hatOmega_edge == 0)
  FP <- sum(Omega_edge == 0 & hatOmega_edge != 0)
  TN <- sum(Omega_edge == 0 & hatOmega_edge == 0)
  TP <- sum(Omega_edge != 0 & hatOmega_edge != 0)
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  F1 <- 2*TP / (FN+FP+2*TP)
  MCC <- (TN*TP-FN*FP) / (sqrt(FN+TN)*sqrt(FN+TP)*sqrt(FP+TN)*sqrt(FP+TP))
  sparsity <- sum(hatOmega_edge == 0) /length(hatOmega_edge)
  result <- c(Fnorm = Fnorm, KL = KL, Snorm = Snorm,
              precision = precision, recall = recall, specificity = specificity,
              F1 = F1, MCC = MCC, sparsity = sparsity)
  return(result)
}

