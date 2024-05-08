#' Performance measures
#'
#' @description
#' Calculate various measures to evaluate the performance of the estimated precision matrix.
#'
#' @param Omega The precision matrix.
#'
#' @param Sigma The covariance matrix (default = NULL).
#'
#' @param hatOmega The estimated precision matrix.
#'
#' @return A list containing the following components: \describe{
#' \item{Fnorm1}{Frobenius (Hilbert-Schmidt) norm between the true and estimated precision
#' matrices.}
#' \item{Fnorm2}{Frobenius (Hilbert-Schmidt) norm between the product of the true
#' covariance matrix and estimated precision matrix and the identity matrix.}
#' \item{KL}{Kullback-Leibler divergence between the true and estimated precision matrices.}
#' \item{Ql}{Quadratic loss between the diagonal elements of the product of the true
#' covariance matrix and estimated precision matrix.}
#' \item{Snorm}{Spectral (operator) norm of the difference between the true and estimated
#' precision matrices.}
#' \item{Precision}{Precision measure, the ratio of true positives to the total predicted
#' positives.}
#' \item{Recall}{Recall measure, also known as Sensitivity, the ratio of true positives to
#' the total actual positives.}
#' \item{Specificity}{Specificity measure, the ratio of true negatives to the total actual
#' negatives.}
#' \item{F1}{F1 score, the harmonic mean of Precision and Recall.}
#' \item{MCC}{Matthews correlation coefficient, a measure of the quality of binary
#' classifications.}
#' }
#'
#' @export

performance <- function(Omega, Sigma = NULL, hatOmega) {
  if (missing(Omega)) {
    if (is.null(Sigma)) {
      stop("Either 'Sigma' or 'Omega' must be provided.")
    } else {
      Omega <- solve(Sigma)
    }
  } else {
    if (is.null(Sigma)) {
      Sigma <- solve(Omega)
    }
  }

  p <- ncol(hatOmega)
  Fnorm1 <- norm(Omega - hatOmega, "F")
  Fnorm2 <- norm(Sigma%*%hatOmega - diag(p), "F")
  KL <- sum(diag(Sigma%*%hatOmega)) - log(det(Sigma%*%hatOmega)) - p
  Ql <- sum(diag(Sigma%*%hatOmega - diag(p))^2)
  Snorm <- svd(Omega - hatOmega)$d[1]
  # temp <- (Omega - hatOmega) %*% (Omega - hatOmega)
  # e1 <- max(Re(eigen(temp)$values))
  # Snorm <- sqrt(e1)
  FN <- sum(Omega != 0 & hatOmega == 0)
  FP <- sum(Omega == 0 & hatOmega != 0)
  TN <- sum(Omega == 0 & hatOmega == 0)
  TP <- sum(Omega != 0 & hatOmega != 0)
  Precision <- TP / (TP + FP)
  Recall <- TP / (TP + FN)
  Specificity <- TN / (TN + FP)
  F1 <- 2*TP / (FN+FP+2*TP)
  # F1 <- 2*Precision*Recall / (Precision+Recall)
  MCC <- (TN*TP-FN*FP) / (sqrt(FN+TN)*sqrt(FN+TP)*sqrt(FP+TN)*sqrt(FP+TP))
  result <- c(Fnorm1 = Fnorm1, Fnorm2 = Fnorm2,
              KL = KL, Ql = Ql, Snorm = Snorm,
              Precision = Precision, "Recall/Sensitivity" = Recall, Specificity = Specificity,
              F1 = F1, MCC = MCC)
  return(result)
}

