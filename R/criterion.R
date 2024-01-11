#' Tuning parameter selection criteria
#'
#' Provide criteria for selecting tuning parameters.
#'
#' @param hatOmega The estimated precision matrix.
#'
#' @param S The sample covariance matrix.
#'
#' @param n An integer specifying the sample size.
#'
#' @param crit A string specifying the tuning parameter selection criterion to use.
#' Available options include: \enumerate{
#' \item "AIC": Akaike information criterion (Akaike, 1998)
#' \item "BIC": Bayesian information criterion (Schwarz, 1978)
#' \item "EBIC": extended Bayesian information criterion (Foygel and Drton, 2010)
#' \item "HBIC": high dimensional Bayesian information criterion (Wang et al., 2013; Fan et al., 2017)
#' }
#'
#' @param ebic.tuning A numeric value (default = 0.5) specifying the tuning parameter to
#' calculate when \code{crit = "EBIC"}.
#'
#' @return A numeric value.
#'
#' @references \itemize{
#' \item Akaike, Hirotogu. (1973).
#' Information Theory and an Extension of the Maximum Likelihood Principle.
#' In B. N. Petrov and F. Csáki (Eds.),
#' \emph{Second International Symposium on Information Theory} (pp. 267-281).
#' Budapest, Hungary: Akadémiai Kiadó.
#' \item Fan, Jianqing and Liu, Han and Ning, Yang and Zou, Hui. (2017)
#' High Dimensional Semiparametric Latent Graphical Model for Mixed Data.
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology},
#' 79(2), 405--421.
#' \item Foygel, Rina and Drton, Mathias. (2010).
#' Extended Bayesian Information Criteria for Gaussian Graphical Models.
#' In J. Lafferty, C. Williams, J. Shawe-Taylor, R. Zemel, and A. Culotta (Eds.),
#' \emph{Advances in Neural Information Processing Systems 23 (NIPS 2010)}.
#' NY, USA: Curran Associates, Inc.
#' \item Schwarz, Gideon. (1978).
#' Estimating the Dimension of a Model.
#' \emph{The Annals of Statistics}, 6(2), 461--464.
#' \item Wang, Lan and Kim, Yongdai and Li, Runze. (2013).
#' Calibrating Nonconvex Penalized Regression in Ultra-High Dimension.
#' \emph{The Annals of Statistics}, 41(5), 2505--2536.
#' }
#' @export

criterion <- function(hatOmega, S, n, crit, ebic.tuning = 0.5) {
  ## Gaussian log-likelihood
  loglik <- (n/2) * (log(det(hatOmega)) - sum(diag(S%*%hatOmega)))
  ## cardinality of the edge set
  edges <- sum(hatOmega[upper.tri(hatOmega)] != 0)
  if (crit == "AIC") {
    result <- -2*loglik + 2*edges
  } else if (crit == "BIC") {
    result <- -2*loglik + log(n)*edges
  } else if (crit == "EBIC") {
    result <- -2*loglik + log(n)*edges + 4*ebic.tuning*log(ncol(S))*edges
  } else if (crit == "HBIC") {
    result <- -2*loglik + log(log(n))*log(ncol(S))*edges
  }
  return(result)
}

