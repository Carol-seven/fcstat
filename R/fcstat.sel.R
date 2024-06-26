#' Selection for precision matrix estimates
#'
#' @description
#' Provide a collection of methods to perform parameter selection for \code{fcstat.est} object.
#'
#' @param est.obj The \code{fcstat.est} object output from \code{fcstat.est()}.
#'
#' @param n An integer (default - NULL) specifying the sample size. This is only required
#' when the input matrix \code{X} in \code{fcstat.est} is a p-by-p sample
#' covariance/correlation matrix with dimension p.
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
#' @import foreach
#' @importFrom stats cov sd
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
#' @autoglobal
#'
#' @export

fcstat.sel <- function(est.obj, n = NULL, crit = "CV", fold = 5, ebic.tuning = 0.5) {

  list2env(est.obj, envir = environment())

  if (is.null(X)) {
    if (is.null(n)) {
      stop("The input X is the p-by-p sample covariance/correlation matrix.
           The selection requires the sample size n.")
    }
  } else {
    n <- nrow(X)
  }
  n.para <- length(hatOmega)

  if (crit == "CV") {

    if(is.null(X)) {
      stop("CV requires the n-by-p data matrix!")
    }

    ## indices for test and training sets
    n.test <- floor(n/fold)
    n.train <- n - n.test
    mat.test <- matrix(NA, n.test, fold)
    mat.train <- matrix(NA, n.train, fold)
    index <- sample(n)
    for (j in 1:fold) {
      ind.test <- ((j-1)*n.test+1):(j*n.test)
      mat.test[,j] <- index[ind.test]
      ind.train <- (1:n)[!(1:n) %in% ind.test]
      mat.train[,j] <- index[ind.train]
    }

    ## initialize the loss
    loss <- matrix(NA, fold, n.para)

    ## loop over each fold
    for (j in 1:fold) {

      ## training and test sets
      X.train <- X[mat.train[,j],]
      X.test <- X[mat.test[,j],]

      ## sample covariance/correlation matrix
      S.test <- eval(parse(text =  paste0(base, "(X.test)")))

      ## compute the precision matrix estimate
      cvlist <- fcstat.est(X = X.train,
                           method = method, base = base,
                           lambda = lambda, gamma = gamma,
                           target = target, initial = initial)

      ## loss: negative log-likelihood
      for (k in 1:n.para) {
        loss[j,k] <- - log(det(cvlist$hatOmega[[k]])) + sum(diag(S.test%*%cvlist$hatOmega[[k]]))
        # determinant(cvlist$hatOmega[[k]])$modulus[1]
      }
    }

    loss[!is.finite(loss)] <- Inf

    ## the mean and sd of the k-fold loss for each parameter grid value
    loss.mean <- apply(loss, 2, mean)
    loss.sd <- apply(loss, 2, sd)

    ## find the optimal parameter
    index <- which.min(loss.mean)

    result <- list(hatOmega_opt = hatOmega[[index]],
                   lambda_opt = lambda[index],
                   gamma_opt = gamma[index],
                   loss_opt = loss.mean[index],
                   lambda = lambda,
                   gamma = gamma,
                   loss.mean = loss.mean,
                   loss.sd = loss.sd)

  } else {

    ## select the optimal parameters among a set of possible values
    eval(parse(text = paste0(
      "score <- foreach(k = 1:n.para, .combine = 'c') %do% {",
      "criterion(hatOmega = hatOmega[[k]], S = S, n = n, crit = '", crit,
      "', ebic.tuning = ", ebic.tuning, ")}"
    )))

    ## find the optimal parameter
    index <- which.min(score)

    result <- list(hatOmega_opt = hatOmega[[index]],
                   lambda_opt = lambda[index],
                   gamma_opt = gamma[index],
                   score_opt = score[index],
                   lambda = lambda,
                   gamma = gamma,
                   score = score)

  }

  class(result) <- c("fcstat.sel")
  return(result)

}

