#' spice_internal
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom glasso glasso
#' @importFrom GLassoElnetFast gelnet
#' @importFrom glassoFast glassoFast
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom Rdpack reprompt
#'
#' @noRd

spice_internal <- function(
    X, method, base = "cov",
    lambda = NULL, nlambda = 20, lambda.min.ratio = 0.01,
    gamma = NULL, ## for elnet, adapt, atan, exp, mcp, scad
    initial = "glasso", ## initial estimator for atan, exp, mcp, scad; adaptive weight for adapt
    pkg = "glasso", ## package option
    cores = 1) {

  if (!method %in% c("glasso", "ridge", "elnet", "clime", "tiger",
                     "adapt", "atan", "exp", "mcp", "scad")) {
    stop('Error in `method`!
         Available options: "glasso", "ridge", "elnet", "clime", "tiger", "adapt", "atan", "exp", "mcp", "scad".')
  }

  ## dimensionality
  p <- ncol(X)
  n <- nrow(X)

  ## sample covariance/correlation matrix
  if (isSymmetric(X)) {
    S <- X
    X <- NULL
  } else {
    if (base == "cov") {
      S <- (n-1)/n*cov(X)
    } else if (base == "cor") {
      S <- cor(X)
    }
  }

  ## lambda grid
  if(is.null(lambda)) {
    S_adjusted <- S - diag(p)
    lambda.max <- max(max(S_adjusted), -min(S_adjusted))
    lambda.min <- lambda.min.ratio*lambda.max
    lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }

  ## gamma grid
  if (all(is.na(gamma))) {
    gamma <- switch(method, "elnet" = seq(0.1, 0.9, 0.1), "adapt" = 0.5,
                    "atan" = 0.005, "exp"  = 0.01,
                    "mcp"  = 3, "scad" = 3.7, NA)
  }

  ## parameter grid combination
  parameter <- expand.grid(lambda = unique(lambda), gamma = unique(gamma))

  npara <- nrow(parameter)

  if (npara > 1 & cores > 1) {

    ## CPU cores
    num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
    if (cores > num.cores) {
      cat("The number of available CPU cores is ", num.cores, "!\n", sep = "")
    }
    if (cores > npara) {
      cores <- npara
    }
    cluster <- makeCluster(cores)
    registerDoParallel(cluster)

    ## compute the precision matrix estimator hatOmega along the parameter grid
    if (method %in% c("glasso", "ridge", "elnet", "clime", "tiger")) {
      hatOmega <- foreach(k = 1:npara, .packages = "spice",
                          .export = c("spice_method")) %dopar% {
        spice_method(method = method, X = X, S = S,
                     lambda = parameter$lambda[k], gamma = parameter$gamma[k],
                     pkg = pkg)
      }
    } else { ## method %in% c("adapt", "atan", "exp", "mcp", "scad")
      Omega <- gen_initial(X = X, S = S, base = base, initial = initial,
                           lambda = parameter$lambda, pkg = pkg)
      lambda_mat <- lapply(1:npara, function(k) {
        spice::ncv_deriv(penalty = method, Omega = Omega[[k]],
                         lambda = parameter$lambda[k], gamma = parameter$gamma[k])
      })
      if (pkg == "glasso") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          glasso::glasso(s = S, rho = lambda_mat[[k]], penalize.diagonal = TRUE, start = "cold")$wi
        }
      } else if (pkg == "GLassoElnetFast") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          GLassoElnetFast::gelnet(S = S, lambda = lambda_mat[[k]], alpha = 1, penalize.diagonal = TRUE)$Theta
        }
      } else if (pkg == "glassoFast") {
        hatOmega <- foreach(k = 1:npara) %dopar% {
          glassoFast::glassoFast(S = S, rho = lambda_mat[[k]], start = "cold")$wi
        }
      }
    }

    stopCluster(cluster)

  } else {

    ## compute the precision matrix estimator hatOmega along the parameter grid
    if (method %in% c("glasso", "ridge", "elnet", "clime", "tiger")) {
      hatOmega <- lapply(1:npara, function(k) {
        spice_method(method = method, X = X, S = S,
                     lambda = parameter$lambda[k], gamma = parameter$gamma[k],
                     pkg = pkg)
      })
    } else { ## method %in% c("adapt", "atan", "exp", "mcp", "scad")
      Omega <- gen_initial(X = X, S = S, base = base, initial = initial,
                           lambda = parameter$lambda, pkg = pkg)
      lambda_mat <- lapply(1:npara, function(k) {
        spice::ncv_deriv(penalty = method, Omega = Omega[[k]],
                         lambda = parameter$lambda[k], gamma = parameter$gamma[k])
      })
      if (pkg == "glasso") {
        hatOmega <- lapply(1:npara, function(k) {
          glasso::glasso(s = S, rho = lambda_mat[[k]], penalize.diagonal = TRUE, start = "cold")$wi
        })
      } else if (pkg == "GLassoElnetFast") {
        hatOmega <- lapply(1:npara, function(k) {
          GLassoElnetFast::gelnet(S = S, lambda = lambda_mat[[k]], alpha = 1, penalize.diagonal = TRUE)$Theta
        })
      } else if (pkg == "glassoFast") {
        hatOmega <- lapply(1:npara, function(k) {
          glassoFast::glassoFast(S = S, rho = lambda_mat[[k]], start = "cold")$wi
        })
      }
    }
  }

  result <- list(hatOmega = hatOmega,
                 lambda = parameter$lambda,
                 gamma = parameter$gamma,
                 X = X,
                 S = S)
  class(result) <- c("spice.internal")
  return(result)

}

