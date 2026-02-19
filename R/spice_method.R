#' spice_method: "glasso", "ridge", "elnet", "clime", "tiger"
#'
#' @import RBGL
#' @import graph
#' @importFrom ADMMsigma ADMMsigma
#' @importFrom CovTools PreEst.glasso
#' @importFrom CVglasso CVglasso
#' @importFrom Glarmadillo glarma
#' @importFrom glasso glasso
#' @importFrom GLassoElnetFast gelnet
#' @importFrom glassoFast glassoFast
#' @importFrom huge huge.glasso
#' @importFrom ADMMsigma RIDGEsigma
#' @importFrom porridge ridgePgen
#' @importFrom rags2ridges ridgeP
#' @importFrom clime clime
#' @importFrom flare sugm
#' @importFrom huge huge.tiger
#' @importFrom Rdpack reprompt
#'
#' @noRd

spice_method <- function(method, X = NULL, S = NULL,
                         lambda = NULL, gamma = NULL,
                         pkg = NULL) {
  if (method == "glasso") {
    if (pkg == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = 1, diagonal = TRUE)$Z
    } else if (pkg == "CovTools") {
      hatOmega <- CovTools::PreEst.glasso(X = X, method = list(type = "fixed", param = lambda))$C
    } else if (pkg == "CVglasso") {
      hatOmega <- CVglasso::CVglasso(S = S, lam = lambda, diagonal = TRUE)$Omega
    } else if (pkg == "Glarmadillo") {
      hatOmega <- Glarmadillo::glarma(s = S, rho = lambda)$Theta
    } else if (pkg == "glasso") {
      hatOmega <- glasso::glasso(s = S, rho = lambda, penalize.diagonal = TRUE, start = "cold")$wi
    } else if (pkg == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 1, penalize.diagonal = TRUE)$Theta
    } else if (pkg == "glassoFast") {
      hatOmega <- glassoFast::glassoFast(S = S, rho = lambda, start = "cold")$wi
    } else if (pkg == "huge") {
      hatOmega <- huge::huge.glasso(x = S, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  } else if (method == "ridge") {
    if (pkg == "ADMMsigma") {
      hatOmega <- ADMMsigma::RIDGEsigma(S = S, lam = lambda)$Omega
    } else if (pkg == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = 0)$Theta
    } else if (pkg == "porridge") {
      hatOmega <- porridge::ridgePgen(S = S, lambda = matrix(lambda, ncol(S), ncol(S)), target = matrix(0, ncol(S), ncol(S)))
    } else if (pkg == "rags2ridges") {
      hatOmega <- rags2ridges::ridgeP(S = S, lambda = lambda, target = matrix(0, ncol(S), ncol(S)))
    }
  } else if (method == "elnet") {
    if (pkg == "ADMMsigma") {
      hatOmega <- ADMMsigma::ADMMsigma(S = S, lam = lambda, alpha = gamma)$Z
    } else if (pkg == "GLassoElnetFast") {
      hatOmega <- GLassoElnetFast::gelnet(S = S, lambda = lambda, alpha = gamma)$Theta
    }
  } else if (method == "clime") {
    if (pkg == "clime") {
      hatOmega <- clime::clime(x = S, lambda = lambda, sigma = TRUE, standardize = FALSE, linsolver = "simplex")$Omegalist[[1]]
    } else if (pkg == "flare") {
      hatOmega <- flare::sugm(data = S, lambda = lambda, method = "clime", verbose = FALSE)$icov[[1]]
    }
  } else if (method == "tiger") {
    if (pkg == "flare") {
      hatOmega <- flare::sugm(data = X, lambda = lambda, method = "tiger", verbose = FALSE)$icov[[1]]
    } else if (pkg == "huge") {
      hatOmega <- huge::huge.tiger(x = X, lambda = lambda, verbose = FALSE)$icov[[1]]
    }
  }
  return(hatOmega)
}

