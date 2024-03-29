% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ledoit_wolf_est.R
\name{ledoit_wolf_est}
\alias{ledoit_wolf_est}
\title{Ledoit-Wolf shrinkage estimator}
\usage{
ledoit_wolf_est(X, method = "lin", res = "cov")
}
\arguments{
\item{X}{A data matrix.}

\item{method}{A character string (default = "lin") specifying the method used in
shrinkage, includes: \enumerate{
\item Linear shrinkage: "lin" \insertCite{ledoit2004well}{fcstat}.
\item Non-linear shrinkage \insertCite{ledoit2015spectrum,ledoit2017numerical}{fcstat}:
\itemize{
\item "nlminb": non-linear shrinkage using the optimization routine "nlminb".
\item "nloptr": non-linear shrinkage using the optimization routine "nloptr".
}
}
See \code{\link[stats]{nlminb}} and \code{\link[nloptr]{nloptr}} for details.}

\item{res}{A character string (default = "cov") specifying the result matrix to be
obtained, either the covariance matrix ("cov") or the correlation matrix ("cor").}
}
\value{
A numeric matrix.
}
\description{
Compute the Ledoit-Wolf shrinkage estimator for the covariance or correlation matrix.
}
\references{
\insertAllCited{}
}
