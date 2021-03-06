% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DPtree.R
\name{pDPTreeRealize}
\alias{pDPTreeRealize}
\title{The disitribution function for realized distribution from D-P tree.}
\usage{
pDPTreeRealize(d, x)
}
\arguments{
\item{d}{A \eqn{2^m} by \eqn{2^m} matrix, m being the approximating level. Normalized measures for all \eqn{2^m} by \eqn{2^m} sub-partititons on copula space
given by the realized distribution from D-P tree, as returned by \code{DPTreeDensity}.}

\item{x}{An array of dimension n by 2. The points on copula space for distribution function evluation. Should be between 0 and 1.}
}
\value{
An array of length n. The values of CDF of the input D-P tree distribution evaluated at the input points.
}
\description{
\code{pDPTreeRealize} returns
the value of distribution function of realized distribution from D-P tree at certain given point on copula space.
}
\examples{
pDPTreeRealize(DPTreePMeanDensity(DPTreePrior(m=2, z=1)),c(0.5,0.5))
}
\references{
\insertRef{DPtree}{DPtree}
}
