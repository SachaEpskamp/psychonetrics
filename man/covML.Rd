\name{covML}
\alias{covML}
\alias{covMLtoUB}
\alias{covUBtoML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum likelihood covariance estimate
}
\description{
These functions complement the base R \code{cov} function by simplifying obtaining maximum likelihood (ML) covariance estimates (denominator n) instead of unbiased (UB) covariance estimates (denominator n-1). The function \code{covML} can be used to obtain ML estimates, the function \code{covUBtoML} transforms from UB to ML estimates, and the function \code{covMLtoUB} transforms from UB to ML estimates.
}
\usage{
covML(x, ...)
covUBtoML(x, n, ...)
covMLtoUB(x, n, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A dataset
}
\item{n}{
The sample size
}
  \item{\dots}{
Arguments sent to the \code{cov} function.
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{
data("StarWars")
Y <- StarWars[,1:10]

# Unbiased estimate:
UB <- cov(Y)

# ML Estimate:
ML <- covML(Y)

# Check:
all(abs(UB - covMLtoUB(ML, nrow(Y))) < sqrt(.Machine$double.eps))
all(abs(ML - covUBtoML(UB, nrow(Y))) < sqrt(.Machine$double.eps))
}
