\name{getVCOV}
\alias{getVCOV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the asymptotic covariance matrix
}
\description{
This function can be used to obtain the estimated asymptotic covariance matrix from a \code{psychonetrics} object.
}
\usage{
getVCOV(model, approximate_SEs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A \code{psychonetrics} model.
}
\item{approximate_SEs}{
Logical, should standard errors be approximated? If true, an approximate matrix inverse of the Fischer information is used to obtain the standard errors.
}
}

\value{
This function returns a matrix.
}
\author{
Sacha Epskamp
}
