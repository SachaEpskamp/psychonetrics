\name{changedata}
\alias{changedata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Change the data of a psychonetrics object
}
\description{
This function can be used to change the data in a psychonetrics object.
}
\usage{
changedata(x, data, covs, nobs, means, groups, missing = "listwise")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{covs}{
A sample variance--covariance matrix, or a list/array of such matrices for multiple groups. IMPORTANT NOTE: psychonetrics expects the maximum likelihood (ML) covariance matrix, which is NOT obtained from \code{\link{cov}} directly. Manually rescale the result of \code{\link{cov}} with \code{(nobs - 1)/nobs} to obtain the ML covariance matrix.
}
  \item{nobs}{
The number of observations used in \code{covs} and \code{means}, or a vector of such numbers of observations for multiple groups.
}
  \item{means}{
A vector of sample means, or a list/matrix containing such vectors for multiple groups.
}
  \item{groups}{
An optional string indicating the name of the group variable in \code{data}.
}
  \item{missing}{
How should missingness be handled in computing the sample covariances and number of observations when \code{data} is used. Can be \code{"listwise"} for listwise deletion, or \code{"pairwise"} for pairwise deletion.
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}