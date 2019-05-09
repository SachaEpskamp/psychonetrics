\name{bifactor}
\alias{bifactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bifactor helper function
}
\description{
This is a convenience function used to specify a bifactor model with \code{\link{lvm}}.
}
\usage{
bifactor(data, lambda, latents, bifactor = "g", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{lambda}{
A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.}
  \item{latents}{
An optional character vector with names of the latent variables.
}
  \item{bifactor}{
Name of the bifactor
}
  \item{\dots}{
Arguments sent to \code{\link{lvm}}
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}