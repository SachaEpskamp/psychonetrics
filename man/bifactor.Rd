\name{bifactor}
\alias{bifactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bi-factor models
}
\description{
Wrapper to \code{lvm} to specify a bi-factor model.
}
\usage{
bifactor(data, lambda, latents, bifactor = "g", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data as used by \code{\link{lvm}}
}
  \item{lambda}{
The factor loadings matrix *without* the bifactor, as used by by \code{\link{lvm}}
}
  \item{latents}{
A vector of names of the latent variables, as used by \code{\link{lvm}}
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
