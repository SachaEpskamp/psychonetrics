\name{latentgrowth}
\alias{latentgrowth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Latent growth helper function
}
\description{
This is a convenience function used to specify a latent growth model with \code{\link{lvm}}.
}
\usage{
latentgrowth(vars, time = seq_len(ncol(vars)), covariates = character(0), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vars}{
An optional character vector encoding the variables used in the analyis. Must equal names of the dataset in \code{data}.
}
  \item{time}{
Vector encoding the time points
}
  \item{covariates}{
An optional character vector encoding covariates used in the analyis. Must equal names of the dataset in \code{data}.
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
