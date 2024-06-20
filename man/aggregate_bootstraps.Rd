\name{aggregate_bootstraps}
\alias{aggregate_bootstraps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aggregate Bootstrapped Models
}
\description{
Aggregates bootstrap results into a \code{psychonetrics_bootstrap} object
}
\usage{
aggregate_bootstraps(sample, bootstraps, remove_problematic = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample}{
The original \code{psychonetrics} object (not bootstrapped)
}
  \item{bootstraps}{
A list of bootstrapped \code{psychonetrics} objects (i.e., using \code{bootstrap = TRUE})
}
\item{remove_problematic}{Remove bootstraps that did not converge (sum of absolute gradient > 1)}
}
\details{
After running this function, the helper functions \code{parameters}, \code{fit}, and \code{CIplot} can be used to investigate bootstrap results.  
}
\value{
An object of the class \code{psychonetrics_bootstrap}
}
\author{
Sacha Epskamp
}

