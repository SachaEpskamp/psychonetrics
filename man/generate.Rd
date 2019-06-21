\name{generate}
\alias{generate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate data from a fitted psychonetrics object
}
\description{
This function will generate new data from the estimated mean and variance-covariance structure of a psychonetrics model.
}
\usage{
generate(x, n = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{n}{
Number of cases to sample per group.
}
}

\author{
Sacha Epskamp
}
\value{
A data frame with simulated data
}