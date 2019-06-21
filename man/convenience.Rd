\name{setestimator}
\alias{setestimator}
\alias{setoptimizer}
\alias{usecpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convenience functions
}
\description{
These functions can be used to change some estimator options.
}
\usage{
setestimator(x, estimator)

setoptimizer(x, optimizer)

usecpp(x, use = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{estimator}{
A string indicating the estimator to be used
}
  \item{optimizer}{
A string indicating the optimizer to be used
}
\item{use}{
Logical indicating if C++ should be used (currently only used in FIML)
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}

