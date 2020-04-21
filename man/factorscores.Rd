\name{factorscores}
\alias{factorscores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute factor scores
}
\description{
Currently, only the \code{lvm} framework with single group and no missing data is supported.
}
\usage{
factorscores(data, model, method = c("bartlett", "regression"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Dataset to compute factor scores for
}
  \item{model}{
A psychonetrics model
}
  \item{method}{
The method to use: \code{"regression"} or \code{"bartlett"}
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
