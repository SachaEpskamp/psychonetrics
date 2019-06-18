\name{parequal}
\alias{parequal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set equality constrains across parameters
}
\description{
This function can be used to set parameters equal
}
\usage{
parequal(x, ..., inds = integer(0), verbose = TRUE, log = TRUE, runmodel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
 \item{\dots}{
Arguments sent to \code{runmodel}
}
  \item{inds}{
Parameter indices of parameters to be constrained equal
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{log}{
Logical, should the log be updated?
}
  \item{runmodel}{
Logical, should the model be updated?
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}
