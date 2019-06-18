\name{unionmodel}
\alias{unionmodel}
\alias{intersectionmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unify models across groups
}
\description{
The \code{unionmodel} will add all parameters to all groups that are free in at least one group, and the \code{intersectionmodel} will constrain all parameters across groups to zero unless they are free to estimate in all groups.
}
\usage{
unionmodel(x, runmodel = FALSE, verbose = TRUE, log = TRUE, identify = TRUE, ...)

intersectionmodel(x, runmodel = FALSE, verbose = TRUE, log = TRUE, identify = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A \code{psychonetrics} model.
}
  \item{runmodel}{
Logical, should the model be updated?
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{log}{
Logical, should the log be updated?
}
  \item{identify}{
Logical, should the model be identified?
}
 \item{\dots}{
Arguments sent to \code{runmodel}
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}
