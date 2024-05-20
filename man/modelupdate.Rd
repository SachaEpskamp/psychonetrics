\name{psychonetrics_update}
\alias{addMIs}
\alias{addSEs}
\alias{addfit}
\alias{identify}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model updating functions
}
\description{
These functions update a psychonetrics model. Typically they are not required.
}
\usage{
addMIs(x, matrices = "all", type = c("normal", "free",
                    "equal"), verbose, analyticFisher = TRUE)

addSEs(x, verbose, approximate_SEs = FALSE)

addfit(x, verbose)

identify(x)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{matrices}{
Optional vector of matrices to include in MIs.
}
  \item{type}{
String indicating which modification indices should be updated.
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{analyticFisher}{
Logical indicating if an analytic Fisher information matrix should be used.
}
\item{approximate_SEs}{
Logical, should standard errors be approximated? If true, an approximate matrix inverse of the Fischer information is used to obtain the standard errors.
}
}
\author{
Sacha Epskamp
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}