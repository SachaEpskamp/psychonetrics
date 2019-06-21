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
addMIs(x, matrices = "all", type = c("normal", "free", "equal"), verbose = TRUE, 
        analyticFisher = TRUE)

addSEs(x)

addfit(x)

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
}
\author{
Sacha Epskamp
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}