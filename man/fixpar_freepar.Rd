\name{fixpar}
\alias{fixpar}
\alias{freepar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameters modification
}
\description{
The \code{fixpar} function can be used to fix a parameter to some value (Typically zero), and the \code{freepar} function can be used to free a parameter from being fixed to a value.
}
\usage{
fixpar(x, matrix, row, col, value = 0, group, verbose, 
        log = TRUE, runmodel = FALSE, ...)

freepar(x, matrix, row, col, start, group, verbose, log =
        TRUE, runmodel = FALSE, startEPC = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{matrix}{
String indicating the matrix of the parameter
}
  \item{row}{
Integer or string indicating the row of the matrix of the parameter
}
  \item{col}{
Integer or string indicating the column of the matrix of the parameter
}
  \item{value}{
Used in \code{fixpar} to indicate the value to which a parameters is constrained
}
  \item{start}{
Used in \code{freepar} to indicate the starting value of the parameter
}
  \item{group}{
Integer indicating the group of the parameter to be constrained
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
\item{startEPC}{
Logical, should the starting value be set at the expected parameter change?
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

