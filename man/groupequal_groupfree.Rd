\name{groupequal}
\alias{groupequal}
\alias{groupfree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Group equality constrains
}
\description{
The \code{groupequal} function constrains parameters equal across groups, and the \code{groupfree} function frees equality constrains across groups.
}
\usage{
groupequal(x, matrix, row, col, verbose = TRUE, log = TRUE, runmodel = FALSE, 
      identify = TRUE, ...)

groupfree(x, matrix, row, col, verbose = TRUE, log = TRUE, runmodel = FALSE, 
      identify = TRUE, ...)
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
  \item{verbose}{
Logical, should messages be printed?
}
  \item{log}{
Logical, should the log be updated?
}
  \item{runmodel}{
Logical, should the model be updated?
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

