\name{compare}
\alias{compare}
\alias{print.psychonetrics_compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model comparison
}
\description{
This function will print a table comparing multiple models on chi-square, AIC and BIC.
}
\usage{
compare(...)

\method{print}{psychonetrics_compare}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
Any number of \code{psychonetrics} models. Can be named to change the rownames of the output.
}
\item{x}{Output of the \code{compare} function.}
}
\value{
A data frame with chi-square values, degrees of freedoms, RMSEAs, AICs, and BICs.
}
\author{
Sacha Epskamp
}