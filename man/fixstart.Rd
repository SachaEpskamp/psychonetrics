\name{fixstart}
\alias{fixstart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Attempt to Fix Starting Values
}
\description{
This function attempts to fix starting values by comparing the analytic gradient to a numerically approximated gradient. Parameters with a difference between the analytic and numeric gradient that exceeds 'maxdiff' will be reduced by a factor of 'reduce' in each iteration until the average absolute difference between analytic and numeric gradients is lower than 'tol'. Only off-diagonal elements in omega, sigma, kappa, lowertri or rho matrices or any element in beta matrices are adjusted.
}
\usage{
fixstart(x, reduce = 0.5, maxdiff = 0.1, tol = 0.01, maxtry = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A 'psychonetrics' model
}
  \item{reduce}{
The factor with which problematic parameters are reduced in each iteration.
}
  \item{maxdiff}{
Maximum difference between analytic and numeric gradient to be considered problematic.
}
  \item{tol}{
Average absolute difference between analytic and numeric gradient that is considered acceptable.
}
  \item{maxtry}{
Maximum number of iterations to attempt to fix starting values.
}
}

\author{
Sacha Epskamp
}
