\name{runmodel}
\alias{runmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run a psychonetrics model
}
\description{
This is the main function  used to run a psychonetrics model.
}
\usage{
runmodel(x, level = c("gradient", "fitfunction"), addfit = TRUE, addMIs = TRUE, 
addSEs = TRUE, addInformation = TRUE, log = TRUE, verbose = TRUE, 
optim.control = list(), maxtry = 5, analyticFisher = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{level}{
Level at which the model should be estimated. Defaults to \code{"gradient"} to indicate the analytic gradient should be used.
}
  \item{addfit}{
Logical, should fit measures be added?
}
  \item{addMIs}{
Logical, should modification indices be added?
}
  \item{addSEs}{
Logical, should standard errors be added?
}
  \item{addInformation}{
Logical, should the Fisher information be added?
}
  \item{log}{
Logical, should the log be updated?
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{optim.control}{
A list with options for \code{optimr}
}
  \item{maxtry}{
Integer indicating how often estimation should be attempted.
}
  \item{analyticFisher}{
Logical, should the analytic Fisher information be used? If \code{FALSE}, numeric information is used instead.
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}
