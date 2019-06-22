\name{parameters}
\alias{parameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print parameter estimates
}
\description{
This function will print a list of parameters of the model
}
\usage{
parameters(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
}
\author{
Sacha Epskamp
}
\value{
Invisibly returns a data frame containing information on all parameters.
}
\examples{
# Load bfi data from psych package:
library("psych")
data(bfi)

# Also load dplyr for the pipe operator:
library("dplyr")

# Let's take the agreeableness items, and gender:
ConsData <- bfi \%>\% 
  select(A1:A5, gender) \%>\% 
  na.omit # Let's remove missingness (otherwise use Estimator = "FIML)

# Define variables:
vars <- names(ConsData)[1:5]

# Let's fit a full GGM:
mod <- ggm(ConsData, vars = vars, omega = "empty")

# Run model:
mod <- mod \%>\% runmodel

# Parameter estimates:
mod \%>\% parameters
}