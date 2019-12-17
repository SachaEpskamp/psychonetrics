\name{fit}
\alias{fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print fit indices
}
\description{
This function will print all fit indices of the model/
}
\usage{
fit(x)
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
Invisibly returns a data frame with fit measure estimates.
}
\examples{
# Load bfi data from psych package:
library("psychTools")
data(bfi)

# Also load dplyr for the pipe operator:
library("dplyr")

# Let's take the agreeableness items, and gender:
ConsData <- bfi \%>\% 
  select(A1:A5, gender) \%>\% 
  na.omit # Let's remove missingness (otherwise use Estimator = "FIML)

# Define variables:
vars <- names(ConsData)[1:5]

# Let's fit an empty GGM:
mod0 <- ggm(ConsData, vars = vars, omega = "empty")

# Run model:
mod0 <- mod0 \%>\% runmodel

# Inspect fit:
mod0 \%>\% fit # Pretty bad fit...
}