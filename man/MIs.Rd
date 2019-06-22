\name{MIs}
\alias{MIs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print modification indices
}
\description{
This function prints a list of modification indices (MIs)
}
\usage{
MIs(x, all = FALSE, matrices, type = c("normal", "equal", "free"), top = 10, 
    verbose = TRUE, nonZero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{all}{
Logical, should all MIs be printed or only the highest?
}
  \item{matrices}{
Optional vector of matrices to include in the output.
}
  \item{type}{
String indicating which kind of modification index should be printed. (\code{"mi"} is the typical MI, \code{"mi_free"} is the modification index free from equality constrains across groups, and \code{"mi_equal"} is the modification index if the parameter is added constrained equal across all groups).
}
  \item{top}{
Number of MIs to include in output if \code{all = FALSE}
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{nonZero}{
Logical, should only MIs be printed of non-zero parameters? Useful to explore violations of group equality.
}
}
\author{
Sacha Epskamp
}
\value{
Invisibly returns a relevant subset of the data frame containing all information on the parameters, or a list of such data frames if multiple types of MIs are requested.
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

# Modification indices:
mod \%>\% MIs
}