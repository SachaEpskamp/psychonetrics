\name{modelsearch}
\alias{modelsearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stepwise model search
}
\description{
This function peforms stepwise model search to find an optimal model that (locally) minimzes some criterion (by default, the BIC). 
}
\usage{
modelsearch(x, criterion = "bic", matrices, prunealpha = 0.001, addalpha = 0.05, 
  verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{criterion}{
String indicating the criterion to minimize. Any criterion from \code{\link{fit}} can be used.
}
  \item{matrices}{
Vector of strings indicating which matrices should be searched. Will default to network structures and factor loadings.
}
  \item{prunealpha}{
Minimal alpha used to consider edges to be removed
}
  \item{addalpha}{
Maximum alpha used to consider edges to be added
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{\dots}{
Arguments sent to \code{\link{runmodel}}
}
}
\details{
The full algorithm is as follows:

1. Evaluate all models in which an edge is removed that has p > prunealpha, or an edge is added that has a modification index with p < addalpha

2. If none of these models improve the criterion, return the previous model and stop the algorithm

3. Update the model to the model that improved the criterion the most

4. Evaluate all other considered models that improved the criterion 

5. If none of these models improve the criterion, go to 1, else go to 3

}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}

\author{
Sacha Epskamp
}



\seealso{
\code{\link{prune}}, \code{\link{stepup}}
}
\examples{
\donttest{
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

# Let's fit a full GGM:
mod <- ggm(ConsData, vars = vars)

# Run model:
mod <- mod \%>\% runmodel

# Model search
mod <- mod \%>\% prune(alpha= 0.01) \%>\% modelsearch
}
}
