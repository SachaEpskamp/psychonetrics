\name{stepup}
\alias{stepup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stepup model search along modification indices
}
\description{
This function automatically peforms step-up search by adding the parameter with the largest modification index until some criterion is reached or no modification indices are significant at alpha.
}
\usage{
stepup(x, alpha = 0.01, criterion = "bic", matrices, mi = c("mi", "mi_free", "mi_equal"), 
    greedyadjust = c("bonferroni", "none", "holm", "hochberg", "hommel", "fdr", "BH", 
    "BY"), stopif, greedy = FALSE, verbose = TRUE, checkinformation = TRUE, 
    singularinformation = c("tryfix", "skip", "continue", "stop"), startEPC = TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{alpha}{
Significance level to use.
}
  \item{criterion}{
String indicating the criterion to minimize. Any criterion from \code{\link{fit}} can be used.
}
  \item{matrices}{
Vector of strings indicating which matrices should be searched. Will default to network structures and factor loadings.
}
  \item{mi}{
String indicating which kind of modification index should be used (\code{"mi"} is the typical MI, \code{"mi_free"} is the modification index free from equality constrains across groups, and \code{"mi_equal"} is the modification index if the parameter is added constrained equal across all groups).
}
  \item{greedyadjust}{
String indicating which p-value adjustment should be used in greedy start. Any method from \code{p.adjust} can be used.
}
  \item{stopif}{
An R expression, using objects from \code{fit}, which will break stepup search if it evaluates to \code{TRUE}. For example, \code{stopif = rmsea < 0.05} will lead to search to stop if rmsea is below 0.05. 
}
  \item{greedy}{
Logical, should a greedy start be used? If \code{TRUE}, the first step adds any parameter that is significant (after adjustement)
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{checkinformation}{
Logical, should the Fisher information be checked for potentially non-identified models?
}
  \item{singularinformation}{
String indicating how to proceed if the information matrix is singular. \code{"tryfix"} 
will adjust starting values to try to fix the proble, \code{"skip"} will lead to the algorithm
to skip the current parameter, \code{"continue"} will ignore the situation, and \code{"stop"}
will break the algorithm and return a list with the last two models.
}
\item{startEPC}{
Logical, should the starting value be set at the expected parameter change?
}
  \item{\dots}{
Arguments sent to \code{\link{runmodel}}
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{prune}}
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
mod <- ggm(ConsData, vars = vars, omega = "full")

# Run model:
mod <- mod \%>\%runmodel \%>\%prune(alpha = 0.05)

# Remove an edge (example):
mod <- mod \%>\%fixpar("omega",1,2) \%>\%runmodel

# Stepup search
mod <- mod \%>\%stepup(alpha = 0.05)
}
}
