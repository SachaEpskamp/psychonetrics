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
runmodel(x, level = c("gradient", "fitfunction"), addfit =
                   TRUE, addMIs = TRUE, addSEs = TRUE, addInformation =
                   TRUE, log = TRUE, verbose, optim.control,
                   analyticFisher = TRUE, warn_improper = FALSE,
                   warn_gradient = TRUE, warn_bounds = TRUE,
                   return_improper = TRUE, bounded = TRUE,
                   approximate_SEs = FALSE,
                   criterion = "ebic.5", beta_min)

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
  \item{analyticFisher}{
Logical, should the analytic Fisher information be used? If \code{FALSE}, numeric information is used instead.
}
\item{return_improper}{
Should a result in which improper computation was used be return? Improper computation can mean that a pseudoinverse of small spectral shift was used in computing the inverse of a matrix.
}
\item{warn_improper}{Logical. Should a warning be given when at some point in the estimation a pseudoinverse was used?
}
\item{warn_gradient}{Logical. Should a warning be given when the average absolute gradient is > 1?}
\item{bounded}{Logical. Should bounded estimation be used (e.g., variances should be positive)?
}
\item{approximate_SEs}{
Logical, should standard errors be approximated? If true, an approximate matrix inverse of the Fischer information is used to obtain the standard errors.
}
\item{warn_bounds}{
Should a warning be given when a parameter is estimated near its bounds?
}
\item{criterion}{
Character string specifying the information criterion for automatic lambda selection in PML/PFIML models. One of \code{"bic"}, \code{"ebic.25"}, \code{"ebic.5"} (default), \code{"ebic.75"}, or \code{"ebic1"}. Only used when the model contains auto-select penalty parameters (\code{penalty_lambda = NA}). See \code{\link{find_penalized_lambda}} for details.
}
\item{beta_min}{
Numeric threshold for zeroing out small penalized parameters during automatic lambda selection. Default is \code{sqrt(log(p)/n)}. See \code{\link{find_penalized_lambda}} for details.
}
}
\details{
For penalized models (PML or PFIML) with auto-select penalty parameters (\code{penalty_lambda = NA}, the default), \code{runmodel} automatically calls \code{\link{find_penalized_lambda}} to select the optimal penalty strength via EBIC grid search before fitting. The \code{criterion} and \code{beta_min} arguments control this search. After the lambda search, use \code{\link{refit}} for post-selection inference with standard errors.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{find_penalized_lambda}} for manual control over the lambda search,
\code{\link{penalize}} for setting penalty parameters,
\code{\link{refit}} for post-selection inference after penalized estimation.
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

# Let's fit a full GGM:
mod <- ggm(ConsData, vars = vars, omega = "full")

# Run model:
mod <- mod \%>\% runmodel
}