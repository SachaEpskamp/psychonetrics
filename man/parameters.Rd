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
parameters(x, standardized = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{standardized}{
Logical. If \code{TRUE}, the completely standardized (\code{std.all}) solution is
computed and two extra columns are added next to the unstandardized estimate:
\code{std} (the standardized point estimate) and \code{se_std} (its standard error,
obtained by the delta method). Standardization rescales every variable to unit
variance using the model-implied standard deviations, so that loadings, regression
coefficients and (co)variances are expressed on a standardized scale (covariances
become correlations; residual/disturbance variances become proportions of variance;
residual/disturbance covariances become residual correlations), matching
\code{lavaan::standardizedSolution(type = "std.all")} (Rosseel, 2012). Currently
implemented for the \code{lvm} and \code{varcov} model families only; for all other
families the \code{std}/\code{se_std} columns are returned as \code{NA} (with a
warning). Default \code{FALSE} (unstandardized estimates only, exactly as before).
}
}
\details{
The standardized standard errors are obtained post-fit by differentiating the
standardization map with respect to the free parameters
(\code{numDeriv::jacobian}) and applying the delta method to the parameter
covariance matrix; a fixed parameter (e.g. a marker loading) has \code{se_std =
NA} because it has no sampling variance.
}
\references{
Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
\emph{Journal of Statistical Software, 48}(2), 1--36.
\doi{10.18637/jss.v048.i02}
}
\author{
Sacha Epskamp
}
\value{
Invisibly returns a data frame containing information on all parameters (including
the \code{std} and \code{se_std} columns when \code{standardized = TRUE}).
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
mod <- ggm(ConsData, vars = vars, omega = "zero")

# Run model:
mod <- mod \%>\% runmodel

# Parameter estimates:
mod \%>\% parameters

# Standardized (std.all) solution with standard errors:
mod \%>\% parameters(standardized = TRUE)
}