\name{dependency_share}
\alias{dependency_share}
\title{
Share of the dependence between composites carried by a composite model
}
\description{
Computes the share of the dependence between composites that is carried by a composite model made with \code{\link{cvm}} or \code{\link{cnm}}: \code{1 - (chisq - df) / (chisq_0 - df_0)}, in which \code{chisq} and \code{df} are the chi-square and degrees of freedom of the model, and \code{chisq_0} and \code{df_0} those of a baseline model. In the baseline model the indicators of different composites are uncorrelated, and all means, variances and covariances within composites are free. This function is new in version 0.18 and experimental (a message is shown once per session): arguments and defaults may still change.
}
\usage{
dependency_share(x, verbose = TRUE)
}
\arguments{
  \item{x}{
A \code{psychonetrics} model made with \code{\link{cvm}} or \code{\link{cnm}}, which has been run with \code{\link{runmodel}}. See Details for the models that are supported.
}
  \item{verbose}{
Logical, should the share be printed to the console?
}
}
\details{
The share is the relative noncentrality index (the comparative fit index without truncation) of the model, with the baseline model described above in place of the default baseline model. The estimate is not restricted to the range from 0 to 1: it is larger than 1 when the chi-square of the model is smaller than its degrees of freedom.

In the baseline model, all composites with a single indicator are combined into one set of indicators with free covariances, as a model with all relations between the composites reproduces the covariances among these indicators exactly. The baseline model is estimated with the same estimator and from the same sample statistics as the model.

The share is available for models with a single group that are estimated with the \code{"ML"} or \code{"FIML"} estimator, in which:
\itemize{
\item all constructs are composites (no \code{lambda}), and every variable is an indicator of one composite;
\item all composites are related: all edges in \code{omega_zeta} (\code{latent = "ggm"}) or all covariances in \code{sigma_zeta} (\code{latent = "cov"}) are free, and there are no effects in \code{beta};
\item all variances and covariances of the indicators within a composite are free (the matrices of the \code{composite} argument), as are the intercepts.
}
The weights can be free or fixed (e.g., \code{sumscores = TRUE}). The function stops with an error for other models, and when the chi-square of the baseline model is not larger than its degrees of freedom.
}
\value{
The share as a single number, returned invisibly.
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{cvm}}, \code{\link{cnm}}, \code{\link{fit}}
}
\examples{
library("dplyr")

# Load data:
data(NA2020)

# Two sleep problems as one composite, and six variables as composites with a
# single indicator (these represent the observed variables themselves):
vars <- c("regular_sleep", "worried_sleep", "sleep_interfere", "happy_health",
          "optimistic_future", "very_happy", "feel_alone", "happy_love_life")
composites <- c("regular_sleep", "sleep_problems", "happy_health",
                "optimistic_future", "very_happy", "feel_alone", "happy_love_life")
weights <- matrix(0, 8, 7)
weights[cbind(1:8, c(1, 2, 2, 3, 4, 5, 6, 7))] <- 1

# Composite network model:
mod <- cnm(NA2020, weights = weights, vars = vars, latents = composites,
           estimator = "FIML") \%>\% runmodel

# Share of the dependence between the composites carried by the model:
dependency_share(mod)

# The value is returned invisibly:
share <- dependency_share(mod, verbose = FALSE)
share
}
