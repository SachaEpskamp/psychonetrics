\name{transmod}
\alias{transmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform between model types
}
\description{
This function allows to transform a model variance--covariance structure from one type to another. Its main uses are to (1) use a Cholesky decomposition to estimate a saturated covariance matrix or GGM, and (2) to transform between conditional (ggm) and marginal associations (cov).
}
\usage{
transmod(x, ..., verbose, keep_computed = FALSE, log = TRUE,
         identify = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A psychonetrics model
}
  \item{\dots}{
Named arguments with the new types to use (e.g., \code{between = "ggm"} or \code{y = "cov"})
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{keep_computed}{
Logical, should the model be stated to be uncomputed adter the transformation? In general, a model does not need to be re-computed as transformed parameters should be at the maximum likelihood estimate.
}
  \item{log}{
Logical, should a logbook entry be made?
}
 \item{identify}{
Logical, should the model be identified after transforming?
}
}
\details{
Transformations are only possible if the model is diagonal (e.g., no partial correlations) or saturated (e.g., all covariances included).
}
\author{
Sacha Epskamp
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

# Model with Cholesky decompositon:
mod <- varcov(ConsData, vars = vars, type = "chol")

# Run model:
mod <- mod \%>\% runmodel

# Transform to GGM:
mod_trans <- transmod(mod, type = "ggm") \%>\% runmodel
# Note: runmodel often not needed

# Obtain thresholded GGM:
getmatrix(mod_trans, "omega", threshold = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
