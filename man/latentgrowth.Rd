\name{latentgrowth}
\alias{latentgrowth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Latnet growth curve model
}
\description{
Wrapper to \code{lvm} to specify a latent growth curve model.
}
\usage{
latentgrowth(vars, time = seq_len(ncol(vars)) - 1, covariates =
                   character(0), covariates_as = c("regression",
                   "covariance"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vars}{
Different from in other psychonetrics models, this must be a *matrix* with each row indicating a variable and each column indicating a measurement. The matrix must be filled with names of the variables in the dataset corresponding to variable i at wave j. NAs can be used to indicate missing waves. The rownames of this matrix will be used as variable names.
}
  \item{time}{
A vector with the encoding of each measurement (e.g., 0, 1, 2, 3).
}
  \item{covariates}{
A vector with strings indicating names of between-person covariate variables in the data
}
\item{covariates_as}{
Should covariates be included as regressions or actual covariates?
}
  \item{\dots}{
Arguments sent to \code{\link{lvm}}
}
}
\details{
See \url{https://github.com/SachaEpskamp/SEM-code-examples/tree/master/Latent_growth_examples/psychonetrics} for examples
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class}). See for an example \url{https://github.com/SachaEpskamp/SEM-code-examples/tree/master/Latent_growth_examples/psychonetrics}.
}
\author{
Sacha Epskamp
}
\examples{
library("dplyr")

# Smoke data cov matrix, based on LISS data panel https://www.dataarchive.lissdata.nl
smoke <- structure(c(47.2361758611759, 43.5366809116809, 41.0057465682466, 
                     43.5366809116809, 57.9789886039886, 47.6992521367521, 
                     41.0057465682466, 
                     47.6992521367521, 53.0669434731935), .Dim = c(3L, 3L), 
                   .Dimnames = list(
                       c("smoke2008", "smoke2009", "smoke2010"), c("smoke2008", 
                   "smoke2009", "smoke2010")))

# Design matrix:
design <- matrix(rownames(smoke),1,3)

# Form model:
mod <- latentgrowth(vars = design, 
                covs = smoke, nobs = 352
)

\dontrun{
# Run model:
mod <- mod \%>\% runmodel

# Evaluate fit:
mod \%>\% fit

# Look at parameters:
mod \%>\% parameters
}
}

