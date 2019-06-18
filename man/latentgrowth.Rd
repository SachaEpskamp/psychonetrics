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
latentgrowth(vars, time = seq_len(ncol(vars)) - 1, covariates = character(0), ...)
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

