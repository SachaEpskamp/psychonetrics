\name{ri_clpm}
\alias{ri_clpm}
\alias{ri_clpm_stationary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random intercept cross-lagged panel models
}
\description{
Function to run random intercept cross-lagged panel models under the lvm framework.
}
\usage{
ri_clpm(data, vars, lambda, 
        type = c("cov","chol","prec","ggm"), 
        verbose = FALSE, ...)
ri_clpm_stationary(x, 
                   stationary = c("intercepts", 
                                  "contemporaneous", 
                                  "innovation", 
                                  "temporal"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}

  \item{stationary}{
  The part of the model to implement stationarity on.
}
  
  \item{data}{
  A data frame encoding the data used in the analysis. Can be missing if covs and nobs are supplied.
}

  \item{vars}{
  Required argument. Different from in other psychonetrics models, this must be a *matrix* with each row indicating a variable and each column indicating a measurement. The matrix must be filled with names of the variables in the dataset corresponding to variable i at wave j. NAs can be used to indicate missing waves. The rownames of this matrix will be used as variable names.
  }
  
  \item{lambda}{
  A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
  }
  
  \item{type}{
  The type of model to model innovation
  }
  
  \item{verbose}{
  Logical, should progress be printed to the console?
  }
  
  \item{...}{
  Arguments sent to \code{lvm}
  }
  
}
\value{
A single psychonetrics object
}
\author{
Sacha Epskamp
}
\examples{
# to be made
}
