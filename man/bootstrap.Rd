\name{bootstrap}
\alias{bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap a psychonetrics model
}
\description{
This function will bootstrap the data (once) and return a new unevaluated psychonetrics object. It requres \code{storedata = TRUE} to be used when forming a model.
}
\usage{
bootstrap(x, replacement = TRUE, proportion = 1, verbose = TRUE, storedata = FALSE, 
          baseline_saturated = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{replacement}{
Logical, should new samples be drawn with replacement?
}
  \item{proportion}{
Proportion of sample to be drawn. Set to lower than $1$ for subsampling.
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{storedata}{
Logical, should the bootstrapped data also be stored?
}
  \item{baseline_saturated}{
Logical, should the baseline and saturated models be included?
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}