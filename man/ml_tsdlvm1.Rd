\name{ml_tsdlvm1}
\alias{ml_tsdlvm1}
\alias{ml_ts_lvgvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-level Lag-1 dynamic latent variable model family of psychonetrics models for time-series data
}
\description{
This function is a wrapper around \code{\link{dlvm1}} that allows for specifying the model using a long format data and similar input as the \code{mlVAR} package.
}
\usage{
ml_tsdlvm1(data, beepvar, idvar, vars, groups, estimator = "FIML", 
  standardize = c("none", "z", "quantile"), ...)

ml_ts_lvgvar(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data to be used. Must be raw data in long format (each row indicates one person at one time point).
}
  \item{beepvar}{
Optional string indicating assessment beep per day. Adding this argument will cause non-consecutive beeps to be treated as missing!
}
  \item{idvar}{
	
String indicating the subject ID
}
  \item{vars}{
Vectors of variables to include in the analysis
}
  \item{groups}{
An optional string indicating the name of the group variable in \code{data}.
}
  \item{estimator}{
Estimator to be used. Must be \code{"FIML"}.
}
  \item{standardize}{
Which standardization method should be used? \code{"none"} (default) for no standardization, \code{"z"} for z-scores, and \code{"quantile"} for a non-parametric transformation to the quantiles of the marginal standard normal distribution.
}
  \item{\dots}{
Arguments sent to \code{\link{dlvm1}}
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
