\name{ml_tsdlvm1}
\alias{ml_tsdlvm1}
\alias{ml_ts_lvgvar}
\alias{ml_gvar}
\alias{ml_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-level Lag-1 dynamic latent variable model family of psychonetrics models for time-series data
}
\description{
This function is a wrapper around \code{\link{dlvm1}} that allows for specifying the model using a long format data and similar input as the \code{mlVAR} package. The \code{ml_ts_lvgvar} simply sets \code{within_latent = "ggm"} and \code{between_latent = "ggm"} by default. The \code{ml_gvar} and \code{ml_var} are simple wrappers with different named defaults for contemporaneous and between-person effects.
}
\usage{
ml_tsdlvm1(data, beepvar, idvar, vars, groups, estimator = "FIML", 
  standardize = c("none", "z", "quantile"), ...)

ml_ts_lvgvar(...)

ml_gvar(..., contemporaneous = c("ggm", "cov", "chol", "prec"), 
        between = c("ggm", "cov", "chol", "prec"))
             
ml_var(..., contemporaneous = c("cov", "chol", "prec", "ggm"), 
        between = c("cov", "chol", "prec", "ggm"))
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

  \item{contemporaneous}{
The type of within-person latent contemporaneous model to be used.
}
  \item{between}{
The type of between-person latent model to be used.
}
  \item{\dots}{
Arguments sent to \code{\link{dlvm1}}
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
