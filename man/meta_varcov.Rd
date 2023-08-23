\name{meta_varcov}
\alias{meta_varcov}
\alias{meta_ggm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance-covariance and GGM meta analysis
}
\description{
Meta analysis of correlation matrices to fit a homogenous correlation matrix or Gaussian graphical model. Based on meta-analytic SEM (Jak and Cheung, 2019).
}
\usage{
meta_varcov(cors, nobs, Vmats, Vmethod = c("individual", "pooled",
                    "metaSEM_individual", "metaSEM_weighted"), Vestimation
                    = c("averaged", "per_study"), type = c("cor", "ggm"),
                    sigma_y = "full", kappa_y = "full", omega_y = "full",
                    lowertri_y = "full", delta_y = "full", rho_y = "full",
                    SD_y = "full", randomEffects = c("chol", "cov",
                    "prec", "ggm", "cor"), sigma_randomEffects = "full",
                    kappa_randomEffects = "full", omega_randomEffects =
                    "full", lowertri_randomEffects = "full",
                    delta_randomEffects = "full", rho_randomEffects =
                    "full", SD_randomEffects = "full", vars,
                    baseline_saturated = TRUE, optimizer, estimator =
                    c("FIML", "ML"), sampleStats, verbose = FALSE)
  
meta_ggm(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cors}{
A list of correlation matrices. Must contain rows and columns with \code{NA}s for variables not included in a study. 
}
  \item{nobs}{
A vector with the number of observations per study.
}
  \item{Vmats}{
Optional list with 'V' matrices (sampling error variance approximations).
}
  \item{Vmethod}{
Which method should be used to apprixomate the sampling error variance?
}
\item{Vestimation}{
How should the sampling error estimates be evaluated?
}
  \item{type}{
What to model? Currently only \code{"cor"} and \code{"ggm"} are supported.
}
  \item{sigma_y}{
Only used when \code{type = "cov"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_y}{
Only used when \code{type = "prec"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_y}{
Only used when \code{type = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_y}{
Only used when \code{type = "chol"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_y}{
Only used when \code{type = "ggm"}. Either \code{"diag"} or \code{"zero"} (not recommended), or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{rho_y}{
Only used when \code{type = "cor"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{SD_y}{
Only used when \code{type = "cor"}. Either \code{"diag"} or \code{"zero"}, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{randomEffects}{
What to model for the random effects? 
}
  \item{sigma_randomEffects}{
Only used when \code{type = "cov"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_randomEffects}{
Only used when \code{randomEffects = "prec"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_randomEffects}{
Only used when \code{randomEffects = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_randomEffects}{
Only used when \code{randomEffects = "chol"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_randomEffects}{
Only used when \code{randomEffects = "ggm"}. Either \code{"diag"} or \code{"zero"}, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{rho_randomEffects}{
Only used when \code{randomEffects = "cor"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{SD_randomEffects}{
Only used when \code{randomEffects = "cor"}. Either \code{"diag"} or \code{"zero"}, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{vars}{
Variables to be included.
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
  \item{optimizer}{
The optimizer to be used. Can be one of \code{"nlminb"} (the default R \code{nlminb} function), \code{"ucminf"} (from the \code{optimr} package), and C++ based optimizers \code{"cpp_L-BFGS-B"}, \code{"cpp_BFGS"}, \code{"cpp_CG"}, \code{"cpp_SANN"}, and \code{"cpp_Nelder-Mead"}. The C++ optimizers are faster but slightly less stable. Defaults to \code{"nlminb"}.
}
  \item{estimator}{
The estimator to be used. Currently implemented are \code{"ML"} for maximum likelihood estimation or \code{"FIML"} for full-information maximum likelihood estimation.
}
 \item{sampleStats}{
An optional sample statistics object. Mostly used internally. 
}
  \item{verbose}{
Logical, should progress be printed to the console?
}
\item{\dots}{
Arguments sent to \code{meta_varcov}
}
}

\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\references{
Jak, S., and Cheung, M. W. L. (2019). Meta-analytic structural equation modeling with moderating effects on SEM parameters. Psychological methods.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
