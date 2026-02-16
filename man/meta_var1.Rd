\name{meta_var1}
\alias{meta_var1}
\alias{meta_gvar}
\title{
Meta-analytic VAR(1) model
}
\description{
Single-stage meta-analytic vector autoregressive model (VAR(1)) for time-series data collected across multiple studies. Pools temporal (lag-1) and contemporaneous network structures across studies using a random-effects framework. The \code{meta_gvar} wrapper sets \code{contemporaneous = "ggm"} for graphical VAR estimation.
}
\usage{
meta_var1(data, covs, nobs,
          vars, studyvar, idvar, dayvar, beepvar,
          contemporaneous = c("cov", "chol", "prec", "ggm"),
          beta = "full",
          omega_zeta = "full", delta_zeta = "full",
          kappa_zeta = "full", sigma_zeta = "full",
          lowertri_zeta = "full",
          randomEffects = c("chol", "cov", "prec", "ggm", "cor"),
          sigma_randomEffects = "full",
          kappa_randomEffects = "full",
          omega_randomEffects = "full",
          lowertri_randomEffects = "full",
          delta_randomEffects = "full",
          rho_randomEffects = "full",
          SD_randomEffects = "full",
          Vmats,
          Vmethod = c("individual", "pooled"),
          Vestimation = c("averaged", "per_study"),
          baseline_saturated = TRUE, optimizer,
          estimator = c("FIML", "ML"),
          sampleStats, verbose = FALSE,
          bootstrap = FALSE, boot_sub, boot_resample)

meta_gvar(...)
}
\arguments{
  \item{data}{
A data frame containing time-series data for multiple studies. Use together with \code{studyvar} to identify studies.
}
  \item{covs}{
A list of pre-computed Toeplitz covariance matrices (one per study). Each matrix should have dimension \code{2p x 2p} where \code{p} is the number of variables, with the first \code{p} rows/columns being lagged variables and the second \code{p} being current variables. Alternative to \code{data}.
}
  \item{nobs}{
A vector with the number of observations per study. Required when \code{covs} is supplied.
}
  \item{vars}{
Character vector of observed variable names. If missing, inferred from data.
}
  \item{studyvar}{
A string indicating the column name in \code{data} that identifies the study. Required when \code{data} is supplied. If not supplied but \code{idvar} is, \code{idvar} is used as \code{studyvar} with a warning.
}
  \item{idvar}{
Optional string indicating the subject ID column within each study. When both \code{studyvar} and \code{idvar} are supplied, data is first split by \code{studyvar}, then \code{idvar} is used within each study to collate time series (as in \code{var1}).
}
  \item{dayvar}{
Optional string indicating the day variable (passed to \code{tsData} per study).
}
  \item{beepvar}{
Optional string indicating the beep variable within day (passed to \code{tsData} per study).
}
  \item{contemporaneous}{
Parameterization of the contemporaneous (residual innovation) covariance structure. One of \code{"cov"}, \code{"chol"}, \code{"prec"}, or \code{"ggm"}.
}
  \item{beta}{
Temporal lag-1 regression matrix specification. Defaults to \code{"full"} (all elements free, including diagonal autoregressive effects).
}
  \item{omega_zeta}{
Contemporaneous partial correlation matrix specification (used when \code{contemporaneous = "ggm"}).
}
  \item{delta_zeta}{
Contemporaneous scaling matrix specification (used when \code{contemporaneous = "ggm"}).
}
  \item{kappa_zeta}{
Contemporaneous precision matrix specification (used when \code{contemporaneous = "prec"}).
}
  \item{sigma_zeta}{
Contemporaneous covariance matrix specification (used when \code{contemporaneous = "cov"}).
}
  \item{lowertri_zeta}{
Contemporaneous Cholesky factor specification (used when \code{contemporaneous = "chol"}).
}
  \item{randomEffects}{
Parameterization of the random effects covariance structure.
}
  \item{sigma_randomEffects}{
Random effects covariance matrix specification (used when \code{randomEffects = "cov"}).
}
  \item{kappa_randomEffects}{
Random effects precision matrix specification (used when \code{randomEffects = "prec"}).
}
  \item{omega_randomEffects}{
Random effects partial correlation matrix specification (used when \code{randomEffects = "ggm"}).
}
  \item{lowertri_randomEffects}{
Random effects Cholesky factor specification (used when \code{randomEffects = "chol"}).
}
  \item{delta_randomEffects}{
Random effects scaling matrix specification (used when \code{randomEffects = "ggm"}).
}
  \item{rho_randomEffects}{
Random effects correlation matrix specification (used when \code{randomEffects = "cor"}).
}
  \item{SD_randomEffects}{
Random effects standard deviation matrix specification (used when \code{randomEffects = "cor"}).
}
  \item{Vmats}{
Optional list with 'V' matrices (sampling error variance approximations).
}
  \item{Vmethod}{
Which method should be used to approximate the sampling error variance? \code{"individual"} or \code{"pooled"}.
}
  \item{Vestimation}{
How should the sampling error estimates be evaluated? \code{"averaged"} or \code{"per_study"}.
}
  \item{baseline_saturated}{
Logical indicating if baseline and saturated models should be included.
}
  \item{optimizer}{
The optimizer to be used. Defaults to \code{"nlminb"}.
}
  \item{estimator}{
The estimator to be used. \code{"ML"} or \code{"FIML"} (default).
}
  \item{sampleStats}{
Optional sample statistics object.
}
  \item{verbose}{
Logical, should progress be printed?
}
  \item{bootstrap}{
Should the data be bootstrapped?
}
  \item{boot_sub}{
Proportion of cases to subsample for bootstrap.
}
  \item{boot_resample}{
Logical, should bootstrap be with replacement?
}
  \item{\dots}{
Arguments sent to \code{meta_var1}.
}
}
\details{
This function implements a single-stage meta-analytic VAR(1) model. For \code{p} observed variables, each study's time-series data is transformed into a Toeplitz covariance structure from which two blocks are extracted:

\itemize{
  \item \bold{Sigma_0}: the \code{p x p} stationary covariance (symmetric, \code{p(p+1)/2} unique elements)
  \item \bold{Sigma_1}: the \code{p x p} lag-1 cross-covariance (non-symmetric, \code{p^2} elements)
}

The VAR(1) structural model implies:
\deqn{\textrm{vec}(\Sigma_0) = (I - \beta \otimes \beta)^{-1} \textrm{vec}(\Sigma_\zeta)}
\deqn{\Sigma_1 = \beta \Sigma_0}

The meta-level mean is \eqn{\mu = [\textrm{vech}(\Sigma_0), \textrm{vec}(\Sigma_1)]} and heterogeneity is modeled as \eqn{\Sigma = \Sigma^{(\textrm{ran})} + V}.

Note: the exogenous (lagged) covariance block from the full Toeplitz matrix is not modeled, as it is a nuisance parameter at the meta-analytic level.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\examples{
\dontrun{
library(mlVAR)
set.seed(42)
Model <- mlVARsim(nPerson = 50, nNode = 3, nTime = 50, lag = 1)

# Fit meta-analytic VAR(1):
mod <- meta_var1(Model$Data,
                 vars = Model$varNames,
                 studyvar = Model$idvar,
                 estimator = "ML")
mod <- mod \%>\% runmodel

# Inspect results:
print(mod)
getmatrix(mod, "beta")

# Fit meta-analytic GVAR (with GGM contemporaneous):
mod_gvar <- meta_gvar(Model$Data,
                      vars = Model$varNames,
                      studyvar = Model$idvar,
                      estimator = "ML")
mod_gvar <- mod_gvar \%>\% runmodel
getmatrix(mod_gvar, "omega_zeta")
}
}
\seealso{
\code{\link{var1}}, \code{\link{meta_varcov}}, \code{\link{meta_lvm}}
}
