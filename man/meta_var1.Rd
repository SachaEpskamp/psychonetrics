\name{meta_var1}
\alias{meta_var1}
\alias{meta_var}
\alias{meta_gvar}
\title{
Meta-analytic VAR(1) model
}
\description{
Single-stage meta-analytic vector autoregressive model (VAR(1)) for time-series data collected across multiple studies. Pools temporal (lag-1) and contemporaneous network structures across studies using a random-effects framework. The \code{meta_gvar} wrapper sets \code{contemporaneous = "ggm"} for graphical VAR estimation, and \code{meta_var()} is the primary short name of \code{meta_var1()} (identical function).
}
\usage{
meta_var1(data, covs, nobs,
          vars, studyvar, idvar, dayvar, beepvar,
          contemporaneous = c("cov", "chol", "prec", "ggm", "cor"),
          beta = "full",
          omega_zeta = "full", delta_zeta = "full",
          kappa_zeta = "full", sigma_zeta = "full",
          lowertri_zeta = "full",
          rho_zeta = "full", SD_zeta = "full",
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

meta_var(...)

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
Parameterization of the contemporaneous (residual innovation) covariance structure. One of \code{"cov"}, \code{"chol"}, \code{"prec"}, \code{"ggm"}, or \code{"cor"}.
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
  \item{rho_zeta}{Only used when \code{contemporaneous = "cor"}. Contemporaneous (innovation) correlation matrix specification (\code{"full"}, \code{"zero"}, or a pattern matrix).}
  \item{SD_zeta}{Only used when \code{contemporaneous = "cor"}. Diagonal matrix of contemporaneous (innovation) standard deviations (\code{"full"}, \code{"diag"}, or a pattern matrix).}
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

\bold{Sample size convention.} When raw \code{data} is supplied, the per-study sample size used in the sampling-error approximation (the 'V' matrices) is the number of rows of the lag-augmented data set, including rows in which the lagged block is missing (e.g., the first measurement of each day or subject). The covariances themselves are computed using pairwise-complete observations. For studies consisting of many short segments (e.g., few beeps per day), this convention slightly overstates the effective number of (current, lagged) pairs for the lagged-covariance elements, and hence slightly understates their sampling variance. When this is a concern, the 'V' matrices can be supplied manually via \code{Vmats}.

\bold{Model size limits.} The number of meta-level outcomes is \code{p(p+1)/2 + p^2}, which grows quickly with the number of variables: the default full random-effects structure has on the order of \code{p^4} parameters, and internal matrices grow with the fourth power of the number of outcomes. Models with more than 11 variables (215 outcomes) exceed internal limits and are refused with an informative error; models close to that limit are very slow and memory-hungry, and require many studies to be estimable (a warning is given when the number of studies is smaller than the number of outcomes). For large numbers of variables, or when the "studies" are in fact subjects of a single intensive longitudinal study, consider the multi-level VAR model \code{\link{ml_var1}} / \code{ml_gvar1} instead.

\bold{Fit indices.} The saturated meta-analytic model (fully free \eqn{\Sigma_0}, \eqn{\Sigma_1} and random effects) has zero degrees of freedom, in which case chi-square based fit indices are not informative and are left blank.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\examples{
\dontrun{
library("dplyr")
set.seed(42)

# Simulate stationary VAR(1) time-series data for several subjects/studies
# directly in base R (no extra packages needed):
nNode  <- 3
nTime  <- 50
nStudy <- 20

# Fixed lag-1 (temporal) coefficient matrix:
beta <- matrix(c( 0.4,  0.1,  0.0,
                  0.0,  0.4,  0.1,
                  0.1,  0.0,  0.4), nNode, nNode, byrow = TRUE)

# Innovation (contemporaneous) covariance:
Sigma <- diag(nNode)

varNames <- paste0("V", 1:nNode)

simVAR1 <- function(id, nTime, beta, Sigma){
  p <- nrow(beta)
  L <- chol(Sigma)
  Y <- matrix(0, nTime, p)
  for (t in 2:nTime){
    Y[t, ] <- beta \%*\% Y[t - 1, ] + as.vector(t(L) \%*\% rnorm(p))
  }
  out <- as.data.frame(Y)
  names(out) <- paste0("V", 1:p)
  out$study <- id
  out
}

Data <- do.call(rbind, lapply(seq_len(nStudy), simVAR1,
                              nTime = nTime, beta = beta, Sigma = Sigma))

# Fit meta-analytic VAR(1):
mod <- meta_var1(Data,
                 vars = varNames,
                 studyvar = "study",
                 estimator = "ML")
mod <- mod \%>\% runmodel

# Inspect results:
print(mod)
getmatrix(mod, "beta")

# Fit meta-analytic GVAR (with GGM contemporaneous):
mod_gvar <- meta_gvar(Data,
                      vars = varNames,
                      studyvar = "study",
                      estimator = "ML")
mod_gvar <- mod_gvar \%>\% runmodel
getmatrix(mod_gvar, "omega_zeta")
}
}
\seealso{
\code{\link{var1}}, \code{\link{meta_varcov}}, \code{\link{meta_lvm}}
}
