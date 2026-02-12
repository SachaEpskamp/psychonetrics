\name{meta_lvm}
\alias{meta_lvm}
\title{
Meta-analytic latent variable model
}
\description{
Single-stage meta-analytic latent variable model (LVM) for correlation matrices. Combines a latent variable model (CFA/SEM) for the mean structure with a random-effects model for between-study heterogeneity. Input is always correlation matrices. Based on the MAGNA framework (Epskamp et al., 2024).
}
\usage{
meta_lvm(cors, nobs, Vmats, Vmethod = c("individual", "pooled",
                   "metaSEM_individual", "metaSEM_weighted"), Vestimation
                   = c("averaged", "per_study"),
                   lambda, beta = "zero",
                   latent = c("cov", "chol", "prec", "ggm"),
                   sigma_zeta = "full", kappa_zeta = "full",
                   omega_zeta = "full", lowertri_zeta = "full",
                   delta_zeta = "full",
                   residual = c("cov", "chol", "prec", "ggm"),
                   sigma_epsilon = "diag", kappa_epsilon = "diag",
                   omega_epsilon = "zero", lowertri_epsilon = "diag",
                   delta_epsilon = "diag",
                   identification = "variance",
                   randomEffects = c("chol", "cov", "prec", "ggm", "cor"),
                   sigma_randomEffects = "full",
                   kappa_randomEffects = "full",
                   omega_randomEffects = "full",
                   lowertri_randomEffects = "full",
                   delta_randomEffects = "full",
                   rho_randomEffects = "full",
                   SD_randomEffects = "full",
                   vars, latents,
                   baseline_saturated = TRUE, optimizer,
                   estimator = c("FIML", "ML"),
                   sampleStats, verbose = FALSE,
                   bootstrap = FALSE, boot_sub, boot_resample)
}
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
Which method should be used to approximate the sampling error variance?
}
\item{Vestimation}{
How should the sampling error estimates be evaluated?
}
  \item{lambda}{
A matrix encoding the factor loading structure. Can be a pattern matrix with \code{TRUE}/\code{FALSE} or 0/1 entries, or a matrix with starting values and higher integer codes for equality constraints.
}
  \item{beta}{
Structural (regression) matrix among latent variables. Defaults to \code{"zero"}.
}
  \item{latent}{
Parameterization of the latent covariance structure. One of \code{"cov"}, \code{"chol"}, \code{"prec"}, or \code{"ggm"}.
}
  \item{sigma_zeta}{
Latent covariance matrix specification (used when \code{latent = "cov"}).
}
  \item{kappa_zeta}{
Latent precision matrix specification (used when \code{latent = "prec"}).
}
  \item{omega_zeta}{
Latent partial correlation matrix specification (used when \code{latent = "ggm"}).
}
  \item{lowertri_zeta}{
Latent Cholesky factor specification (used when \code{latent = "chol"}).
}
  \item{delta_zeta}{
Latent scaling matrix specification (used when \code{latent = "ggm"}).
}
  \item{residual}{
Parameterization of the residual covariance structure. One of \code{"cov"}, \code{"chol"}, \code{"prec"}, or \code{"ggm"}.
}
  \item{sigma_epsilon}{
Residual covariance matrix specification (used when \code{residual = "cov"}). Defaults to \code{"diag"}. Diagonal elements are always fixed (derived from the correlation constraint).
}
  \item{kappa_epsilon}{
Residual precision matrix specification (used when \code{residual = "prec"}).
}
  \item{omega_epsilon}{
Residual partial correlation matrix specification (used when \code{residual = "ggm"}).
}
  \item{lowertri_epsilon}{
Residual Cholesky factor specification (used when \code{residual = "chol"}).
}
  \item{delta_epsilon}{
Residual scaling matrix specification (used when \code{residual = "ggm"}).
}
  \item{identification}{
How to identify the model. Only \code{"variance"} is supported for correlation input.
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
  \item{vars}{
Character vector of observed variable names. If missing, names are taken from the correlation matrices.
}
  \item{latents}{
Character vector of latent variable names.
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
}
\details{
This function implements a single-stage meta-analytic latent variable model. The model specifies that the mean of the vectorized correlation matrices follows an LVM-implied correlation structure:

\deqn{\mu = \textrm{vechs}(\Lambda (I-B)^{-1} \Sigma_\zeta (I-B)^{-\top} \Lambda' + \Sigma_\varepsilon)}

where diagonal elements of \eqn{\Sigma_\varepsilon} are derived from the constraint that \eqn{\mu} represents correlations (diagonal of implied covariance = 1).

Between-study heterogeneity is modeled as:

\deqn{\Sigma = \Sigma^{(\textrm{ran})} + V}

where \eqn{V} is the known sampling error covariance and \eqn{\Sigma^{(\textrm{ran})}} captures true between-study variation.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\references{
Epskamp, S. et al. (2024). Meta-Analytic Gaussian Network Aggregation. \emph{Psychometrika}.

Jak, S., and Cheung, M. W. L. (2019). Meta-analytic structural equation modeling with moderating effects on SEM parameters. \emph{Psychological Methods}.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\examples{
\dontrun{
# Generate simulated data
set.seed(42)
nStudies <- 10
nvar <- 6
nlat <- 2

# True factor loading matrix
lambda_true <- matrix(0, nvar, nlat)
lambda_true[1:3, 1] <- c(0.7, 0.8, 0.6)
lambda_true[4:6, 2] <- c(0.7, 0.8, 0.6)

# Generate correlation matrices per study
cors <- list()
nobs <- sample(100:400, nStudies)
for (i in 1:nStudies) {
  sigma <- lambda_true \%*\% t(lambda_true) + diag(1 - rowSums(lambda_true^2))
  dat <- MASS::mvrnorm(nobs[i], rep(0, nvar), sigma)
  cors[[i]] <- cor(dat)
  colnames(cors[[i]]) <- rownames(cors[[i]]) <- paste0("V", 1:nvar)
}

# Fit meta-analytic CFA
mod <- meta_lvm(cors = cors, nobs = nobs, lambda = lambda_true != 0)
mod <- mod \%>\% runmodel

# Inspect results
print(mod)
parameters(mod)
getmatrix(mod, "lambda")
}
}
\seealso{
\code{\link{meta_varcov}}, \code{\link{lvm}}
}
