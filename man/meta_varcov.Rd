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
meta_varcov(cors, nobs, data, covs, studyvar, groups, groupvar,
                   corinput, Vmats, Vmethod = c("individual", "pooled",
                   "OS_individual", "OS_pooled",
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
                   c("FIML", "ML"), sampleStats, verbose = FALSE,
                   bootstrap = FALSE, boot_sub, boot_resample)
  
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
  \item{data}{
Optional data frame. When supplied together with \code{studyvar}, correlation matrices and sample sizes are computed internally per study. Cannot be used together with \code{cors}, \code{covs}, or \code{nobs}.
}
  \item{covs}{
A list of covariance matrices. Alternative to \code{cors}; when supplied, \code{corinput} defaults to \code{FALSE}.
}
  \item{studyvar}{
A string indicating the column name in \code{data} that identifies the study. Required when \code{data} is supplied.
}
  \item{groups}{
Deprecated. Use \code{groupvar} instead. Multi-group support is not yet included for meta-analytic models.
}
  \item{groupvar}{
Not yet supported for meta-analytic models. Supplying this argument will produce an error.
}
  \item{corinput}{
Logical. Defaults to \code{TRUE} when \code{cors} is used, \code{FALSE} when \code{covs} is used. Controls whether the input is treated as correlation matrices.
}
  \item{Vmats}{
Optional list with 'V' matrices (sampling error variance approximations). May not be combined with an explicitly supplied \code{Vmethod} (which would otherwise silently be ignored).
}
  \item{Vmethod}{
Which method should be used to approximate the sampling error variance? One of \code{"individual"} (the default; the sampling-error covariance matrix is approximated separately for each study), \code{"pooled"} (a single pooled approximation is used for all studies), \code{"OS_individual"} or \code{"OS_pooled"} (the asymptotic sampling covariance matrix of the sample correlations following Olkin and Siotani, 1976, evaluated per study or at the sample-size-weighted pooled correlation matrix; for covariance input the corresponding Wishart expression is used), \code{"metaSEM_individual"}, or \code{"metaSEM_weighted"}. The two \code{metaSEM_*} options compute the sampling-error covariances via \code{metaSEM::asyCov} and therefore require the Suggests package \code{metaSEM} to be installed (per study or using sample-size weighting, respectively). See Details for the difference between the default methods and the \code{OS_*} methods.
}
\item{Vestimation}{
How should the sampling error estimates be evaluated? Either \code{"averaged"} (the default; a single sampling-error covariance matrix is averaged over the studies) or \code{"per_study"} (a separate sampling-error covariance matrix is retained for each study).
}
  \item{type}{
What to model? Currently only \code{"cor"} and \code{"ggm"} are supported.
}
  \item{sigma_y}{
Reserved covariance-matrix specification. Not used for the currently supported values of \code{type} (\code{"cor"} and \code{"ggm"}); retained for internal use and forward compatibility.
}
  \item{kappa_y}{
Reserved precision-matrix specification. Not used for the currently supported values of \code{type} (\code{"cor"} and \code{"ggm"}); retained for internal use and forward compatibility.
}
  \item{omega_y}{
Only used when \code{type = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_y}{
Reserved Cholesky-factor specification. Not used for the currently supported values of \code{type} (\code{"cor"} and \code{"ggm"}); retained for internal use and forward compatibility.
}
  \item{delta_y}{
Only used when \code{type = "ggm"}. Either \code{"diag"} or \code{"zero"} (not recommended), or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{rho_y}{
Only used when \code{type = "cor"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{SD_y}{
Only used when \code{type = "cor"}. Either \code{"diag"} or \code{"zero"}, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{randomEffects}{
What to model for the random effects? 
}
  \item{sigma_randomEffects}{
Only used when \code{type = "cov"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_randomEffects}{
Only used when \code{randomEffects = "prec"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_randomEffects}{
Only used when \code{randomEffects = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_randomEffects}{
Only used when \code{randomEffects = "chol"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_randomEffects}{
Only used when \code{randomEffects = "ggm"}. Either \code{"diag"} or \code{"zero"}, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{rho_randomEffects}{
Only used when \code{randomEffects = "cor"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{SD_randomEffects}{
Only used when \code{randomEffects = "cor"}. Either \code{"diag"} or \code{"zero"}, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{vars}{
Variables to be included.
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
  \item{optimizer}{
The optimizer to be used. Can be one of \code{"nlminb"} (the default R \code{nlminb} function), \code{"ucminf"} (from the \code{optimr} package), \code{"nloptr_TNEWTON"} (preconditioned truncated Newton via \code{nloptr}), and \code{"LBFGS++"} (pure C++ L-BFGS-B). Defaults to \code{"nlminb"}.
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
\item{bootstrap}{
Should the data be bootstrapped? If \code{TRUE} the data are resampled and a bootstrap sample is created. These must be aggregated using \code{\link{aggregate_bootstraps}}! Can be \code{TRUE} or \code{FALSE}. Can also be \code{"nonparametric"} (which sets \code{boot_sub = 1} and \code{boot_resample = TRUE}) or \code{"case"} (which sets \code{boot_sub = 0.75} and \code{boot_resample = FALSE}).
}
\item{boot_sub}{
Proportion of cases to be subsampled (\code{round(boot_sub * N)}).
}
\item{boot_resample}{
Logical, should the bootstrap be with replacement (\code{TRUE}) or without replacement (\code{FALSE})
}
\item{\dots}{
Arguments sent to \code{meta_varcov}
}
}

\details{
The default sampling-error approximations \code{Vmethod = "individual"} and \code{"pooled"} invert the Fisher information of a saturated correlation model (with unit diagonal constraint). This is the asymptotic covariance matrix of the maximum-likelihood estimates of a constrained correlation model, which is smaller than the sampling covariance of the raw sample correlations that enter the meta-level data (for a single correlation: (1 - rho^2)^2 / ((1 + rho^2) n) instead of (1 - rho^2)^2 / n). The alternative options \code{Vmethod = "OS_individual"} and \code{"OS_pooled"} instead use the asymptotic sampling covariance matrix of the sample correlations themselves (Olkin and Siotani, 1976), which matches \code{metaSEM::asyCov} and may be preferred when the absolute size of the random-effects (heterogeneity) variances is of interest: with the default methods part of the sampling variability can be absorbed into the random-effects variances, inflating them somewhat. The pooled correlation estimates themselves are typically nearly identical between these choices. The default is kept at \code{"individual"} for backward compatibility with published results.

Note that the saturated meta-analytic model (a fully free pooled correlation structure with fully free random effects) has zero degrees of freedom, in which case chi-square based fit indices (e.g., RMSEA, CFI) are not informative and are left blank.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\references{
Jak, S., and Cheung, M. W. L. (2019). Meta-analytic structural equation modeling with moderating effects on SEM parameters. Psychological methods.

Olkin, I., and Siotani, M. (1976). Asymptotic distribution of functions of a correlation matrix. In S. Ikeda (Ed.), Essays in probability and statistics (pp. 235-251). Shinko Tsusho.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\examples{
\donttest{
# Simulate a list of correlation matrices for several studies:
set.seed(123)
nStudy <- 8
nNode  <- 4
trueRho <- diag(nNode)
trueRho[lower.tri(trueRho)] <- trueRho[upper.tri(trueRho)] <- 0.3
nobs <- sample(100:300, nStudy)
cors <- lapply(seq_len(nStudy), function(i){
  X <- MASS::mvrnorm(nobs[i], rep(0, nNode), trueRho)
  R <- cor(X)
  colnames(R) <- rownames(R) <- paste0("V", 1:nNode)
  R
})

# Meta-analytic GGM using the list of correlation matrices and sample sizes:
mod <- meta_varcov(cors = cors, nobs = nobs, type = "ggm")
mod <- runmodel(mod)

# Pooled network:
getmatrix(mod, "omega_y")
}
}
