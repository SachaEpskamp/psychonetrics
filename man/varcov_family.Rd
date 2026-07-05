\name{varcov}
\alias{varcov}
\alias{cholesky}
\alias{precision}
\alias{prec}
\alias{ggm}
\alias{corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance-covariance family of psychonetrics models
}
\description{
This is the family of models that models only a variance-covariance matrix with mean structure. The \code{type} argument can be used to define what model is used: \code{type = "cov"} (default) models a variance-covariance matrix directly, \code{type = "chol"} (alias: \code{cholesky()}) models a Cholesky decomposition, \code{type = "prec"} (alias: \code{precision()}) models a precision matrix, \code{type = "ggm"} (alias: \code{ggm()}) models a Gaussian graphical model (Epskamp, Rhemtulla and Borsboom, 2017), and \code{type = "cor"} (alias: \code{corr()}) models a correlation matrix.
}
\usage{
varcov(data, type = c("cov", "chol", "prec", "ggm", "cor"),
                   sigma = "full", kappa = "full", omega = "full",
                   lowertri = "full", delta = "diag", rho = "full", SD =
                   "full", mu, tau, vars, ordered = character(0), groups,
                   groupvar, covs, cors, means, nobs, missing = "auto", equal =
                   "none", baseline_saturated = TRUE, estimator =
                   "default", likelihood = c("normal", "wishart"),
                   fixed_x = character(0), sampling_weights,
                   optimizer, storedata = FALSE, WLS.W,
                   sampleStats, meanstructure, corinput, verbose = FALSE,
                   covtype = c("choose", "ML", "UB"), standardize =
                   c("none", "z", "quantile"), fullFIML = FALSE,
                   bootstrap = FALSE, boot_sub, boot_resample,
                   penalty_lambda = NA, penalty_alpha = 1,
                   penalize_matrices)
cholesky(\dots)
precision(\dots)
prec(\dots)
ggm(\dots)
corr(\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{type}{
The type of model used. See description.
}
  \item{sigma}{
Only used when \code{type = "cov"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa}{
Only used when \code{type = "prec"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega}{
Only used when \code{type = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri}{
Only used when \code{type = "chol"}. Either \code{"full"} to estimate every element freely, \code{"diag"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
\item{delta}{
Only used when \code{type = "ggm"}. Either \code{"diag"} or \code{"zero"} (not recommended), or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
\item{rho}{
Only used when \code{type = "cor"}. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
\item{SD}{
Only used when \code{type = "cor"}. Either \code{"diag"} or \code{"zero"} (not recommended), or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{mu}{
Optional vector encoding the mean structure. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free means, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{tau}{
Optional list encoding the thresholds per variable.
}
  \item{vars}{
An optional character vector encoding the variables used in the analysis. Must equal names of the dataset in \code{data}.
}
  \item{groups}{
Deprecated. Use \code{groupvar} instead. An optional string indicating the name of the group variable in \code{data}.
}
  \item{groupvar}{
An optional string indicating the name of the group variable in \code{data}. Replaces the deprecated \code{groups} argument; if both are supplied, \code{groupvar} takes precedence with a warning.
}
  \item{covs}{
A sample variance--covariance matrix, or a list/array of such matrices for multiple groups. Make sure \code{covtype} argument is set correctly to the type of covariances used.
}
  \item{cors}{
A sample correlation matrix, or a list/array of such matrices for multiple groups. When supplied, \code{corinput} defaults to \code{TRUE} and the matrix is used in place of \code{covs}. Requires \code{nobs}.
}
  \item{means}{
A vector of sample means, or a list/matrix containing such vectors for multiple groups.
}
  \item{nobs}{
The number of observations used in \code{covs} and \code{means}, or a vector of such numbers of observations for multiple groups.
}
\item{covtype}{
  If 'covs' is used, this is the type of covariance (maximum likelihood or unbiased) the input covariance matrix represents. Set to \code{"ML"} for maximum likelihood estimates (denominator n) and \code{"UB"} to unbiased estimates (denominator n-1). The default will try to find the type used, by investigating which is most likely to result from integer valued datasets.
}
  \item{missing}{
How should missingness be handled in computing the sample covariances and number of observations when \code{data} is used. Can be \code{"auto"} (default) for automatic detection, \code{"listwise"} for listwise deletion, or \code{"pairwise"} for pairwise deletion. When \code{"auto"}, the function checks for missing data and switches ML to FIML, PML to PFIML, or defaults to listwise for LS estimators.
}
  \item{equal}{
A character vector indicating which matrices should be constrained equal across groups. 
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
  \item{estimator}{
The estimator to be used. Currently implemented are \code{"ML"} for maximum likelihood estimation, \code{"FIML"} for full-information maximum likelihood estimation, \code{"PML"} for penalized maximum likelihood estimation, \code{"PFIML"} for penalized full-information maximum likelihood estimation, \code{"ULS"} for unweighted least squares estimation, \code{"WLS"} for weighted least squares estimation, and \code{"DWLS"} for diagonally weighted least squares estimation. \code{"WLSMV"} is accepted as a synonym for \code{"DWLS"} (both use DWLS estimation with a mean-and-variance adjusted scaled test statistic). When \code{missing = "auto"} (default), \code{"ML"} is automatically switched to \code{"FIML"} and \code{"PML"} to \code{"PFIML"} if missing data is detected. Defaults to \code{"ML"} for continuous data and \code{"DWLS"} when ordinal variables are specified via \code{ordered}.
}
  \item{likelihood}{
The Gaussian likelihood scaling for maximum-likelihood estimation. \code{"normal"} (the default) uses the \eqn{n} (biased) sample-covariance denominator, so the chi-square is \eqn{N \hat{F}} and the parameter covariance is \eqn{\mathrm{Info}^{-1}/N}. \code{"wishart"} uses the \eqn{n-1} (unbiased) sample covariance per group, so the chi-square uses \eqn{(n_g - 1)} multipliers and the standard errors are inflated by \eqn{\sqrt{n_g/(n_g-1)}}, matching \code{lavaan}'s \code{likelihood = "wishart"} (Rosseel, 2012). Only available for complete-data maximum likelihood (\code{estimator = "ML"}); it errors for \code{FIML}, least-squares and ordinal estimators and for raw time-series input. Default \code{"normal"} reproduces the behaviour of previous versions exactly.
}
  \item{fixed_x}{
Character vector of exogenous variable names whose means and mutual (co)variances are fixed to their sample values and excluded from the free-parameter count and the degrees-of-freedom statistic count, matching \code{lavaan}'s \code{sem(..., fixed.x = TRUE)} (Rosseel, 2012). The model is thereby conditioned on these variables, which relaxes the distributional (multivariate-normality) assumption on them; the cross-covariances between the fixed-x and the endogenous variables remain free. Only supported for \code{type = "cov"} and complete-data \code{estimator = "ML"}. The reported log-likelihood is the conditional log-likelihood (of the endogenous variables given x). Default \code{character(0)} (no fixed.x; behaviour unchanged). The baseline and saturated reference models are conditioned on the exogenous block as well, so the incremental fit indices (CFI/TLI/NFI/\ldots) and the absolute fit measures (chi-square, df, RMSEA, AIC, BIC, log-likelihood) and all estimates and standard errors match \code{lavaan}.
}
  \item{sampling_weights}{
Optional single column name in \code{data} giving (survey) SAMPLING weights, enabling pseudo-maximum-likelihood estimation. The moments are computed as weighted means/covariances (weights normalized so that they sum, over all groups, to the total sample size -- \code{lavaan}'s \code{sampling.weights.normalization = "total"}), and the model is fit with robust (\code{MLR}) Huber-White sandwich standard errors and a Yuan-Bentler-Mplus scaled test statistic, exactly matching \code{lavaan}'s \code{sampling.weights=}. Point estimates, standard errors, the scaled chi-square and the robust fit indices agree with \code{lavaan} to numerical precision (single- and multi-group). Requires complete-data continuous input (raw \code{data}); not supported for ordinal data, missing data (FIML), or \code{likelihood = "wishart"}. References: Asparouhov (2005); Rosseel (2012).
}
  \item{optimizer}{
The optimizer to be used. Can be one of \code{"nlminb"} (the default R \code{nlminb} function), \code{"ucminf"} (from the \code{optimr} package), \code{"nloptr_TNEWTON"} (preconditioned truncated Newton via \code{nloptr}), and \code{"LBFGS++"} (pure C++ L-BFGS-B). Defaults to \code{"nlminb"}.
}
  \item{storedata}{
Logical, should the raw data be stored? Needed for bootstrapping (see \code{bootstrap}).
}
  \item{standardize}{
Which standardization method should be used? \code{"none"} (default) for no standardization, \code{"z"} for z-scores, and \code{"quantile"} for a non-parametric transformation to the quantiles of the marginal standard normal distribution.
}
\item{WLS.W}{
Optional WLS weights matrix.
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally. 
}
  \item{verbose}{
Logical, should progress be printed to the console?
}
\item{ordered}{
A vector with strings indicating the variables that are ordered categorical, or set to \code{TRUE} to model all variables as ordered categorical.
}
\item{meanstructure}{
Logical, should the meanstructure be modeled explicitly?
}
\item{corinput}{
Logical, is the input a correlation matrix?
}
  \item{fullFIML}{
Logical, should row-wise FIML be used? Not recommended!
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
\item{penalty_lambda}{Numeric penalty strength for penalized ML estimation (PML/PFIML). Default is \code{NA}, which triggers automatic lambda selection via EBIC grid search when \code{estimator = "PML"} or \code{"PFIML"} (see \code{\link{find_penalized_lambda}}). Set to a specific numeric value (e.g., \code{0.1}) for manual lambda, or \code{0} for no penalty.}
\item{penalty_alpha}{Elastic net mixing parameter: 1 = LASSO (default), 0 = ridge.}
\item{penalize_matrices}{Character vector of matrix names to penalize. If missing, defaults are selected based on the model type.}
\item{\dots}{
Arguments sent to \code{varcov}
}
}
\details{
The model used in this family is:

\eqn{\mathrm{var}(\boldsymbol{y} ) = \boldsymbol{\Sigma}}{var(y) = sigma}

\eqn{\mathcal{E}( \boldsymbol{y} ) = \boldsymbol{\mu}}{E(y) = mu}

in which the covariance matrix can further be modeled in four ways. With \code{type = "chol"} as Cholesky decomposition:

\eqn{\boldsymbol{\Sigma} = \boldsymbol{L}\boldsymbol{L}^{\top}}{sigma = lowertri * lowertri'},

with \code{type = "prec"} as Precision matrix:

\eqn{\boldsymbol{\Sigma} = \boldsymbol{K}^{-1}}{sigma = kappa^(-1)},

with \code{type = "ggm"} as Gaussian graphical model:

\eqn{\boldsymbol{\Sigma} = \boldsymbol{\Delta}(\boldsymbol{I} - \boldsymbol{\Omega})^(-1) \boldsymbol{\Delta}}{sigma = delta * (I - omega)^(-1) * delta},

and finally with \code{type = "cor"} as a correlation matrix scaled by standard deviations:

\eqn{\boldsymbol{\Sigma} = \boldsymbol{D}\boldsymbol{R}\boldsymbol{D}}{sigma = SD * rho * SD},

in which \eqn{\boldsymbol{D}} is a diagonal matrix of standard deviations and \eqn{\boldsymbol{R}} a correlation matrix.
}
\value{
An object of the class psychonetrics 
}
\references{
Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017). Generalized network psychometrics: Combining network and latent variable models. Psychometrika, 82(4), 904-927.

Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. \emph{Journal of Statistical Software, 48}(2), 1--36. \doi{10.18637/jss.v048.i02}
}
\author{
Sacha Epskamp
}


\seealso{
\code{\link{lvm}}, \code{\link{var1}}, \code{\link{dlvm1}},
\code{\link{penalize}} for manual penalty control,
\code{\link{find_penalized_lambda}} for automatic lambda selection,
\code{\link{refit}} for post-selection inference after penalized estimation.
}
\examples{
# Load bfi data from psych package:
library("psychTools")
data(bfi)

# Also load dplyr for the pipe operator:
library("dplyr")

# Let's take the agreeableness items, and gender:
ConsData <- bfi \%>\% 
  select(A1:A5, gender) \%>\% 
  na.omit # Let's remove missingness (or use missing = "auto" default which auto-selects FIML)

# Define variables:
vars <- names(ConsData)[1:5]

# Saturated estimation:
mod_saturated <- ggm(ConsData, vars = vars)

# Run the model:
mod_saturated <- mod_saturated \%>\% runmodel

# We can look at the parameters:
mod_saturated \%>\% parameters

# Labels:
labels <- c(
  "indifferent to the feelings of others",
  "inquire about others' well-being",
  "comfort others",
  "love children",
  "make people feel at ease")
  
# Plot CIs:
CIplot(mod_saturated,  "omega", labels = labels, labelstart = 0.2)


\donttest{
# We can also fit an empty network:
mod0 <- ggm(ConsData, vars = vars, omega = "zero")

# Run the model:
mod0 <- mod0 \%>\% runmodel

# We can look at the modification indices:
mod0 \%>\% MIs

# To automatically add along modification indices, we can use stepup:
mod1 <- mod0 \%>\% stepup

# Let's also prune all non-significant edges to finish:
mod1 <- mod1 \%>\% prune

# Look at the fit:
mod1 \%>\% fit

# Compare to original (baseline) model:
compare(baseline = mod0, adjusted = mod1)

# We can also look at the parameters:
mod1 \%>\% parameters

# Or obtain the network as follows:
getmatrix(mod1, "omega")

# Penalized GGM estimation with automatic lambda selection:
mod_pml <- ggm(ConsData, vars = vars, estimator = "PML")
mod_pml <- mod_pml \%>\% runmodel

# Check selected lambda:
mod_pml@optim$lambda_search

# Obtain the sparse network:
getmatrix(mod_pml, "omega")

# Wishart likelihood (n-1 denominator), matching lavaan likelihood = "wishart":
mod_wis <- varcov(ConsData, vars = vars, type = "cov", likelihood = "wishart")
mod_wis <- mod_wis \%>\% runmodel

# fixed.x: condition on A1 and A2 as exogenous covariates (their means and
# mutual covariances are fixed to the sample values, excluded from npar/df):
mod_fx <- varcov(ConsData, vars = vars, type = "cov",
                 fixed_x = c("A1", "A2"))
mod_fx <- mod_fx \%>\% runmodel
}
}
