\name{panelvar}
\alias{panelvar}
\alias{panelgvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lag-1 panel vector autoregression family of psychonetrics models
}
\description{
This is the family of models that models a lag-1 vector autoregression (VAR) process with random intercepts on panel data of observed variables: a within-person stationary VAR(1) process (temporal matrix \code{beta} and contemporaneous structure \code{sigma_zeta_within}) combined with a between-person (random intercept) structure (\code{sigma_zeta_between}). The \code{panelvar} function models the contemporaneous and between-person structures as variance--covariance matrices by default, and the \code{panelgvar} wrapper models them as Gaussian graphical models (GGMs), yielding the panel graphical VAR model.

This model family is the observed-variable special case of \code{\link{dlvm1}} (a factor loading matrix fixed to identity, zero residual variances and observed means). Since psychonetrics 0.16.2 it is implemented as its own model framework: results are identical to routing through \code{dlvm1} with dummy matrices, but estimation skips all factor-loading algebra and is faster. The mean structure is a single \code{mu} matrix of observed stationary means (in older versions these appeared as \code{nu}). Latent panel models (\code{lambda} not identity) are still handled by \code{\link{dlvm1}} and its wrapper \code{panellvgvar}.

All functions accept both wide-format data (with a design matrix for \code{vars}) and long-format data (with a character vector for \code{vars} plus \code{idvar} and optionally \code{beepvar}). The format is auto-detected from the type of \code{vars}.
}
\usage{
panelvar(data, vars, within_latent = c("cov", "chol", "prec",
                   "ggm", "cor"), between_latent = c("cov", "chol",
                   "prec", "ggm", "cor"), datatype = c("auto", "wide",
                   "long"), idvar, beepvar, standardize = c("none",
                   "z", "quantile", "z_per_wave"), beta = "full",
                   omega_zeta_within = "full", delta_zeta_within =
                   "diag", kappa_zeta_within = "full",
                   sigma_zeta_within = "full", lowertri_zeta_within =
                   "full", omega_zeta_between = "full",
                   delta_zeta_between = "diag", kappa_zeta_between =
                   "full", sigma_zeta_between = "full",
                   lowertri_zeta_between = "full", mu, groups,
                   groupvar, covs, cors, means, nobs, corinput, start
                   = "version2", covtype = c("choose", "ML", "UB"),
                   missing = "auto", equal = "none",
                   baseline_saturated = TRUE, estimator = "ML",
                   optimizer, storedata = FALSE, verbose = FALSE,
                   sampleStats, baseline =
                   c("stationary_random_intercept", "stationary",
                   "independence", "none"), bootstrap = FALSE, boot_sub,
                   boot_resample, within, between,
                   penalty_lambda = NA, penalty_alpha = 1,
                   penalize_matrices,
                   rho_zeta_within = "full", SD_zeta_within = "diag",
                   rho_zeta_between = "full", SD_zeta_between = "diag",
                   temporal = c("raw","PDC"), PDC = "full")

panelgvar(data, vars, within_latent = c("ggm","chol","cov","prec","cor"),
          between_latent = c("ggm","chol","cov","prec","cor"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{vars}{ Required argument.
For wide-format data, this must be a *matrix* with each row indicating a variable and each column indicating a measurement. The matrix must be filled with names of the variables in the dataset corresponding to variable i at wave j. NAs can be used to indicate missing waves. The rownames of this matrix will be used as variable names. For long-format data, this should be a character vector of variable names in the dataset.
}
  \item{within_latent}{
The type of within-person contemporaneous model to be used. Since there are no latent variables in panelvar, the shorter alias \code{within} is preferred.
}
  \item{between_latent}{
The type of between-person model to be used. Since there are no latent variables in panelvar, the shorter alias \code{between} is preferred.
}
  \item{datatype}{
Data format: \code{"auto"} (default) auto-detects based on whether \code{vars} is a matrix (wide) or character vector (long). \code{"wide"} forces wide-format interpretation. \code{"long"} forces long-format interpretation with automatic reshape to wide.
}
  \item{idvar}{
Optional string indicating the subject/cluster ID variable in \code{data}. Required when \code{datatype = "long"} or when \code{vars} is a character vector.
}
  \item{beepvar}{
Optional string indicating the time point / measurement occasion variable. If missing when \code{datatype = "long"}, a sequential measurement variable is created per subject.
}
  \item{standardize}{
Standardization method: \code{"none"} (default), \code{"z"} for global z-scores per variable across all waves (uses the overall mean and standard deviation for each variable computed across all waves of data), \code{"quantile"} for quantile normalization across all waves, \code{"z_per_wave"} for independent z-scores per wave column (note: this imposes stationarity as means become 0 and variances become 1 at each wave).
}
  \item{beta}{
A model matrix encoding the temporal relationships (transpose of temporal network). A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix. Can also be \code{"full"} for a full temporal network or \code{"zero"} for an empty temporal network.
}
  \item{omega_zeta_within}{
Only used when \code{within_latent = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta_within}{
Only used when \code{within_latent = "ggm"}. Can be \code{"diag"}, \code{"zero"} (not recommended), or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta_within}{
Only used when \code{within_latent = "prec"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_zeta_within}{
Only used when \code{within_latent = "cov"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta_within}{
Only used when \code{within_latent = "chol"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_zeta_between}{
Only used when \code{between_latent = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta_between}{
Only used when \code{between_latent = "ggm"}. Can be \code{"diag"}, \code{"zero"} (not recommended), or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta_between}{
Only used when \code{between_latent = "prec"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_zeta_between}{
Only used when \code{between_latent = "cov"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta_between}{
Only used when \code{between_latent = "chol"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{mu}{
Optional vector encoding the observed stationary means. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free means, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{groups}{
Deprecated. Use \code{groupvar} instead. An optional string indicating the name of the group variable in \code{data}.
}
  \item{groupvar}{
An optional string indicating the name of the group variable in \code{data}. Replaces the deprecated \code{groups} argument; if both are supplied, \code{groupvar} takes precedence with a warning.
}
  \item{covs}{
A sample variance--covariance matrix, or a list/array of such matrices for multiple groups. Make sure the \code{covtype} argument is set correctly to the type of covariances used.
}
  \item{cors}{
Not supported in \code{panelvar()}; supplying correlation input produces an error.
}
  \item{means}{
A vector of sample means, or a list/matrix containing such vectors for multiple groups.
}
  \item{nobs}{
The number of observations used in \code{covs} and \code{means}, or a vector of such numbers of observations for multiple groups.
}
  \item{corinput}{
Logical. Not supported in \code{panelvar()} and will produce an error if set to \code{TRUE}.
}
\item{start}{Start value specification. Can be either a string or a psychonetrics model. If it is a string, \code{"version2"} indicates the latest version of start value computation, \code{"version1"} indicates start values as they were computed up to version 0.11, and \code{"simple"} indicate simple starting values. If this is a psychonetrics model the starting values will be based on the output. This can be useful, for example, if you first estimate a model with matrices set to a Cholesky decomposition, then use those values as start values for estimating Gaussian graphical models.}
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
The estimator to be used. Currently implemented are \code{"ML"} for maximum likelihood estimation, \code{"FIML"} for full-information maximum likelihood estimation, \code{"PML"} for penalized maximum likelihood estimation, \code{"PFIML"} for penalized full-information maximum likelihood estimation, \code{"ULS"} for unweighted least squares estimation, \code{"WLS"} for weighted least squares estimation, and \code{"DWLS"} for diagonally weighted least squares estimation. When \code{missing = "auto"} (default), \code{"ML"} is automatically switched to \code{"FIML"} and \code{"PML"} to \code{"PFIML"} if missing data is detected.
}
  \item{optimizer}{
The optimizer to be used. Can be one of \code{"nlminb"} (the default R \code{nlminb} function), \code{"ucminf"} (from the \code{optimr} package), \code{"nloptr_TNEWTON"} (preconditioned truncated Newton via \code{nloptr}), and \code{"LBFGS++"} (pure C++ L-BFGS-B). Defaults to \code{"nlminb"}.
}
  \item{storedata}{
Logical, should the raw data be stored? Needed for bootstrapping (see \code{\link{bootstrap}}).
}
  \item{verbose}{
Logical, should progress be printed to the console?
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally.
}
\item{baseline}{What baseline model should be used? \code{"stationary_random_intercept"} includes both within- and between person variances constrained equal across time (default), \code{"stationary"} only includes within-person variances constrained equal across time, \code{"independence"} (default up to version 0.11) includes a variance for every variable at every time point (not constrained equal across time), and \code{"none"} includes no baseline model.}
\item{bootstrap}{
Should the data be bootstrapped? If \code{TRUE} the data are resampled and a bootstrap sample is created. These must be aggregated using \code{\link{aggregate_bootstraps}}! Can be \code{TRUE} or \code{FALSE}. Can also be \code{"nonparametric"} (which sets \code{boot_sub = 1} and \code{boot_resample = TRUE}) or \code{"case"} (which sets \code{boot_sub = 0.75} and \code{boot_resample = FALSE}).
}
\item{boot_sub}{
Proportion of cases to be subsampled (\code{round(boot_sub * N)}).
}
\item{boot_resample}{
Logical, should the bootstrap be with replacement (\code{TRUE}) or without replacement (\code{FALSE})
}
\item{within}{
Preferred name for \code{within_latent} (there are no latent variables in panelvar). If both \code{within} and \code{within_latent} are explicitly specified, \code{within_latent} takes precedence and a warning is issued.
}
\item{between}{
Preferred name for \code{between_latent} (there are no latent variables in panelvar). If both \code{between} and \code{between_latent} are explicitly specified, \code{between_latent} takes precedence and a warning is issued.
}
\item{penalty_lambda}{Numeric penalty strength for penalized ML estimation (PML/PFIML). \code{NA} (default) triggers automatic selection via EBIC-based grid search when a penalized estimator is used; set to a specific numeric value to use a fixed penalty strength (0 = no penalty). See \code{\link{find_penalized_lambda}} and \code{\link{penalize}}.}
\item{penalty_alpha}{Elastic net mixing parameter: 1 = LASSO (default), 0 = ridge.}
\item{penalize_matrices}{Character vector of matrix names to penalize. If missing, defaults are selected based on the model type.}
\item{rho_zeta_within}{Only used when \code{within_latent = "cor"}. Correlation matrix specification (\code{"full"}, \code{"zero"}, or a pattern matrix; see \code{omega_zeta_within}).}
\item{SD_zeta_within}{Only used when \code{within_latent = "cor"}. Diagonal standard-deviation matrix specification (\code{"diag"}, \code{"full"}, or a pattern matrix; see \code{delta_zeta_within}).}
\item{rho_zeta_between}{Only used when \code{between_latent = "cor"}. Correlation matrix specification (\code{"full"}, \code{"zero"}, or a pattern matrix).}
\item{SD_zeta_between}{Only used when \code{between_latent = "cor"}. Diagonal standard-deviation matrix specification (\code{"diag"}, \code{"full"}, or a pattern matrix).}
\item{temporal}{Parameterization of the temporal effects: \code{"raw"} (default) models \code{beta} directly; \code{"PDC"} models the partial directed correlations directly, with the \code{PDC} matrix encoding from = row and to = column (the transposed orientation relative to \code{beta}, matching \code{getmatrix(..., "PDC")}). A (pattern-constrained) PDC model is fit-equivalent to the corresponding raw model, but standard errors, confidence intervals, equality constraints and modification indices then apply directly on the PDC scale.}
\item{PDC}{Only used when \code{temporal = "PDC"}. Model matrix for the partial directed correlations (from = row, to = column): \code{"full"}, \code{"zero"}, \code{"diag"}, or a pattern matrix with 0s indicating elements fixed to zero, 1s indicating free elements, and higher integers indicating equality constraints.}
\item{\dots}{
Arguments sent to \code{panelvar}.
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}

\author{
Sacha Epskamp
}

\seealso{
\code{\link{dlvm1}} for the latent variable generalization of this model family.
}

\examples{
library("dplyr")

# Smoke data cov matrix, based on LISS data panel https://www.dataarchive.lissdata.nl
smoke <- structure(c(47.2361758611759, 43.5366809116809, 41.0057465682466,
                     43.5366809116809, 57.9789886039886, 47.6992521367521,
                     41.0057465682466,
                     47.6992521367521, 53.0669434731935), .Dim = c(3L, 3L),
                   .Dimnames = list(
                       c("smoke2008", "smoke2009", "smoke2010"), c("smoke2008",
                   "smoke2009", "smoke2010")))

# Design matrix:
design <- matrix(rownames(smoke),1,3)

# Form model:
mod <- panelvar(vars = design,
                covs = smoke, nobs = 352
)

\donttest{
# Run model:
mod <- mod \%>\% runmodel

# Evaluate fit:
mod \%>\% fit
}

# panelgvar and panelvar also accept long-format data.
# For long-format data, pass vars as a character vector
# and provide idvar (and optionally beepvar):
\donttest{
# Simulate some long-format panel data:
longdata <- data.frame(
  id = rep(1:50, each = 3),
  time = rep(1:3, 50),
  x = rnorm(150),
  y = rnorm(150)
)

# Use panelgvar with long-format data:
mod_long <- panelgvar(data = longdata,
                      vars = c("x", "y"),
                      idvar = "id",
                      beepvar = "time")
}
}
