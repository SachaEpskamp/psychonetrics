\name{dlvm1}
\alias{dlvm1}
\alias{dlvm}
\alias{panel_lvgvar}
\alias{panellvgvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lag-1 dynamic latent variable model family of psychonetrics models for panel data
}
\description{
This is the family of models that models a dynamic factor model on panel data. All functions accept both wide-format data (with a design matrix for \code{vars}) and long-format data (with a character vector for \code{vars} plus \code{idvar} and optionally \code{beepvar}). The format is auto-detected from the type of \code{vars}. There are four covariance structures that can be modeled in different ways: \code{within_latent}, \code{between_latent} for the within-person and between-person latent (contemporaneous) models respectively, and \code{within_residual}, \code{between_residual} for the within-person and between-person residual models respectively. The \code{panellvgvar} wrapper automatically models all latent networks as GGMs. The old name \code{panel_lvgvar} is deprecated in favor of \code{panellvgvar}. \code{dlvm()} is the primary short name of \code{dlvm1()} (identical function; the historical "1" referred to the lag order). The observed-variable special case of this family (a factor loading matrix fixed to identity and zero residual variances) is the panel VAR model, implemented as its own faster framework in \code{\link{panelvar}} and \code{\link{panelgvar}}. The mean structure is stationary by default (item intercepts \code{nu} invariant across waves, latent means fixed to zero); \code{stationary_means = FALSE} makes the intercepts wave-specific, giving the trend-stationary case (see the details below).
}
\usage{
dlvm1(data, vars, datatype = c("auto", "wide", "long"),
                   idvar, beepvar,
                   standardize = c("none", "z", "quantile", "z_per_wave"),
                   lambda, within_latent = c("cov", "chol",
                   "prec", "ggm", "cor"), within_residual = c("cov", "chol",
                   "prec", "ggm", "cor"), between_latent = c("cov", "chol",
                   "prec", "ggm", "cor"), between_residual = c("cov", "chol",
                   "prec", "ggm", "cor"), beta = "full", omega_zeta_within =
                   "full", delta_zeta_within = "diag", kappa_zeta_within
                   = "full", sigma_zeta_within = "full",
                   lowertri_zeta_within = "full", omega_epsilon_within =
                   "zero", delta_epsilon_within = "diag",
                   kappa_epsilon_within = "diag", sigma_epsilon_within =
                   "diag", lowertri_epsilon_within = "diag",
                   omega_zeta_between = "full", delta_zeta_between =
                   "diag", kappa_zeta_between = "full",
                   sigma_zeta_between = "full", lowertri_zeta_between =
                   "full", omega_epsilon_between = "zero",
                   delta_epsilon_between = "diag", kappa_epsilon_between
                   = "diag", sigma_epsilon_between = "diag",
                   lowertri_epsilon_between = "diag", nu, mu_eta,
                   identify = TRUE, identification = c("loadings",
                   "variance"), latents, groups, groupvar, covs,
                   cors, means, nobs, corinput, start
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
                   rho_epsilon_within = "zero", SD_epsilon_within = "diag",
                   rho_zeta_between = "full", SD_zeta_between = "diag",
                   rho_epsilon_between = "zero", SD_epsilon_between = "diag",
                   beta_epsilon = "zero",
                   temporal_latent = c("raw","PDC"),
                   temporal_residual = c("raw","PDC"),
                   PDC = "full", PDC_epsilon = "zero",
                   stationary_means = TRUE)

dlvm(...)

panellvgvar(...)

panel_lvgvar(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{vars}{ Required argument.
For wide-format data, this must be a *matrix* with each row indicating a variable and each column indicating a measurement. The matrix must be filled with names of the variables in the dataset corresponding to variable i at wave j. NAs can be used to indicate missing waves. The rownames of this matrix will be used as variable names. For long-format data, this should be a character vector of variable names in the dataset.
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
  \item{lambda}{ Required argument.
A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{within_latent}{
The type of within-person latent contemporaneous model to be used.
}
  \item{within_residual}{
The type of within-person residual model to be used.
}
  \item{between_latent}{
The type of between-person latent model to be used.
}
  \item{between_residual}{
The type of between-person residual model to be used.
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
  \item{omega_epsilon_within}{
Only used when \code{within_residual = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon_within}{
Only used when \code{within_residual = "ggm"}. Can be \code{"diag"}, \code{"zero"} (not recommended), or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon_within}{
Only used when \code{within_residual = "prec"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon_within}{
Only used when \code{within_residual = "cov"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon_within}{
Only used when \code{within_residual = "chol"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
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
  \item{omega_epsilon_between}{
Only used when \code{between_residual = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon_between}{
Only used when \code{between_residual = "ggm"}. Can be \code{"diag"}, \code{"zero"} (not recommended), or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon_between}{
Only used when \code{between_residual = "prec"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon_between}{
Only used when \code{between_residual = "cov"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon_between}{
Only used when \code{between_residual = "chol"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{nu}{
Optional vector encoding the intercepts of the observed variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector. When \code{stationary_means = FALSE} this is instead a \code{nVar} by \code{nTime} matrix (variables in the rows, waves in the columns), or an array/list of such matrices for multiple groups; elements at variable--wave combinations that are not in the \code{vars} design matrix are always fixed to zero.
}
  \item{mu_eta}{
Optional vector encoding the means of the latent variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector. The latent means are fixed to zero by \code{identify()} for identification, with a stationary as well as with a wave-specific mean structure.
}
  \item{identify}{
Logical, should the model be automatically identified?
}
  \item{identification}{
Type of identification used. \code{"loadings"} to fix the first factor loadings to 1, and \code{"variance"} to fix the diagonal of the latent variable model matrix (sigma_zeta, lowertri_zeta, delta_zeta or kappa_zeta) to 1.
}
  \item{latents}{
An optional character vector with names of the latent variables.
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
A sample correlation matrix, or a list/array of such matrices for multiple groups. When supplied, the matrix is treated as a covariance matrix with a warning (appropriate for standardized data). Requires \code{nobs}. Note: \code{corinput = TRUE} is not supported in \code{dlvm1}.
}
  \item{means}{
A vector of sample means, or a list/matrix containing such vectors for multiple groups. This argument always refers to the sample means; the mean \emph{structure} of the model is set with \code{stationary_means}.
}
  \item{nobs}{
The number of observations used in \code{covs} and \code{means}, or a vector of such numbers of observations for multiple groups.
}
  \item{corinput}{
Logical. Not supported in \code{dlvm1()} and will produce an error if set to \code{TRUE}.
}
\item{start}{Start value specification. Can be either a string or a psychonetrics model. If it is a string, \code{"version2"} indicates the latest version of start value computation, \code{"version1"} indicates start values as they were computed up to version 0.11, and \code{"simple"} indicate simple starting values. If this is a psychonetrics model the starting values will be based on the output. This can be useful, for example, if you first estimate a model with matrices set to a Cholesky decomposition, then use those values as start values for estimating Gaussian graphical models.}
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
\item{covtype}{
  If 'covs' is used, this is the type of covariance (maximum likelihood or unbiased) the input covariance matrix represents. Set to \code{"ML"} for maximum likelihood estimates (denominator n) and \code{"UB"} to unbiased estimates (denominator n-1). The default will try to find the type used, by investigating which is most likely to result from integer valued datasets.
}
\item{baseline}{What baseline model should be used? \code{"stationary_random_intercept"} includes both within- and between person variances constrained equal across time (default), \code{"stationary"} only includes within-person variances constrained equal across time, \code{"independence"} (default up to version 0.11) includes a variance for every variable at every time point (not constrained equal across time), and \code{"none"} includes no baseline model. The baseline model always uses the same mean structure as the model itself (see the \code{stationary_means} argument and the details below).}
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
Optional alias for \code{within_latent}. If both \code{within} and \code{within_latent} are explicitly specified, \code{within_latent} takes precedence and a warning is issued.
}
\item{between}{
Optional alias for \code{between_latent}. If both \code{between} and \code{between_latent} are explicitly specified, \code{between_latent} takes precedence and a warning is issued.
}
\item{penalty_lambda}{Numeric penalty strength for penalized ML estimation (PML/PFIML). \code{NA} (default) triggers automatic selection via EBIC-based grid search when a penalized estimator is used; set to a specific numeric value to use a fixed penalty strength (0 = no penalty). See \code{\link{find_penalized_lambda}} and \code{\link{penalize}}.}
\item{penalty_alpha}{Elastic net mixing parameter: 1 = LASSO (default), 0 = ridge.}
\item{penalize_matrices}{Character vector of matrix names to penalize. If missing, defaults are selected based on the model type.}
\item{rho_zeta_within}{Only used when \code{within_latent = "cor"}. Correlation matrix specification (\code{"full"}, \code{"zero"}, or a pattern matrix; see \code{omega_zeta_within}).}
\item{SD_zeta_within}{Only used when \code{within_latent = "cor"}. Diagonal standard-deviation matrix specification (\code{"diag"}, \code{"full"}, or a pattern matrix; see \code{delta_zeta_within}).}
\item{rho_epsilon_within}{Only used when \code{within_residual = "cor"}. Correlation matrix specification (\code{"zero"}, \code{"full"}, or a pattern matrix).}
\item{SD_epsilon_within}{Only used when \code{within_residual = "cor"}. Diagonal standard-deviation matrix specification (\code{"diag"}, \code{"full"}, or a pattern matrix).}
\item{rho_zeta_between}{Only used when \code{between_latent = "cor"}. Correlation matrix specification (\code{"full"}, \code{"zero"}, or a pattern matrix).}
\item{SD_zeta_between}{Only used when \code{between_latent = "cor"}. Diagonal standard-deviation matrix specification (\code{"diag"}, \code{"full"}, or a pattern matrix).}
\item{rho_epsilon_between}{Only used when \code{between_residual = "cor"}. Correlation matrix specification (\code{"zero"}, \code{"full"}, or a pattern matrix).}
\item{SD_epsilon_between}{Only used when \code{between_residual = "cor"}. Diagonal standard-deviation matrix specification (\code{"diag"}, \code{"full"}, or a pattern matrix).}
\item{beta_epsilon}{Residual (within-person) temporal effects: a model matrix for the lag-1 autoregression of the residuals, \code{epsilon_t = beta_epsilon epsilon_(t-1) + u_t}. \code{"zero"} (default) gives white-noise residuals (the model as it was before this argument existed). \code{"diag"} gives each indicator its own residual autoregression, which models item-specific carry-over (autocorrelated uniquenesses) with a geometrically decaying structure; \code{"full"} or a pattern matrix are also allowed. When \code{beta_epsilon} is non-zero, \code{sigma_epsilon_within} (and its typed variants) parameterizes the \emph{innovation} covariance of the residual process, and the stationary residual covariance is available as the \code{"sigma_epsilon_within_stationary"} matrix in \code{\link{getmatrix}}. At least three waves are required for a non-zero \code{beta_epsilon} to be identified, and at least four are recommended (the residual autoregression, the stable between-person residual variance in \code{sigma_epsilon_between}, and the latent dynamics compete for the same lagged covariance information).}
\item{temporal_latent}{Parameterization of the latent temporal effects: \code{"raw"} (default) models \code{beta} directly; \code{"PDC"} models the latent partial directed correlations directly, with the \code{PDC} matrix encoding from = row and to = column (the transposed orientation relative to \code{beta}). A (pattern-constrained) PDC model is fit-equivalent to the corresponding raw model, but standard errors, confidence intervals, equality constraints and modification indices then apply directly on the PDC scale.}
\item{temporal_residual}{Parameterization of the residual temporal effects (\code{beta_epsilon}): \code{"raw"} (default) or \code{"PDC"} (models the residual partial directed correlations \code{PDC_epsilon}, standardized with respect to the residual innovation covariance, directly).}
\item{PDC}{Only used when \code{temporal_latent = "PDC"}. Model matrix for the latent partial directed correlations (from = row, to = column): \code{"full"}, \code{"zero"}, \code{"diag"}, or a pattern matrix.}
\item{PDC_epsilon}{Only used when \code{temporal_residual = "PDC"}. Model matrix for the residual partial directed correlations (from = row, to = column): \code{"zero"} (default), \code{"diag"}, \code{"full"}, or a pattern matrix.}
\item{stationary_means}{Logical. \code{TRUE} (default, and the behavior of all earlier versions) has item intercepts \code{nu} that are invariant across waves and latent means fixed to zero, so that the implied means are equal across waves. \code{FALSE} makes the item intercepts wave-specific (the \code{nu} matrix becomes \code{nVar} by \code{nTime}), which saturates the mean structure and adds \code{nVar * (nTime - 1)} parameters; the latent means stay fixed to zero. See the details below.}
\item{\dots}{
Arguments sent to \code{dlvm1}.
}
}
\details{
\bold{Stationary versus wave-specific means.} By default the model is fully stationary in the means: the item intercepts \code{nu} are invariant across waves and the latent means are fixed to zero, so the implied means are equal at every wave. With \code{stationary_means = FALSE} the item intercepts become wave-specific, \eqn{\mu_t = \nu_t + \Lambda \alpha} with \eqn{\alpha = 0}. This is the classical \emph{trend-stationary} case: the level shifts are period effects exogenous to the dynamics, the within-person latent deviation process remains covariance-stationary, and wave 1 is still drawn from that stationary distribution. The implied covariance matrix is therefore \emph{identical} to that of the stationary-means model; only the mean part of the model (and of its Jacobian) changes.

\bold{Identification.} The latent means (\code{mu_eta}) are fixed to zero with a stationary as well as with a wave-specific mean structure, and the wave-specific intercepts are free: one per variable--wave combination present in the \code{vars} design matrix (variable--wave combinations that are absent have no intercept parameter). Wave-specific intercepts and wave-specific latent means are not jointly identified, so one of the two has to carry the level shifts. Placing them in the intercepts saturates the mean structure, which is what keeps the baseline model (an observed-variable \code{\link{panelvar}} model, which cannot have latent means) nested in the fitted model, so that incremental fit indices remain interpretable. The estimated intercepts then equal the observed per-wave means, and the likelihood-ratio test of \code{stationary_means = TRUE} against \code{FALSE} is a test of equal means over time with \code{nVar * (nTime - 1)} degrees of freedom. The flip side is that a saturated mean structure cannot separate period effects on the constructs from mean-level measurement drift (scalar non-invariance), and cannot be tested for misfit: \code{nu} absorbs both. Users who want the more restrictive alternative (invariant intercepts with free latent means at waves 2 and later, \code{nLat(nTime - 1)} parameters, testable when \code{nLat < nVar}) can specify it manually with the \code{nu} and \code{mu_eta} arguments plus \code{identify = FALSE}, keeping in mind that the baseline model is then no longer nested in the fitted model. In the degenerate observed-variable case (\code{lambda} the identity matrix with zero residual variances) the model with \code{stationary_means = FALSE} is identical to \code{\link{panelvar}} with \code{stationary_means = FALSE}, with \code{nu} playing the role of \code{mu}.

\bold{No free innovation variances.} Only the means can be freed over time, not the innovation (co)variances: a wave-varying innovation covariance has no stationary distribution for wave 1 to be drawn from, so the stationarity machinery on which this family is built does not apply. Wave-specific innovation structures are covered by the RI-CLPM / RI-CLPN framework instead, see \code{\link{ri_clpm}}.

\bold{Baseline model.} The baseline model always uses the same mean structure as the fitted model: with \code{stationary_means = FALSE} the baseline (a \code{\link{panelvar}} model with empty networks) gets wave-specific means as well, so that both models have exactly the same (saturated) mean structure. A baseline model must be nested in the fitted model, otherwise it can fit better than the fitted model on the mean side (a free-means baseline against an equal-means target would beat it), which makes incremental fit indices meaningless or negative. A corollary is that incremental fit indices (CFI, TLI, NFI, ...) are \bold{not} directly comparable between a \code{stationary_means = TRUE} and a \code{stationary_means = FALSE} fit of the same data, because the two use different baselines by design; compare those two models with a likelihood-ratio test (\code{\link{compare}}) instead. The saturated model already has a saturated mean structure and is unchanged.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}

\author{
Sacha Epskamp
}

\seealso{
\code{\link{panelvar}} for the observed-variable (panel VAR) special case of this family.
}

\examples{
library("dplyr")

# Smoke data cov matrix, based on LISS data panel https://www.dataarchive.lissdata.nl
smoke <- structure(c(47.2361758611759, 43.5366809116809, 41.0057465682466, 
                     43.5366809116809, 57.9789886039886, 47.6992521367521, 
                     41.0057465682466, 
                     47.6992521367521, 53.0669434731935), dim = c(3L, 3L),
                   dimnames = list(
                       c("smoke2008", "smoke2009", "smoke2010"), c("smoke2008", 
                   "smoke2009", "smoke2010")))

# Design matrix:
design <- matrix(rownames(smoke),1,3)

# Form dlvm1 model with a single latent variable per wave:
mod <- dlvm1(vars = design, lambda = matrix(1),
             covs = smoke, nobs = 352
)

\donttest{
# Run model:
mod <- mod \%>\% runmodel

# Evaluate fit:
mod \%>\% fit
}
}
