\name{ml_var1}
\alias{ml_var1}
\alias{ml_gvar1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-level lag-1 vector auto-regression family (random intercepts)
}
\description{
This family models a multi-level (random intercept) lag-1 vector auto-regression
(VAR) / graphical VAR for intensive longitudinal data (many repeated measures nested
in subjects), estimated by a two-level summary-statistics pseudo maximum likelihood.
Each subject follows a within-person stationary VAR(1) process (temporal matrix
\code{beta} and contemporaneous structure \code{sigma_zeta_within}) around a
person-specific random intercept, and the random intercepts have a between-person
covariance structure (\code{sigma_zeta_between}). The \code{ml_var1} function models
the contemporaneous and between-person structures as variance--covariance matrices by
default, and the \code{ml_gvar1} wrapper models them as Gaussian graphical models
(GGMs), yielding the three multi-level GVAR networks: the temporal network
(\code{beta} / partial directed correlations), the contemporaneous network
(\code{omega_zeta_within}) and the between-person network (\code{omega_zeta_between}).

This model targets the same population model as \code{\link{panelvar}} /
\code{panelgvar} (a multi-level VAR with random intercepts), but is designed for
intensive longitudinal / experience-sampling (ESM) data, where the number of repeated
measures per subject is large and the full-likelihood \code{panelvar} route is
infeasible. See Details for the estimator and its statistical caveats.

\strong{This estimator is experimental} (a message is shown once per session).
}
\usage{
ml_var1(data, idvar, vars, dayvar, beepvar, groups, groupvar,
        within_latent = c("cov", "chol", "prec", "ggm", "cor"),
        between_latent = c("cov", "chol", "prec", "ggm", "cor"),
        beta = "full",
        omega_zeta_within = "full", delta_zeta_within = "diag",
        kappa_zeta_within = "full", sigma_zeta_within = "full",
        lowertri_zeta_within = "full", rho_zeta_within = "full",
        SD_zeta_within = "full",
        omega_zeta_between = "full", delta_zeta_between = "diag",
        kappa_zeta_between = "full", sigma_zeta_between = "full",
        lowertri_zeta_between = "full", rho_zeta_between = "full",
        SD_zeta_between = "full",
        mu, equal = "none", baseline_saturated = TRUE, estimator = "ML",
        optimizer, storedata = FALSE,
        standardize = c("none", "z", "quantile"),
        sampleStats, verbose = FALSE, toeplitz = TRUE)

ml_gvar1(..., within_latent = "ggm", between_latent = "ggm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the (long-format) data used in the analysis. Must be a raw dataset with one row per measurement occasion, containing the subject id variable (\code{idvar}) and the variables in \code{vars}.
}
  \item{idvar}{
Required string indicating the subject/cluster ID variable in \code{data}. Each subject's rows form one within-person time series; the random intercepts vary over subjects. (Use \code{\link{var1}} for single-subject time-series models.)
}
  \item{vars}{
An optional character vector encoding the variables used in the analysis. Must equal names in \code{data}. If missing, all columns except \code{idvar}, \code{dayvar}, \code{beepvar} and the group variable are used.
}
  \item{dayvar}{
Optional string indicating assessment day. Adding this argument makes sure that the first measurement of a day is not regressed on the last measurement of the previous day. IMPORTANT: only add this if the data has multiple observations per day.
}
  \item{beepvar}{
Optional string indicating assessment beep per day. Adding this argument will cause non-consecutive beeps to be treated as missing (the corresponding lag pairs are dropped).
}
  \item{groups}{
Deprecated. Use \code{groupvar} instead. An optional string indicating the name of the group variable in \code{data}.
}
  \item{groupvar}{
An optional string indicating the name of the group variable in \code{data} for multi-group analysis. Replaces the deprecated \code{groups} argument.
}
  \item{within_latent}{
The type of within-person contemporaneous model to be used: \code{"cov"} (covariance, default), \code{"chol"} (Cholesky), \code{"prec"} (precision), \code{"ggm"} (Gaussian graphical model, giving the contemporaneous network \code{omega_zeta_within}), or \code{"cor"} (correlation).
}
  \item{between_latent}{
The type of between-person (random-intercept) model to be used: \code{"cov"} (default), \code{"chol"}, \code{"prec"}, \code{"ggm"} (giving the between-person network \code{omega_zeta_between}), or \code{"cor"}.
}
  \item{beta}{
A model matrix encoding the temporal relationships (transpose of the temporal network; from = row, to = column, as in \code{\link{var1}}). A 0 encodes a fixed to zero element, a 1 a free to estimate element, and higher integers equality constraints. For multiple groups, this can be a list or array. Can also be \code{"full"} (default) for a full temporal network or \code{"zero"} for an empty temporal network.
}
  \item{omega_zeta_within}{
Only used when \code{within_latent = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constraints. For multiple groups, this can be a list or array.
}
  \item{delta_zeta_within}{
Only used when \code{within_latent = "ggm"}. Scaling (delta) matrix. Can be \code{"diag"}, \code{"zero"} (not recommended), or a model matrix.
}
  \item{kappa_zeta_within}{
Only used when \code{within_latent = "prec"}. Precision matrix. Can be \code{"full"}, \code{"diag"}, or a model matrix.
}
  \item{sigma_zeta_within}{
Only used when \code{within_latent = "cov"}. Covariance matrix. Can be \code{"full"}, \code{"diag"}, or a model matrix.
}
  \item{lowertri_zeta_within}{
Only used when \code{within_latent = "chol"}. Cholesky (lower triangular) matrix. Can be \code{"full"}, \code{"diag"}, or a model matrix.
}
  \item{rho_zeta_within}{
Only used when \code{within_latent = "cor"}. Correlation matrix specification (\code{"full"}, \code{"zero"}, or a pattern matrix).
}
  \item{SD_zeta_within}{
Only used when \code{within_latent = "cor"}. Diagonal standard-deviation matrix specification (\code{"full"}, \code{"diag"}, or a pattern matrix).
}
  \item{omega_zeta_between}{
Only used when \code{between_latent = "ggm"}. Between-person network. Can be \code{"full"}, \code{"zero"}, or a model matrix.
}
  \item{delta_zeta_between}{
Only used when \code{between_latent = "ggm"}. Scaling (delta) matrix. Can be \code{"diag"}, \code{"zero"} (not recommended), or a model matrix.
}
  \item{kappa_zeta_between}{
Only used when \code{between_latent = "prec"}. Precision matrix. Can be \code{"full"}, \code{"diag"}, or a model matrix.
}
  \item{sigma_zeta_between}{
Only used when \code{between_latent = "cov"}. Covariance matrix. Can be \code{"full"}, \code{"diag"}, or a model matrix.
}
  \item{lowertri_zeta_between}{
Only used when \code{between_latent = "chol"}. Cholesky (lower triangular) matrix. Can be \code{"full"}, \code{"diag"}, or a model matrix.
}
  \item{rho_zeta_between}{
Only used when \code{between_latent = "cor"}. Correlation matrix specification (\code{"full"}, \code{"zero"}, or a pattern matrix).
}
  \item{SD_zeta_between}{
Only used when \code{between_latent = "cor"}. Diagonal standard-deviation matrix specification (\code{"full"}, \code{"diag"}, or a pattern matrix).
}
  \item{mu}{
Optional vector encoding the stationary means. Set elements to 0 to indicate fixed to zero constraints, 1 to indicate free means, and higher integers to indicate equality constraints. For multiple groups, this can be a list or matrix.
}
  \item{equal}{
A character vector indicating which matrices should be constrained equal across groups (e.g. \code{"beta"}).
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT recommended to be used manually.
}
  \item{estimator}{
The estimator to be used. Only \code{"ML"} (the two-level summary-statistics pseudo maximum likelihood) is supported in this version. For full-information maximum likelihood on the same estimand use \code{\link{panelvar}} / \code{\link{dlvm1}}.
}
  \item{optimizer}{
The optimizer to be used. Can be one of \code{"nlminb"} (the default R \code{nlminb} function), \code{"ucminf"}, \code{"nloptr_TNEWTON"}, and \code{"LBFGS++"}.
}
  \item{storedata}{
Logical, should the raw data be stored?
}
  \item{standardize}{
Which standardization method should be used on the raw variables before the lag embedding? \code{"none"} (default), \code{"z"} for z-scores, and \code{"quantile"} for a non-parametric quantile transformation. Standardization is grand (over all observations), not within-person, so that the between-person variance survives.
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally (for the recursive baseline/saturated calls).
}
  \item{verbose}{
Logical, should progress be printed to the console?
}
  \item{toeplitz}{
Internal. If \code{TRUE} (default) the stationary Toeplitz VAR(1) within structure is imposed; \code{FALSE} selects the free-within (saturated) mode used internally to build the reference model. Not intended for direct use.
}
\item{\dots}{
Arguments sent to \code{ml_var1} (used by \code{ml_gvar1}).
}
}
\details{
For each subject the observed \eqn{p}-variate series is augmented into complete lag
pairs \eqn{z_{it} = (y_{i,t-1}', y_{it}')'} (the \eqn{2p}-variate "Toeplitz" embedding
also used by \code{\link{var1}}; pairs are broken across days/beeps via \code{dayvar} /
\code{beepvar}). These lag pairs are then treated, \emph{within subject}, as
conditionally independent given the subject's random intercept \eqn{b_i}:
\deqn{z_{it} \mid b_i \sim N\!\left(\mu_z + (b_i', b_i')',\ \Sigma_W\right), \qquad b_i \sim N(0, \Sigma_B),}{z_it | b_i ~ N(mu_z + (b_i', b_i')', Sigma_W), b_i ~ N(0, Sigma_B),}
with \eqn{\mu_z = (\mu', \mu')'}{mu_z = (mu', mu')'} the stationary mean repeated across
the two halves, the within-person block a stationary Toeplitz structure
\eqn{\Sigma_W = [[\Sigma_0, \Sigma_0 \beta'], [\beta \Sigma_0, \Sigma_0]]}{Sigma_W =
[[Sigma0, Sigma0 beta'], [beta Sigma0, Sigma0]]} with
\eqn{\mathrm{vec}(\Sigma_0) = (I - \beta \otimes \beta)^{-1} \mathrm{vec}(\Sigma_\zeta)}{vec(Sigma0) = (I - beta (x) beta)^-1 vec(Sigma_zeta)},
and the between-person block \eqn{\Sigma_{Bf} = 1_{2 \times 2} \otimes \Sigma_B}{Sigma_Bf
= 1_(2x2) (x) Sigma_B} (a rank-\eqn{p}, \eqn{2p \times 2p} matrix). Minus twice the
log-likelihood is the two-level Gaussian objective of McDonald & Goldstein (1989) /
Muthen (1990) evaluated at \eqn{(\mu_z, \Sigma_W, \Sigma_{Bf})}{(mu_z, Sigma_W,
Sigma_Bf)} with the subject as cluster.

Because consecutive lag pairs share \eqn{y_{it}}{y_it}, this is a \emph{working-independence
pseudo (composite) likelihood}, exactly the move \code{\link{var1}} makes at a single
subject. Its key advantage is that the data are compressed \emph{once} into the two-level
sufficient statistics of the lag pairs, so every fit / gradient / Fisher-information
evaluation works only with \eqn{2p \times 2p} matrices, and the per-evaluation cost is
\emph{independent of the number of time points}. \code{ml_var1} therefore uses the same
data augmentation as \code{\link{var1}} and the same two-level estimator as
\code{\link{ml_lvm}} (\code{estimator = "ML"}).

\code{ml_var1} shares the estimand with \code{\link{panelvar}} / \code{panelgvar} (a
multi-level VAR(1) with random intercepts) but not the estimator: \code{panelvar} forms
the full \eqn{(T p) \times (T p)} block-Toeplitz likelihood, which is statistically ideal
for panel data (few waves, many subjects) but grows with the number of waves \eqn{T} and
becomes infeasible for intensive longitudinal / ESM data with dozens to hundreds of
measurements per subject. \emph{Use \code{panelvar} / \code{panelgvar} for few-wave panel
data, and \code{ml_var1} / \code{ml_gvar1} for many repeated measures where \code{panelvar}
is infeasible.} In this first version, incomplete lag pairs (rows with any missing value,
including the structural lag rows at day/beep breaks) are dropped (complete-pairs, listwise
analysis); a full-information treatment of the augmented rows is planned for a later version.

\bold{Statistical caveats.} Because the estimator treats consecutive lag pairs as
independent (working-independence pseudo-likelihood, as in \code{\link{var1}}), two caveats
apply. First, under positive temporal autocorrelation the between-person covariance
\eqn{\Sigma_B}{Sigma_B} carries a small upward finite-length bias of order \eqn{O(1/T)}: in
simulation (\eqn{p = 4}, \eqn{J = 60}, \eqn{\beta}{beta} spectral radius 0.5) the
between-person variances were biased about +9\% at \eqn{T = 25} and +3\% at \eqn{T = 50},
and negligibly by \eqn{T = 100}; a full-likelihood \code{\link{panelvar}} / \code{panelgvar}
fit does not show this bias. Second, the Fisher-based standard errors are pseudo-likelihood
standard errors and can be anti-conservative -- most noticeably for the contemporaneous
(within) network, where reported standard errors were about 75--80\% of the empirical
standard deviation, an under-statement that grows with temporal dependence and does not
vanish as \eqn{T} grows. Point estimates of the temporal network \code{beta} and the
contemporaneous network \code{omega_zeta_within} are consistent and closely track the
full-likelihood \code{panelvar} estimates.

The saturated (unrestricted) and baseline (independence) reference models used by
\code{\link{runmodel}} for the fit measures are estimated within the \code{ml_var1}
framework itself (a free two-level model over the augmented \eqn{2p} variables and a
two-level independence model, respectively). The number of independent observations used in
sample-size-dependent fit measures (BIC, RMSEA) is the number of \emph{subjects}. The model
layer is currently implemented in R only (no C++ acceleration), but estimation is still fast
because every evaluation touches only \eqn{2p \times 2p} matrices.
}

\value{
An object of the class psychonetrics (\link{psychonetrics-class}).
}

\references{
McDonald, R. P., & Goldstein, H. (1989). Balanced versus unbalanced designs for linear structural relations in two-level data. British Journal of Mathematical and Statistical Psychology, 42(2), 215-232.

Muthen, B. O. (1990). Mean and covariance structure analysis of hierarchical data. UCLA Statistics Series #62.

Epskamp, S. (2020). Psychometric network models from time-series and panel data. Psychometrika, 85(1), 206-231. \doi{10.1007/s11336-020-09697-3}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\seealso{
\code{\link{var1}} for the single-subject (or pooled) VAR(1)/GVAR, \code{\link{panelvar}} / \code{panelgvar} for the full-likelihood random-intercept VAR on panel data, \code{\link{ml_lvm}} for the two-level latent variable model that supplies the estimator, and \code{\link{dlvm1}} for the dynamic latent variable model.
}

\examples{
\dontrun{
library("dplyr")

# Simulate a tiny two-level VAR(1) data set (random intercepts):
set.seed(1)
p <- 3        # nodes
J <- 40       # subjects
TT <- 25      # measurements per subject

# Temporal (beta) and contemporaneous (within) structure:
beta <- diag(0.3, p); beta[1, 2] <- 0.2; beta[3, 1] <- -0.15
SigZ <- diag(p) * 0.8; SigZ[1, 2] <- SigZ[2, 1] <- 0.2
Sig0 <- matrix(solve(diag(p^2) - kronecker(beta, beta), as.vector(SigZ)), p, p)
SigB <- diag(p) * 0.4; SigB[2, 3] <- SigB[3, 2] <- 0.1

# Generate:
long <- do.call(rbind, lapply(seq_len(J), function(i){
  bi  <- MASS::mvrnorm(1, rep(0, p), SigB)          # random intercept
  eta <- MASS::mvrnorm(1, rep(0, p), Sig0)
  Y <- matrix(NA_real_, TT, p)
  Y[1, ] <- eta + bi
  for (t in 2:TT){
    eta <- as.vector(beta \%*\% eta) + MASS::mvrnorm(1, rep(0, p), SigZ)
    Y[t, ] <- eta + bi
  }
  data.frame(id = i, beep = seq_len(TT), Y)
}))
names(long)[3:5] <- paste0("V", 1:p)

# Fit the multi-level graphical VAR (GGM within and between):
mod <- ml_gvar1(long, vars = paste0("V", 1:p), idvar = "id", beepvar = "beep")
mod <- runmodel(mod)

# The three multi-level GVAR networks:
getmatrix(mod, "beta")                # temporal network
getmatrix(mod, "omega_zeta_within")   # contemporaneous network
getmatrix(mod, "omega_zeta_between")  # between-person network
}
}
