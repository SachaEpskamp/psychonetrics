\name{ml_lvm}
\alias{ml_lvm}
\alias{ml_lnm}
\alias{ml_rnm}
\alias{ml_lrnm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-level latent variable model family
}
\description{
This family is the two-level random intercept variant of the \code{\link{lvm}} model family. It is mostly a special case of the \code{\link{dlvm1}} family, with the addition of structural effects rather than temporal effects in the \code{beta} matrix.
}
\usage{
ml_lnm(...)
ml_rnm(...)
ml_lrnm(...)
ml_lvm(data, lambda, clusters, within_latent = c("cov",
                   "chol", "prec", "ggm"), within_residual = c("cov",
                   "chol", "prec", "ggm"), between_latent = c("cov",
                   "chol", "prec", "ggm"), between_residual = c("cov",
                   "chol", "prec", "ggm"), beta_within = "zero",
                   beta_between = "zero", omega_zeta_within = "full",
                   delta_zeta_within = "full", kappa_zeta_within =
                   "full", sigma_zeta_within = "full",
                   lowertri_zeta_within = "full", omega_epsilon_within =
                   "zero", delta_epsilon_within = "diag",
                   kappa_epsilon_within = "diag", sigma_epsilon_within =
                   "diag", lowertri_epsilon_within = "diag",
                   omega_zeta_between = "full", delta_zeta_between =
                   "full", kappa_zeta_between = "full",
                   sigma_zeta_between = "full", lowertri_zeta_between =
                   "full", omega_epsilon_between = "zero",
                   delta_epsilon_between = "diag", kappa_epsilon_between
                   = "diag", sigma_epsilon_between = "diag",
                   lowertri_epsilon_between = "diag", nu, nu_eta,
                   identify = TRUE, identification = c("loadings",
                   "variance"), vars, latents, groups, equal = "none",
                   baseline_saturated = TRUE, estimator = c("default",
                   "FIML", "ML"),
                   optimizer, storedata = FALSE, verbose =
                   FALSE, standardize = c("none", "z", "quantile"),
                   sampleStats, bootstrap = FALSE, boot_sub,
                   boot_resample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Must be a raw dataset.
}
  \item{lambda}{
A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix. Could also be the result of \code{\link{simplestructure}}.
}
  \item{clusters}{
A string indicating the variable in the dataset that describes group membership.
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
  \item{beta_within}{
A model matrix encoding the within-cluster structural. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix. Defaults to \code{"zero"}.
}
  \item{beta_between}{
A model matrix encoding the between-cluster structural. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix. Defaults to \code{"zero"}.
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
Optional vector encoding the intercepts of the observed variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{nu_eta}{
Optional vector encoding the intercepts of the latent variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
 \item{identify}{
Logical, should the model be automatically identified?
}
  \item{identification}{
Type of identification used. \code{"loadings"} to fix the first factor loadings to 1, and \code{"variance"} to fix the diagonal of the latent variable model matrix (sigma_zeta, lowertri_zeta, delta_zeta or kappa_zeta) to 1.
}
  \item{vars}{
An optional character vector with names of the variables used.
}
  \item{latents}{
An optional character vector with names of the latent variables.
}
  \item{groups}{
An optional string indicating the name of the group variable in \code{data}.
}
  \item{equal}{
A character vector indicating which matrices should be constrained equal across groups. 
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
  \item{estimator}{
The estimator to be used. \code{"ML"} and \code{"FIML"} maximize the \emph{same} (exact) two-level maximum-likelihood objective and converge to the same solution; they differ only in how that likelihood is \emph{computed}, exactly as for the single-level model families in this package, where \code{"ML"} fits the model to summary statistics and \code{"FIML"} fits it row by row to the raw data. \code{"ML"} (experimental) evaluates the likelihood from two-level sufficient statistics (pooled within-cluster covariance and per-cluster-size between moments; McDonald & Goldstein, 1989), which for complete data is much faster and does not slow down as clusters grow. \code{"FIML"} (wide-format full-information maximum likelihood) stacks each cluster's members into one wide observation; it was the only estimator before version 0.15.31. Both handle within-cluster missing data (\code{"ML"} through a per-pattern, per-cluster likelihood; see Details). This is NOT the approximate "MUML" / pseudo-balanced estimator: the sufficient-statistics evaluation is exact for unbalanced designs, and the two estimators' fit-function values coincide to numerical precision (verified to ~1e-13). \code{"default"} selects \code{"ML"} when the data are complete and the largest cluster contains more than 5 units, and \code{"FIML"} otherwise; with missing data the default stays on the long-standing \code{"FIML"} path, so the (R-only) missing-data \code{"ML"} estimator is used only when \code{estimator = "ML"} is requested explicitly.
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
  \item{standardize}{
Which standardization method should be used? \code{"none"} (default) for no standardization, \code{"z"} for z-scores, and \code{"quantile"} for a non-parametric transformation to the quantiles of the marginal standard normal distribution.
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally. 
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
Arguments sent to 'ml_lvm'
}
}
\details{
The \code{"ML"} and \code{"FIML"} estimators are two computational routes to the \emph{same} exact two-level maximum-likelihood estimate, not two different statistical methods and not an exact-versus-approximate pair: they maximize the same log-likelihood and (up to optimizer tolerance) return identical estimates, standard errors and fit. The naming follows the convention used throughout psychonetrics -- \code{"ML"} fits to summary statistics while \code{"FIML"} fits row by row to the raw data -- and does not by itself indicate complete versus missing data, since both routes can accommodate within-cluster missingness. Choose \code{"ML"} for speed on (near-)complete data with non-trivial cluster sizes and \code{"FIML"} (or \code{"default"}) otherwise.

With \code{estimator = "FIML"}, every cluster is treated as one wide-format observation (members stacked side by side) and the model is estimated with full-information maximum likelihood over the missingness patterns. This handles missing data, but its cost grows quickly with the cluster size.

With \code{estimator = "ML"} (experimental), the same two-level random-intercept likelihood is instead evaluated from sufficient statistics (McDonald & Goldstein, 1989; Muthen, 1990) when the data are complete. Writing \eqn{\Sigma_W}{Sigma_W} and \eqn{\Sigma_B}{Sigma_B} for the implied within- and between-cluster covariance matrices and \eqn{\mu}{mu} for the implied mean vector, minus twice the log-likelihood equals
\deqn{(N - J)\left[\ln|\Sigma_W| + \mathrm{tr}(\Sigma_W^{-1} S_{PW})\right] + \sum_s m_s \left[\ln|\Sigma_s| + n_s \mathrm{tr}(\Sigma_s^{-1} A_s)\right] + N p \ln(2\pi),}{(N - J)[ln|Sigma_W| + tr(Sigma_W^-1 S_PW)] + sum_s m_s [ln|Sigma_s| + n_s tr(Sigma_s^-1 A_s)] + N p ln(2 pi),}
in which \eqn{N} is the number of units, \eqn{J} the number of clusters, \eqn{S_{PW}}{S_PW} the pooled within-cluster covariance matrix, and, per distinct cluster size \eqn{n_s} with \eqn{m_s} clusters, \eqn{\Sigma_s = \Sigma_W + n_s \Sigma_B}{Sigma_s = Sigma_W + n_s Sigma_B} and \eqn{A_s} the second-moment matrix of the corresponding cluster means around \eqn{\mu}{mu}. The two estimators optimize numerically identical objectives (the fit function values coincide at identical parameter values), so they converge to the same solution; the same likelihood is also used by \code{lavaan} for two-level models (Rosseel, 2012). The sufficient statistics are computed once, making \code{"ML"} much faster than \code{"FIML"} when clusters are large, and its cost does not depend on the cluster sizes. For complete data, standard errors and modification indices for \code{"ML"} use the analytic EXPECTED Fisher information. Note that \code{lavaan} uses the OBSERVED information by default for multilevel models, so standard errors reported by \code{lavaan} typically differ by a few percent; with \code{lavaan::sem(..., information = "expected")} the standard errors agree.

When the data contain within-cluster missing values (assumed missing at random), \code{estimator = "ML"} switches to a per-pattern, per-cluster evaluation of the same two-level likelihood: for each missingness pattern the within-cluster precision and log-determinant are obtained by a symmetric (Schur) inverse update of \eqn{\Sigma_W^{-1}}{Sigma_W^-1}, and these are accumulated per cluster into the quantities entering minus twice the log-likelihood,
\deqn{\sum_P f_P \ln|\Sigma_W^{(P)}| + \sum_j \ln|I + \Sigma_B A_j| + \sum_i (y_i - \mu)'\Sigma_W^{(P_i)-1}(y_i - \mu) - \sum_j p_j'(I + \Sigma_B A_j)^{-1}\Sigma_B\, p_j + n_{\mathrm{obs}} \ln(2\pi),}{sum_P f_P ln|Sigma_W^(P)| + sum_j ln|I + Sigma_B A_j| + sum_i (y_i - mu)' Sigma_W^(P_i)^-1 (y_i - mu) - sum_j p_j' (I + Sigma_B A_j)^-1 Sigma_B p_j + n_obs ln(2 pi),}
with \eqn{A_j = \sum_{i \in j}\Sigma_W^{(P_i)-1}}{A_j = sum_i Sigma_W^(P_i)^-1}, \eqn{p_j = \sum_{i \in j}\Sigma_W^{(P_i)-1}(y_i - \mu)}{p_j = sum_i Sigma_W^(P_i)^-1 (y_i - mu)} (observed coordinates only), \eqn{f_P} the number of units in pattern \eqn{P}, and \eqn{n_{\mathrm{obs}}}{n_obs} the number of observed values. This is the same two-level missing-data likelihood that \code{lavaan::sem(..., cluster = , missing = "ml")} maximizes (Rosseel, 2012), and the two agree to optimizer tolerance. The missing-data path is evaluated in R (no C++ acceleration) and reduces exactly to the sufficient-statistics likelihood above when no values are missing; its standard errors use the numeric (observed-information-based) Fisher information.

The saturated (unrestricted) and baseline (independence) reference models used by \code{\link{runmodel}} for the fit measures are themselves estimated as two-level models: the saturated model freely estimates the means and the full within- and between-cluster covariance matrices, and the baseline model only estimates means and the variances at both levels. These are the same reference models \code{lavaan} uses for two-level SEM, so the saturated and baseline log-likelihoods, the model and baseline chi-square, the degrees of freedom, the CFI and the TLI agree with \code{lavaan} (up to optimization tolerances), provided the \code{lavaan} model uses the same parameterization (\code{ml_lvm} shares the factor loadings across the within and between levels, so the matching \code{lavaan} model must constrain its level-1 and level-2 loadings to be equal). One convention difference remains: psychonetrics evaluates sample-size-dependent fit measures, such as the RMSEA (with its confidence interval) and the BIC, using the number of independent observations, which in these models is the number of \emph{clusters} \eqn{J}, whereas \code{lavaan} uses the total number of level-1 \emph{units} \eqn{N}. Because the RMSEA equals \eqn{\sqrt{\max((\chi^2/n)/\mathit{df} - 1/n,\, 0)}}{sqrt(max((chisq/n)/df - 1/n, 0))} with \eqn{n}{n} the number of independent observations, and the chi-square and degrees of freedom otherwise agree, the RMSEA reported by psychonetrics is larger than the value reported by \code{lavaan} by exactly a factor \eqn{\sqrt{N/J}}{sqrt(N/J)}. This applies to both estimators of the \code{ml_lvm} family.
}

\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}

\references{
McDonald, R. P., & Goldstein, H. (1989). Balanced versus unbalanced designs for linear structural relations in two-level data. British Journal of Mathematical and Statistical Psychology, 42(2), 215-232.

Muthen, B. O. (1990). Mean and covariance structure analysis of hierarchical data. UCLA Statistics Series #62.

Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. \doi{10.18637/jss.v048.i02}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

\examples{
\donttest{
# Simulate two-level data with 2 latent variables on both levels:
set.seed(1)
J <- 50 # clusters
nj <- 10 # units per cluster
p <- 6
Lambda <- matrix(0, p, 2)
Lambda[1:3,1] <- c(1, 0.8, 0.7)
Lambda[4:6,2] <- c(1, 0.9, 0.6)
PsiW <- matrix(c(1,0.3,0.3,1), 2)
PsiB <- matrix(c(0.5,0.2,0.2,0.4), 2)
data <- do.call(rbind, lapply(seq_len(J), function(j){
  bj <- MASS::mvrnorm(1, rep(0,2), PsiB)
  ej <- rnorm(p, sd = 0.4) # between-level residuals
  t(replicate(nj, as.vector(Lambda \%*\% (bj + MASS::mvrnorm(1, rep(0,2), PsiW))) +
      ej + rnorm(p, sd = 0.7)))
}))
data <- as.data.frame(data)
names(data) <- paste0("y", 1:p)
data$cluster <- rep(seq_len(J), each = nj)

# Lambda model matrix (shared across levels):
lambda <- matrix(0, p, 2)
lambda[1:3,1] <- lambda[4:6,2] <- 1

# Two-level ML estimator (fast, complete data):
mod <- ml_lvm(data, lambda = lambda, clusters = "cluster", estimator = "ML")
mod <- runmodel(mod)
mod

# The FIML estimator (slower; supports missing data) converges to the
# same solution:
mod_fiml <- ml_lvm(data, lambda = lambda, clusters = "cluster",
                   estimator = "FIML")
mod_fiml <- runmodel(mod_fiml)
}
}

