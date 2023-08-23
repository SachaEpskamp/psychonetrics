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
                    baseline_saturated = TRUE, estimator = c("FIML",
                    "MUML"), optimizer, storedata = FALSE, verbose =
                    FALSE, standardize = c("none", "z", "quantile"),
                    sampleStats)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Must be a raw dataset.
}
  \item{lambda}{
A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix. Could also be the result of \code{\link{simplestructure}}.
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
A model matrix encoding the within-cluster structural. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix. Defaults to \code{"zero"}.
}
  \item{beta_between}{
A model matrix encoding the between-cluster structural. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix. Defaults to \code{"zero"}.
}
  \item{omega_zeta_within}{
Only used when \code{within_latent = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta_within}{
Only used when \code{within_latent = "ggm"}. Can be \code{"diag"}, \code{"zero"} (not recommended), or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta_within}{
Only used when \code{within_latent = "prec"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_zeta_within}{
Only used when \code{within_latent = "cov"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta_within}{
Only used when \code{within_latent = "chol"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_epsilon_within}{
Only used when \code{within_residual = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon_within}{
Only used when \code{within_residual = "ggm"}. Can be \code{"diag"}, \code{"zero"} (not recommended), or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon_within}{
Only used when \code{within_residual = "prec"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon_within}{
Only used when \code{within_residual = "cov"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon_within}{
Only used when \code{within_residual = "chol"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_zeta_between}{
Only used when \code{between_latent = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta_between}{
Only used when \code{between_latent = "ggm"}. Can be \code{"diag"}, \code{"zero"} (not recommended), or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta_between}{
Only used when \code{between_latent = "prec"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_zeta_between}{
Only used when \code{between_latent = "cov"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta_between}{
Only used when \code{between_latent = "chol"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_epsilon_between}{
Only used when \code{between_residual = "ggm"}. Can be \code{"full"}, \code{"zero"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon_between}{
Only used when \code{between_residual = "ggm"}. Can be \code{"diag"}, \code{"zero"} (not recommended), or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon_between}{
Only used when \code{between_residual = "prec"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon_between}{
Only used when \code{between_residual = "cov"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon_between}{
Only used when \code{between_residual = "chol"}. Can be \code{"full"}, \code{"diag"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{nu}{
Optional vector encoding the intercepts of the observed variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{nu_eta}{
Optional vector encoding the intercepts of the latent variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
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
Estimator used. Currently only \code{"FIML"} is supported.
}
 \item{optimizer}{
The optimizer to be used. Usually either \code{"nlminb"} (with box constrains) or \code{"ucminf"} (ignoring box constrains), but any optimizer supported by \code{optimr} can be used.
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
\item{\dots}{
Arguments sent to 'ml_lvm'
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

