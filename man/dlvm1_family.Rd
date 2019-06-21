\name{dlvm1}
\alias{dlvm1}
\alias{panelgvar}
\alias{panelvar}
\alias{panel_lvgvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lag-1 dynamic latent variable model family of psychonetrics models for panel data
}
\description{
This is the family of models that models a dynamic factor model on panel data. There are four covariance structures that can be modeled in different ways: \code{within_latent}, \code{between_latent} for the within-person and between-person latent (contemporaneous) models respectively, and \code{within_residual}, \code{between_residual} for the within-person and between-person residual models respectively. The \code{panelgvar} wrapper function sets the \code{lambda} to an identity matrix, all residual variances to zero, and models within-person and between-person latent (contemporaneous) models as GGMs. The \code{panelvar} wrapper does the same but models contemporaneous relations as a variance-covariance matrix. Finally, the \code{panel_lvgvar} wrapper automatically models all latent networks as GGMs.
}
\usage{
dlvm1(data, vars, lambda, within_latent = c("cov", "chol", "prec", "ggm"), 
      within_residual = c("cov", "chol", "prec", "ggm"), between_latent = c("cov", 
      "chol", "prec", "ggm"), between_residual = c("cov", "chol", "prec", "ggm"), 
      beta = "full", omega_zeta_within = "full", delta_zeta_within = "full", 
      kappa_zeta_within = "full", sigma_zeta_within = "full", lowertri_zeta_within = 
      "full", omega_epsilon_within = "empty", delta_epsilon_within = "empty", 
      kappa_epsilon_within = "empty", sigma_epsilon_within = "empty", 
      lowertri_epsilon_within = "empty", omega_zeta_between = "full", 
      delta_zeta_between = "full", kappa_zeta_between = "full", sigma_zeta_between = 
      "full", lowertri_zeta_between = "full", omega_epsilon_between = "empty", 
      delta_epsilon_between = "empty", kappa_epsilon_between = "empty", 
      sigma_epsilon_between = "empty", lowertri_epsilon_between = "empty", tau, mu_eta, 
      identify = TRUE, identification = c("loadings", "variance"), latents, groups, covs, 
      means, nobs, missing = "listwise", equal = "none", baseline_saturated = TRUE, 
      estimator = "ML", optimizer = "default", storedata = FALSE, sampleStats)

panelgvar(data, vars, ...)

panelvar(data, vars, ...)

panel_lvgvar(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{vars}{ Required argument.
Different from in other psychonetrics models, this must be a *matrix* with each row indicating a variable and each column indicating a measurement. The matrix must be filled with names of the variables in the dataset corresponding to variable i at wave j. NAs can be used to indicate missing waves. The rownames of this matrix will be used as variable names.
}
  \item{lambda}{ Required argument.
A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
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
A model matrix encoding the temporal relationships (transpose of temporal network). A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix. Can also be \code{"full"} for a full temporal network or \code{"empty"} for an empty temporal network.
}
  \item{omega_zeta_within}{
Only used when \code{within_latent = "ggm"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta_within}{
Only used when \code{within_latent = "ggm"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta_within}{
Only used when \code{within_latent = "prec"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_zeta_within}{
Only used when \code{within_latent = "cov"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta_within}{
Only used when \code{within_latent = "chol"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_epsilon_within}{
Only used when \code{within_residual = "ggm"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon_within}{
Only used when \code{within_residual = "ggm"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon_within}{
Only used when \code{within_residual = "prec"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon_within}{
Only used when \code{within_residual = "cov"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon_within}{
Only used when \code{within_residual = "chol"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_zeta_between}{
Only used when \code{between_latent = "ggm"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta_between}{
Only used when \code{between_latent = "ggm"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta_between}{
Only used when \code{between_latent = "prec"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_zeta_between}{
Only used when \code{between_latent = "cov"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta_between}{
Only used when \code{between_latent = "chol"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_epsilon_between}{
Only used when \code{between_residual = "ggm"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon_between}{
Only used when \code{between_residual = "ggm"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon_between}{
Only used when \code{between_residual = "prec"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon_between}{
Only used when \code{between_residual = "cov"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon_between}{
Only used when \code{between_residual = "chol"}. Can be \code{"full"}, \code{"empty"}, or a typical model matrix with 0s indicating parameters constrained to zero, 1s indicating free parameters, and higher integers indicating equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{tau}{
Optional vector encoding the intercepts of the observed variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{mu_eta}{
Optional vector encoding the means of the latent variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
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
An optional string indicating the name of the group variable in \code{data}.
}
  \item{covs}{
A sample variance--covariance matrix, or a list/array of such matrices for multiple groups. IMPORTANT NOTE: psychonetrics expects the maximum likelihood (ML) covariance matrix, which is NOT obtained from \code{\link{cov}} directly. Manually rescale the result of \code{\link{cov}} with \code{(nobs - 1)/nobs} to obtain the ML covariance matrix.
}
  \item{means}{
A vector of sample means, or a list/matrix containing such vectors for multiple groups.
}
  \item{nobs}{
The number of observations used in \code{covs} and \code{means}, or a vector of such numbers of observations for multiple groups.
}
  \item{missing}{
How should missingness be handled in computing the sample covariances and number of observations when \code{data} is used. Can be \code{"listwise"} for listwise deletion, or \code{"pairwise"} for pairwise deletion.
}
  \item{equal}{
A character vector indicating which matrices should be constrained equal across groups. 
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
  \item{estimator}{
The estimator to be used. Currently implemented are \code{"ML"} for maximum likelihood estimation, \code{"FIML"} for full-information maximum likelihood estimation, \code{"ULS"} for unweighted least squares estimation, \code{"WLS"} for weighted least squares estimation, and \code{"DWLS"} for diagonally weighted least squares estimation.
}
  \item{optimizer}{
The optimizer to be used. Usually either \code{"nlminb"} (with box constrains) or \code{"ucminf"} (ignoring box constrains), but any optimizer supported by \code{optimr} can be used.
}
  \item{storedata}{
Logical, should the raw data be stored? Needed for bootstrapping (see \code{\link{bootstrap}}).
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally. 
}
\item{\dots}{
Arguments sent to \code{dlvm1}.
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}

\author{
Sacha Epskamp
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
}
