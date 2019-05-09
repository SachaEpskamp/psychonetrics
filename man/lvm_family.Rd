\name{lvm}
\alias{lvm}
\alias{lnm}
\alias{rnm}
\alias{lrnm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Continuous latent variable family of psychonetrics models
}
\description{
This is the family of models that models the data as a structural equation model (SEM), allowing the latent and residual variance-covariance matrices to be further modeled as networks. The \code{latent} and \code{residual} arguments can be used to define what latent and residual models are used respectively: \code{"cov"} (default) models a variance-covariance matrix directly, \code{"chol"} models a Cholesky decomposition, \code{"prec"} models a precision matrix, and \code{"ggm"} models a Gaussian graphical model (Epskamp, Rhemtulla and Borsboom, 2017). The wrapper \code{lnm()} sets \code{latent = "ggm"} for the latent network model (LNM), the wrapper  \code{rnm()} sets \code{residual = "ggm"} for the residual network model (RNM), and the wrapper  \code{lrnm()} combines the LNM and RNM.
}
\usage{
lvm(data, lambda, latent = c("cov", "chol", "prec", "ggm"), residual = c("cov", "chol", "prec", "ggm"), sigma_zeta = "full", kappa_zeta = "full", omega_zeta = "full", lowertri_zeta = "full", delta_zeta = "full", sigma_epsilon = "empty", kappa_epsilon = "empty", omega_epsilon = "empty", lowertri_epsilon = "empty", delta_epsilon = "empty", beta = "empty", tau, tau_eta, identify = TRUE, identification = c("loadings", "variance"), vars, latents, groups, covs, means, nobs, missing = "listwise", equal = "none", baseline_saturated = TRUE, estimator = "ML", optimizer = "default", storedata = FALSE, WLS.V, sampleStats)

lnm(\dots)
rnm(\dots)
lrnm(\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{lambda}{
A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{latent}{
The type of latent model used. See description.
}
  \item{residual}{
The type of residual model used. See description.
}
  \item{sigma_zeta}{
Only used when \code{latent = "cov"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta}{
Only used when \code{latent = "prec"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_zeta}{
Only used when \code{latent = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta}{
Only used when \code{latent = "chol"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta}{
Only used when \code{latent = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon}{
Only used when \code{residual = "cov"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon}{
Only used when \code{residual = "prec"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_epsilon}{
Only used when \code{residual = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon}{
Only used when \code{residual = "chol"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon}{
Only used when \code{residual = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{beta}{
A model matrix encoding the structural relations between latent variables. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{tau}{
Optional vector encoding the intercepts of the observed variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{tau_eta}{
Optional vector encoding the intercepts of the latent variables. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{identify}{
Logical, should the model be automatically identified?
}
  \item{identification}{
Type of identification used. \code{"loadings"} to fix the first factor loadings to 1, and \code{"variance"} to fix the diagonal of the latent variable model matrix (sigma_zeta, lowertri_zeta, delta_zeta or kappa_zeta) to 1.
}
  \item{vars}{
An optional character vector encoding the variables used in the analyis. Must equal names of the dataset in \code{data}.
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
  \item{WLS.V}{
The weights matrix used in WLS estimation (experimental)
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally. 
}
\item{\dots}{
Arguments sent to \code{varcov}
}
}
\details{
The model used in this family is:

\eqn{\mathrm{var}( \pmb{y} ) = \pmb{\Lambda} (\pmb{I} - \pmb{B})^{-1} \pmb{\Sigma}_{\zeta}  (\pmb{I} - \pmb{B})^{-1\top}  \pmb{\Lambda}^{\top} +  \pmb{\Sigma}_{\varepsilon} }{var(y) = lambda * (I - beta)^(-1) * sigma_zeta * (I - beta)^(-1') * lambda' + sigma_epsilon}

\eqn{\mathcal{E}( \pmb{y} ) = \pmb{\tau} +  \pmb{\Lambda}  (\pmb{I} - \pmb{B})^{-1} \pmb{\tau}_eta}{E(y) = tau + lambda * (I - beta)^(-1) * tau_eta}

in which the latent covariance matrix can further be modeled in three ways. With \code{latent = "chol"} as Cholesky decomposition:

\eqn{\pmb{\Sigma}_{\zeta} = \pmb{L}_{\zeta}\pmb{L}_{\zeta}}{sigma_zeta = lowertri_zeta * lowertri_zeta},

with \code{latent = "prec"} as Precision matrix:

\eqn{\pmb{\Sigma}_{\zeta} = \pmb{K}_{\zeta}^{-1}}{sigma_zeta = kappa_zeta^(-1)},

and finally with \code{latent = "ggm"} as Gaussian graphical model:

\eqn{\pmb{\Sigma}_{\zeta} = \pmb{\Delta}_{\zeta}(\pmb{I} - \pmb{\Omega}_{\zeta})^(-1) \pmb{\Delta}_{\zeta}}{sigma_zeta = delta_zeta * (I - omega_zeta)^(-1) * delta_zeta}.

Likewise, the residual covariance matrix can also further be modeled in three ways. With \code{residual = "chol"} as Cholesky decomposition:

\eqn{\pmb{\Sigma}_{\varepsilon} = \pmb{L}_{\varepsilon}\pmb{L}_{\varepsilon}}{sigma_epsilon = lowertri_epsilon * lowertri_epsilon},

with \code{latent = "prec"} as Precision matrix:

\eqn{\pmb{\Sigma}_{\varepsilon} = \pmb{K}_{\varepsilon}^{-1}}{sigma_epsilon = kappa_epsilon^(-1)},

and finally with \code{latent = "ggm"} as Gaussian graphical model:

\eqn{\pmb{\Sigma}_{\varepsilon} = \pmb{\Delta}_{\varepsilon}(\pmb{I} - \pmb{\Omega}_{\varepsilon})^(-1) \pmb{\Delta}_{\varepsilon}}{sigma_epsilon = delta_epsilon * (I - omega_epsilon)^(-1) * delta_epsilon}.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\references{
Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017). Generalized network psychometrics: Combining network and latent variable models. Psychometrika, 82(4), 904-927.
}
\author{
Sacha Epskamp
}

\examples{

}
