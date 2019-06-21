\name{varcov}
\alias{varcov}
\alias{cholesky}
\alias{precision}
\alias{prec}
\alias{ggm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance-covariance family of psychonetrics models
}
\description{
This is the family of models that models only a variance-covariance matrix with mean structure. The \code{type} argument can be used to define what model is used: \code{type = "cov"} (default) models a variance-covariance matrix directly, \code{type = "chol"} (alias: \code{cholesky()}) models a Cholesky decomposition, \code{type = "prec"} (alias: \code{precision()}) models a precision matrix, and \code{type = "ggm"} (alias: \code{ggm()}) models a Gaussian graphical model (Epskamp, Rhemtulla and Borsboom, 2017).
}
\usage{
varcov(data, type = c("cov", "chol", "prec", "ggm"), sigma = "full", kappa = "full", 
      omega = "full", lowertri = "full", delta = "full", mu, vars, groups, covs, means, 
      nobs, missing = "listwise", equal = "none", baseline_saturated = TRUE, 
      estimator = "ML", optimizer = "default", storedata = FALSE, WLS.V, sampleStats)
cholesky(\dots)
precision(\dots)
prec(\dots)
ggm(\dots)
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
Only used when \code{type = "cov"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa}{
Only used when \code{type = "prec"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega}{
Only used when \code{type = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri}{
Only used when \code{type = "chol"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
\item{delta}{
Only used when \code{type = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{mu}{
Optional vector encoding the mean structure. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free means, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{vars}{
An optional character vector encoding the variables used in the analyis. Must equal names of the dataset in \code{data}.
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
Logical, should the raw data be stored? Needed for bootstrapping (see \code{bootstrap}).
}
\item{WLS.V}{
Optional WLS weights matrix.
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

\eqn{\mathrm{var}(\boldsymbol{y} ) = \boldsymbol{\Sigma}}{var(y) = sigma}

\eqn{\mathcal{E}( \boldsymbol{y} ) = \boldsymbol{\mu}}{E(y) = mu}

in which the covariance matrix can further be modeled in three ways. With \code{type = "chol"} as Cholesky decomposition:

\eqn{\boldsymbol{\Sigma} = \boldsymbol{L}\boldsymbol{L}}{sigma = lowertri * lowertri},

with \code{type = "prec"} as Precision matrix:

\eqn{\boldsymbol{\Sigma} = \boldsymbol{K}^{-1}}{sigma = kappa^(-1)},

and finally with \code{type = "ggm"} as Gaussian graphical model:

\eqn{\boldsymbol{\Sigma} = \boldsymbol{\Delta}(\boldsymbol{I} - \boldsymbol{\Omega})^(-1) \boldsymbol{\Delta}}{sigma = delta * (I - omega)^(-1) * delta}.
}
\value{
An object of the class psychonetrics 
}
\references{
Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017). Generalized network psychometrics: Combining network and latent variable models. Psychometrika, 82(4), 904-927.
}
\author{
Sacha Epskamp
}


\seealso{
\code{\link{lvm}}, \code{\link{var1}}, \code{\link{dlvm1}}
}
\examples{
# Load bfi data from psych package:
library("psych")
data(bfi)

# Also load dplyr for the pipe operator:
library("dplyr")

# Let's take the agreeableness items, and gender:
ConsData <- bfi \%>\% 
  select(A1:A5, gender) \%>\% 
  na.omit # Let's remove missingness (otherwise use Estimator = "FIML)

# Define variables:
vars <- names(ConsData)[1:5]

# Let's fit an empty GGM:
mod0 <- ggm(ConsData, vars = vars, omega = "empty")

# Run the model:
mod0 <- mod0 \%>\% runmodel

\donttest{
# We can look at the modification indices:
mod0 \%>\% MIs

# To automatically add along modification indices, we can use stepup:
mod1 <- mod0 \%>\% stepup(greedy = TRUE, alpha = 0.005, greedyadjust = "fdr")

# Let's also prune all non-significant edges to finish:
mod1 <- mod1 \%>\% prune(alpha = 0.005)

# Look at the fit:
mod1 \%>\% fit

# Compare to original (baseline) model:
compare(baseline = mod0, adjusted = mod1)

# We can also look at the parameters:
mod1 \%>\% parameters

# Or obtain the network as follows:
getmatrix(mod1, "omega")

# We may also be interested in the stability of our search algorithm.
# We can bootstrap our data and repeat the search as follows:
mod_boot <- ggm(ConsData, vars = vars, omega = "empty", storedata = TRUE) \%>\%
  bootstrap \%>\% # bootstrap data
  runmodel \%>\% # Run model
  stepup(greedy = TRUE, alpha = 0.005, greedyadjust = "fdr") \%>\% # Search algorithm 1
  prune(alpha = 0.005) # Search algorithm 2

# Which may give some different results:
getmatrix(mod_boot, "omega")

# This can be repeated (ideally 100 - 1000 times):
bootstraps <- replicate(10,simplify = FALSE,expr = {
  mod_boot <- ggm(ConsData, vars = vars, omega = "empty", storedata = TRUE) \%>\%
    bootstrap \%>\% # bootstrap data
    runmodel \%>\% # Run model
    stepup(greedy = TRUE, alpha = 0.005, greedyadjust = "fdr") \%>\% # Search algorithm 1
    prune(alpha = 0.005) # Search algorithm 2
  getmatrix(mod_boot, "omega")
})

# Now we can look at, for example, the inclusion probability:
inclusionProportion <- 1/length(bootstraps) * Reduce("+",lapply(bootstraps,function(x)1*(x!=0)))
inclusionProportion
}
}
