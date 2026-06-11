\name{BlumeCapel}
\alias{BlumeCapel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Blume-Capel model
}
\description{
The Blume-Capel model is the full \emph{Spin} distribution model. In addition to the linear node potentials (\code{tau}) and pairwise interactions (\code{omega}) of the \code{\link{Ising}} model, it estimates a quadratic node potential \code{delta} (the Blume-Capel single-ion / crystal-field term). The Ising model is the special case obtained by fixing all \code{delta} parameters to zero, so the two models are nested. By convention the Blume-Capel model is defined for the ternary \code{c(-1, 0, 1)} response coding (a warning is raised for any other response set), but any number of ordered numeric response options is allowed (the same set for every variable, and not necessarily integers).
}
\usage{
BlumeCapel(data, omega = "full", tau, delta, beta, beta_model =
                 c("beta", "log_beta"), vars, groups, covs, means,
                 nobs, covtype = c("choose", "ML", "UB"), responses,
                 missing = "listwise", equal = "none",
                 baseline_saturated = TRUE, estimator = "default",
                 optimizer, storedata = FALSE, WLS.W, sampleStats,
                 identify = TRUE, verbose = FALSE, maxNodes = 20,
                 maxStates = 2^maxNodes, min_sum = -Inf,
                 bootstrap = FALSE, boot_sub,
                 boot_resample, penalty_lambda = NA,
                penalty_alpha = 1, penalize_matrices)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{omega}{
The network structure. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions nNode x nNode with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constraints. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{tau}{
Optional vector encoding the threshold/intercept (linear node potential) structure. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{delta}{
Optional vector encoding the quadratic node potential (Blume-Capel) structure. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free parameters, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such a vector. When missing, every \code{delta} parameter is freely estimated (fixing all of them to zero would instead yield the \code{\link{Ising}} model).
}
  \item{beta}{
Optional scalar encoding the inverse temperature. 1 indicate free beta parameters, and higher integers to indicate equality constraints. For multiple groups, this argument can be a list or array with each element/column encoding such scalars.
}
\item{beta_model}{How should beta be modeled? Set \code{beta_model = "log_beta"} to model the log of beta rather than beta directly.}
  \item{vars}{
An optional character vector encoding the variables used in the analysis. Must equal names of the dataset in \code{data}.
}
  \item{groups}{
An optional string indicating the name of the group variable in the data.
}
  \item{covs}{
A sample variance--covariance matrix, or a list/array of such matrices for multiple groups. Make sure \code{covtype} argument is set correctly to the type of covariances used.
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
  \item{responses}{
A vector of the response options used, encoded identically across all variables. Defaults to (and is conventionally) \code{c(-1, 0, 1)}; a warning is raised for any other set. Automatically detected from the data when not supplied; required when \code{covs} is used. Any number of distinct response options is allowed, and the values need not be integers.
}
  \item{missing}{
How should missingness be handled in computing the sample covariances and number of observations when \code{data} is used. Can be \code{"listwise"} for listwise deletion, or \code{"pairwise"} for pairwise deletion. NOT RECOMMENDED TO BE USED YET IN THIS MODEL.
}
  \item{equal}{
A character vector indicating which matrices should be constrained equal across groups.
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
   \item{estimator}{
The estimator to be used. Only \code{"ML"} (maximum likelihood) estimation is currently supported for the Blume-Capel model.
}
  \item{optimizer}{
The optimizer to be used. Can be one of \code{"nlminb"} (the default R \code{nlminb} function), \code{"ucminf"} (from the \code{optimr} package), \code{"nloptr_TNEWTON"} (preconditioned truncated Newton via \code{nloptr}), and \code{"LBFGS++"} (pure C++ L-BFGS-B). Defaults to \code{"nlminb"}.
}
  \item{storedata}{
Logical, should the raw data be stored? Needed for bootstrapping (see \code{bootstrap}).
}
\item{WLS.W}{
Optional WLS weights matrix. CURRENTLY NOT USED.
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally.
}
  \item{identify}{
Logical, should the model be identified?
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{maxNodes}{
The maximum number of nodes allowed in the analysis. Used to set the default of \code{maxStates} (\code{2^maxNodes}); for binary data this reproduces the historical node limit. It is not recommended to set this higher.
}
\item{maxStates}{
The maximum number of response patterns the exact ML estimator may enumerate. Exact ML estimation sums over every possible response pattern when computing the partition function, expected values and (expected) Hessian, so the cost grows as \code{length(responses)^nNode}. The function stops with an error when \code{length(responses)^nNode} exceeds \code{maxStates}. The default \code{2^maxNodes} reproduces the historical binary node limit while accounting for more than two response options. Raise this only if the computation is feasible (note that the number of states---and hence the run time---grows very quickly with both the number of nodes and the number of response options).
}
\item{min_sum}{
 The minimum sum score that is artificially possible in the dataset. Defaults to -Inf. Set this only if you know a lower sum score is not possible in the data, for example due to selection bias.
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
\item{penalty_lambda}{Numeric penalty strength for penalized ML estimation (PML/PFIML). \code{NA} (default) triggers automatic selection via EBIC-based grid search when a penalized estimator is used; set to a specific numeric value to use a fixed penalty strength (0 = no penalty). See \code{\link{find_penalized_lambda}} and \code{\link{penalize}}.}
\item{penalty_alpha}{Elastic net mixing parameter: 1 = LASSO (default), 0 = ridge.}
\item{penalize_matrices}{Character vector of matrix names to penalize. If missing, defaults are selected based on the model type (the network \code{omega}).}
}
\details{
The Blume-Capel model is built on the \emph{Spin} distribution:

\eqn{\Pr(\boldsymbol{Y} = \boldsymbol{y})  =  \frac{\exp\left( -\beta H\left(\boldsymbol{y}; \boldsymbol{\tau}, \boldsymbol{\delta}, \boldsymbol{\Omega}\right)\right)}{Z(\boldsymbol{\tau}, \boldsymbol{\delta}, \boldsymbol{\Omega})}}{P(Y=y) = 1/Z exp(-beta H(y))}

With Hamiltonian:

\eqn{H\left(\boldsymbol{y}; \boldsymbol{\tau}, \boldsymbol{\delta}, \boldsymbol{\Omega}\right) = -\sum_{i=1}^{m} \tau_i y_{i} + \sum_{i=1}^{m} \delta_i y_{i}^2 - \sum_{i=2}^{m} \sum_{j=1}^{i-1} \omega_{ij} y_i y_j.}{
H(y) = - sum_i tau_i y_i + sum_i delta_i y_i^2 - sum_<i,j> omega_ij y_i y_j
}

and Z representing the partition function or normalizing constant. Equivalently,

\eqn{\Pr(\boldsymbol{Y} = \boldsymbol{y}) \propto \exp\left( \sum_i \tau_i y_i - \sum_i \delta_i y_i^2 + \sum_{i<j} \omega_{ij} y_i y_j \right).}{P(Y=y) prop. exp(sum tau_i y_i - sum delta_i y_i^2 + sum omega_ij y_i y_j).}

For the ternary states \code{c(-1, 0, 1)} this is the usual Blume-Capel parameterization: \code{tau} controls the tendency toward the positive versus the negative state, \code{delta} controls the tendency toward the zero/middle state versus the active/extreme states (positive \code{delta} favours the middle category), and \code{omega} controls pairwise alignment. For more general ordered numeric states the same parameterization defines an ordinal spin model.

The \code{\link{Ising}} model is nested within the Blume-Capel model by fixing all \code{delta} parameters to zero. Exact ML estimation enumerates every response pattern, so the cost grows as \code{length(responses)^nNode} (see \code{maxStates}); the model is therefore limited to relatively small networks.
}
\value{
An object of the class psychonetrics
}
\references{
Epskamp, S., Maris, G., Waldorp, L. J., & Borsboom, D. (2018). Network Psychometrics. In: Irwing, P., Hughes, D., & Booth, T. (Eds.), The Wiley Handbook of Psychometric Testing, 2 Volume Set: A Multidisciplinary Reference on Survey, Scale and Test Development. New York: Wiley.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{Ising}}
}
\examples{
\donttest{
library("dplyr")

# Simulate a small Blume-Capel data set (requires the development version of
# IsingSampler that supports the 'delta' argument):
\dontrun{
library("IsingSampler")
nNode  <- 5
graph  <- 0.4 * (matrix(1, nNode, nNode) - diag(nNode)) *
            (matrix(c(0,1,0,0,1, 1,0,1,0,0, 0,1,0,1,0, 0,0,1,0,1, 1,0,0,1,0), 5, 5))
data <- IsingSampler(1000, graph = graph, thresholds = rep(0, nNode),
                     beta = 1, delta = 1, responses = c(-1, 0, 1))

# Fit the Blume-Capel model:
mod <- BlumeCapel(data, responses = c(-1, 0, 1)) \%>\% runmodel

# Estimated network:
getmatrix(mod, "omega")

# Estimated quadratic node potentials:
getmatrix(mod, "delta")
}
}
}
