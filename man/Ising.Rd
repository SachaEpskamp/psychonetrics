\name{Ising}
\alias{Ising}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ising model
}
\description{
This is the family of Ising models fit to datasets with two or more ordered integer response options. The classical case is a dichotomous dataset (e.g., encoded with -1 and 1 or with 0 and 1), but the same parameterization---one threshold (\code{tau}) per variable and a pairwise network (\code{omega})---is also used for any number of ordered response options (e.g., \code{c(-1,0,1)} or \code{seq(-5,5)}; the values need not be integers), encoded identically across all variables. Only the partition function changes: it sums over all \code{length(responses)^nNode} response patterns rather than only \code{2^nNode}. Note that the input matters (see also https://arxiv.org/abs/1811.02916) in this model! Models based on a dataset that is encoded with -1 and 1 are not entirely equivalent to models based on datasets encoded with 0 and 1 (non-equivalences occur in multi-group settings with equality constrains).
}
\usage{
Ising(data, omega = "full", tau, beta, beta_model =
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
The network structure. Either \code{"full"} to estimate every element freely, \code{"zero"} to set all elements to zero, or a matrix of the dimensions nNode x nNode with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{tau}{
Optional vector encoding the threshold/intercept structure. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free intercepts, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
}
  \item{beta}{
Optional scalar encoding the inverse temperature. 1 indicate free beta parameters, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such scalers.
}
\item{beta_model}{How should beta be modeled? Set \code{beta_model = "log_beta"} to model the log of beta rather than beta directly.}
  \item{vars}{
An optional character vector encoding the variables used in the analyis. Must equal names of the dataset in \code{data}.
}
  \item{groups}{
An optional character vector encoding the variables used in the analyis. Must equal names of the dataset in \code{data}.
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
A vector of the response options used, encoded identically across all variables (e.g., \code{c(-1,1)}, \code{c(0,1)}, \code{c(-1,0,1)} or \code{seq(-5,5)}). Automatically detected from the data when not supplied; required when \code{covs} is used. Any number of distinct response options is allowed, and the values need not be integers.
}
  \item{missing}{
How should missingness be handled in computing the sample covariances and number of observations when \code{data} is used. Can be \code{"listwise"} for listwise deletion, or \code{"pairwise"} for pairwise deletion. NOT RECOMMENDED TO BE USED YET IN ISING MODEL.
}
  \item{equal}{
A character vector indicating which matrices should be constrained equal across groups. 
}
  \item{baseline_saturated}{
A logical indicating if the baseline and saturated model should be included. Mostly used internally and NOT Recommended to be used manually.
}
   \item{estimator}{
The estimator to be used. Currently implemented are \code{"ML"} for maximum likelihood estimation, \code{"FIML"} for full-information maximum likelihood estimation, \code{"ULS"} for unweighted least squares estimation, \code{"WLS"} for weighted least squares estimation, and \code{"DWLS"} for diagonally weighted least squares estimation. Only ML estimation is currently supported for the Ising model. 
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
The maximum number of response patterns the exact ML estimator may enumerate. Exact ML estimation of the Ising model sums over every possible response pattern when computing the partition function, expected values and (expected) Hessian, so the cost grows as \code{length(responses)^nNode}. The function stops with an error when \code{length(responses)^nNode} exceeds \code{maxStates}. The default \code{2^maxNodes} reproduces the historical binary node limit while accounting for more than two response options. Raise this only if the computation is feasible (note that the number of states---and hence the run time---grows very quickly with both the number of nodes and the number of response options).
}
\item{min_sum}{
 The minimum sum score that is artifically possible in the dataset. Defaults to -Inf. Set this only if you know a lower sum score is not possible in the data, for example due to selection bias.
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
\item{penalize_matrices}{Character vector of matrix names to penalize. If missing, defaults are selected based on the model type.}
}
\details{
The Ising Model takes the following form:

\eqn{\Pr(\boldsymbol{Y} = \boldsymbol{y})  =  \frac{\exp\left( -\beta H\left(\boldsymbol{y}; \boldsymbol{\tau}, \boldsymbol{\Omega}\right)\right)}{Z(\boldsymbol{\tau}, \boldsymbol{\Omega})}}{P(Y=y) = 1/Z exp(-beta H(y))}

With Hamiltonian:

\eqn{H\left(\boldsymbol{y}; \boldsymbol{\tau}, \boldsymbol{\Omega}\right) = -\sum_{i=1}^{m} \tau_i y_{i} - \sum_{i=2}^{m} \sum_{j=1}^{i-1} \omega_{ij} y_i y_j.}{
H(y) = - sum_i tau_i y_i - sum_<i,j> omega_ij y_i y_j
}

And Z representing the partition function or normalizing constant.

The responses \eqn{y_i} need not be dichotomous: they may take any number of ordered values (the same set for every variable, and not necessarily integers), supplied through \code{responses} or detected automatically from the data. The Hamiltonian above is unchanged; only the partition function \eqn{Z} (and the expected values and Hessian derived from it) sums over all \code{length(responses)^nNode} response patterns instead of only \code{2^nNode}. With two response options the model is the usual Ising model and reproduces the same estimates as before. With more than two response options the model remains a maximum-entropy distribution matching the means (via \code{tau}) and the pairwise products (via \code{omega}); it does not model higher-order interactions or the within-variable second moments, and the enumeration cost grows steeply with the number of nodes and response options (see \code{maxStates}).

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


\examples{
\donttest{
library("dplyr")
data("Jonas")

# Variables to use:
vars <- names(Jonas)[1:10]

# Arranged groups to put unfamiliar group first (beta constrained to 1):
Jonas <- Jonas[order(Jonas$group),]

# Form saturated model:
model1 <- Ising(Jonas, vars = vars, groups = "group")

# Run model:
model1 <- model1 \%>\% runmodel
# Note: SEs are approximated automatically because there are zeroes
# in the crosstables of the "Knows Jonas" group, which makes the
# Fisher information matrix singular. A warning is shown when
# this fallback occurs.

# Prune-stepup to find a sparse model:
model1b <- model1 \%>\% prune(alpha = 0.05) \%>\%  stepup(alpha = 0.05)

# Equal networks:
model2 <- model1 \%>\% groupequal("omega") \%>\% runmodel

# Prune-stepup to find a sparse model:
model2b <- model2 \%>\% prune(alpha = 0.05) \%>\% stepup(mi = "mi_equal", alpha = 0.05)

# Equal thresholds:
model3 <- model2 \%>\% groupequal("tau") \%>\% runmodel

# Prune-stepup to find a sparse model:
model3b <- model3 \%>\% prune(alpha = 0.05) \%>\% stepup(mi = "mi_equal", alpha = 0.05)

# Equal beta:
model4 <- model3 \%>\% groupequal("beta") \%>\% runmodel

# Prune-stepup to find a sparse model:
model4b <- model4 \%>\% prune(alpha = 0.05) \%>\% stepup(mi = "mi_equal", alpha = 0.05)

# Compare all models:
compare(
  `1. all parameters free (dense)` = model1,
  `2. all parameters free (sparse)` = model1b,
  `3. equal networks (dense)` = model2,
  `4. equal networks (sparse)` = model2b,
  `5. equal networks and thresholds (dense)` = model3,
  `6. equal networks and thresholds (sparse)` = model3b,
  `7. all parameters equal (dense)` = model4,
  `8. all parameters equal (sparse)` = model4b
) \%>\% arrange(BIC)
}
}
