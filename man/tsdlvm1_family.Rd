\name{tsdlvm1}
\alias{tsdlvm1}
\alias{ts_lvgvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lag-1 dynamic latent variable model family of psychonetrics models for time-series data
}
\description{
This is the family of models that models a dynamic factor model on time-series. There are two covariance structures that can be modeled in different ways: \code{contemporaneous} for the contemporaneous model and \code{residual} for the residual model. These can be set to \code{"cov"} for covariances, \code{"prec"} for a precision matrix, \code{"ggm"} for a Gaussian graphical model and \code{"chol"} for a Cholesky decomposition. The \code{ts_lvgvar} wrapper function sets \code{contemporaneous = "ggm"} for the graphical VAR model.
}
\usage{
tsdlvm1(data, lambda, contemporaneous = c("cov", "chol",
                    "prec", "ggm"), residual = c("cov", "chol", "prec",
                    "ggm"), beta = "full", omega_zeta = "full", delta_zeta
                    = "full", kappa_zeta = "full", sigma_zeta = "full",
                    lowertri_zeta = "full", omega_epsilon = "empty",
                    delta_epsilon = "empty", kappa_epsilon = "empty",
                    sigma_epsilon = "empty", lowertri_epsilon = "empty",
                    nu, mu_eta, identify = TRUE, identification =
                    c("loadings", "variance"), latents, beepvar, dayvar,
                    idvar, vars, groups, covs, means, nobs, missing =
                    "listwise", equal = "none", baseline_saturated = TRUE,
                    estimator = "ML", optimizer = "default", storedata =
                    FALSE, sampleStats, covtype = c("choose", "ML", "UB"),
                    centerWithin = FALSE, standardize = c("none", "z",
                    "quantile"))

ts_lvgvar(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{lambda}{
A model matrix encoding the factor loading structure. Each row indicates an indicator and each column a latent. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{contemporaneous}{
The type of contemporaneous model used. See description.
}
  \item{residual}{
The type of residual model used. See description.
}
  \item{beta}{
A model matrix encoding the temporal relationships (transpose of temporal network) between latent variables. A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.  Can also be \code{"full"} for a full temporal network or \code{"empty"} for an empty temporal network.
}
  \item{omega_zeta}{
Only used when \code{contemporaneous = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_zeta}{
Only used when \code{contemporaneous = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_zeta}{
Only used when \code{contemporaneous = "prec"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_zeta}{
Only used when \code{contemporaneous = "cov"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_zeta}{
Only used when \code{contemporaneous = "chol"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{omega_epsilon}{
Only used when \code{residual = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{delta_epsilon}{
Only used when \code{residual = "ggm"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to set all elements to zero, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{kappa_epsilon}{
Only used when \code{residual = "prec"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{sigma_epsilon}{
Only used when \code{residual = "cov"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{lowertri_epsilon}{
Only used when \code{residual = "chol"}. Either \code{"full"} to estimate every element freely, \code{"empty"} to only include diagonal elements, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element, 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.
}
  \item{nu}{
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
  \item{beepvar}{
Optional string indicating assessment beep per day. Adding this argument will cause non-consecutive beeps to be treated as missing!
}
  \item{dayvar}{
Optional string indicating assessment day. Adding this argument makes sure that the first measurement of a day is not regressed on the last measurement of the previous day. IMPORTANT: only add this if the data has multiple observations per day.
}
  \item{idvar}{
Optional string indicating the subject ID
}
  \item{vars}{
An optional character vector encoding the variables used in the analyis. Must equal names of the dataset in \code{data}.
}
  \item{groups}{
An optional string indicating the name of the group variable in \code{data}.
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
  \item{standardize}{
Which standardization method should be used? \code{"none"} (default) for no standardization, \code{"z"} for z-scores, and \code{"quantile"} for a non-parametric transformation to the quantiles of the marginal standard normal distribution.
}
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally. 
}
  \item{centerWithin}{
Logical, should data be within-person centered?
}
\item{covtype}{
  If 'covs' is used, this is the type of covariance (maximum likelihood or unbiased) the input covariance matrix represents. Set to \code{"ML"} for maximum likelihood estimates (denominator n) and \code{"UB"} to unbiased estimates (denominator n-1). The default will try to find the type used, by investigating which is most likely to result from integer valued datasets.
}
\item{\dots}{
Arguments sent to \code{tsdlvm1}
}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}

\author{
Sacha Epskamp
}
\examples{

# Note: this example is wrapped in a dontrun environment because the data is not 
# available locally.
\dontrun{
# Obtain the data from:
#
# Epskamp, S., van Borkulo, C. D., van der Veen, D. C., Servaas, M. N., Isvoranu, A. M., 
# Riese, H., & Cramer, A. O. (2018). Personalized network modeling in psychopathology: 
# The importance of contemporaneous and temporal connections. Clinical Psychological 
# Science, 6(3), 416-427.
# 
# Available here: https://osf.io/c8wjz/
tsdata <- read.csv("Supplementary2_data.csv")

# Encode time variable in a way R understands:
tsdata$time <- as.POSIXct(tsdata$time, tz = "Europe/Amsterdam")

# Extract days:
tsdata$Day <- as.Date(tsdata$time, tz = "Europe/Amsterdam")

# Variables to use:
vars <- c("relaxed", "sad", "nervous", "concentration", "tired", "rumination", 
          "bodily.discomfort")

# Create lambda matrix (in this case: one factor):
Lambda <- matrix(1,7,1)

# Estimate dynamical factor model:
model <- tsdlvm1(
  tsdata, 
  lambda = Lambda,
  vars = vars, 
  dayvar = "Day",
  estimator = "FIML"
)

# Run model:
model <- model \%>\% runmodel

# Look at fit:
model \%>\% print
model \%>\% fit # Pretty bad fit
}
}