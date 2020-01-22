\name{var1}
\alias{var1}
\alias{gvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lag-1 vector autoregression family of psychonetrics models
}
\description{
This is the family of models that models time-series data using a lag-1 vector autoregressive model (VAR; Epskamp,Waldorp, Mottus, Borsboom, 2018). The model is fitted to the Toeplitz matrix, but unlike typical SEM software the block of covariances of the lagged variables is not used in estimating the temporal and contemporaneous relationships (the block is modeled completely separately using a cholesky decomposition, and does not enter the model elsewise). The \code{contemporaneous} argument can be used to define what contemporaneous model is used: \code{contemporaneous = "cov"} (default) models a variance-covariance matrix, \code{contemporaneous = "chol"} models a Cholesky decomposition, \code{contemporaneous = "prec"} models a precision matrix, and \code{contemporaneous = "ggm"} (alias: \code{gvar()}) models a Gaussian graphical model, also then known as a graphical VAR model.
}
\usage{
var1(data, contemporaneous = c("cov", "chol", "prec",
                    "ggm"), beta = "full", omega_zeta = "full", delta_zeta
                    = "full", kappa_zeta = "full", sigma_zeta = "full",
                    lowertri_zeta = "full", mu, beepvar, dayvar, idvar,
                    vars, groups, covs, means, nobs, missing = "listwise",
                    equal = "none", baseline_saturated = TRUE, estimator =
                    "ML", optimizer = "default", storedata = FALSE,
                    covtype = c("choose", "ML", "UB"), sampleStats)

gvar(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{contemporaneous}{
The type of contemporaneous model used. See description.
}
  \item{beta}{
A model matrix encoding the temporal relationships (transpose of temporal network). A 0 encodes a fixed to zero element, a 1 encoding a free to estimate element, and higher integers encoding equality constrains. For multiple groups, this argument can be a list or array with each element/slice encoding such a matrix.  Can also be \code{"full"} for a full temporal network or \code{"empty"} for an empty temporal network.
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
  \item{mu}{
Optional vector encoding the mean structure. Set elements to 0 to indicate fixed to zero constrains, 1 to indicate free means, and higher integers to indicate equality constrains. For multiple groups, this argument can be a list or array with each element/column encoding such a vector.
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
  \item{sampleStats}{
An optional sample statistics object. Mostly used internally. 
}
\item{covtype}{
  If 'covs' is used, this is the type of covariance (maximum likelihood or unbiased) the input covariance matrix represents. Set to \code{"ML"} for maximum likelihood estimates (denominator n) and \code{"UB"} to unbiased estimates (denominator n-1). The default will try to find the type used, by investigating which is most likely to result from integer valued datasets.
}
\item{\dots}{
Arguments sent to \code{var1}
}
}
\details{
This will be updated in a later version.
}
\value{
An object of the class psychonetrics
}
\references{
Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, D. (2018). The Gaussian graphical model in cross-sectional and time-series data. Multivariate Behavioral Research, 53(4), 453-480.
}
\author{
Sacha Epskamp
}


\seealso{
\code{\link{lvm}}, \code{\link{varcov}}, \code{\link{dlvm1}}
}


\examples{
library("dplyr")
library("graphicalVAR")

beta <- matrix(c(
  0,0.5,
  0.5,0
),2,2,byrow=TRUE)
kappa <- diag(2)
simData <- graphicalVARsim(50, beta, kappa)

# Form model:
model <- gvar(simData)

# Evaluate model:
model <- model \%>\% runmodel

# Parameter estimates:
model \%>\% parameters

# Note: this example is wrapped in a dontrun environment because the data is not 
# available locally.
\dontrun{
# Longer example:
#
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

# Estimate, prune with FDR, and perform stepup search:
model_FDRprune <- gvar(
  tsdata, 
  vars = vars, 
  dayvar = "Day",
  estimator = "FIML"
  ) \%>\% 
  runmodel \%>\% 
  prune(adjust = "fdr", recursive = FALSE) \%>\% 
  stepup(criterion = "bic")

# Estimate with greedy stepup search:
model_stepup <- gvar(
  tsdata, 
  vars = vars, 
  dayvar = "Day",
  estimator = "FIML",
  omega_zeta = "empty",
  beta = "empty"
) \%>\% 
  runmodel \%>\% 
  stepup(greedy = TRUE, greedyadjust = "bonferroni", criterion = "bic")

# Compare models:
compare(
  FDRprune = model_FDRprune,
  stepup = model_stepup
)
# Very similar but not identical. Stepup is prefered here according to AIC and BIC

# Stepup results:
temporal <- getmatrix(model_stepup, "PDC") # PDC = Partial Directed Correlations
contemporaneous <- getmatrix(model_stepup, "omega_zeta")

# Average layout:
library("qgraph")
L <- averageLayout(temporal, contemporaneous)

# Labels:
labs <- gsub("\\\\.","\\n",vars)

# Plot:
layout(t(1:2))
qgraph(temporal, layout = L, theme = "colorblind", directed=TRUE, diag=TRUE,
       title = "Temporal", vsize = 12, mar = rep(6,4), asize = 5,
       labels = labs)
qgraph(contemporaneous, layout = L, theme = "colorblind", 
       title = "Contemporaneous", vsize = 12, mar = rep(6,4), asize = 5,
       labels = labs)
}
}
