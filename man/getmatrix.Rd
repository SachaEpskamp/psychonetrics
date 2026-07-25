\name{getmatrix}
\alias{getmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract an estimated matrix
}
\description{
This function will extract an estimated matrix, and will either return a single matrix for single group models or a list of such matrices for multiple group models. 
}
\usage{
getmatrix(x, matrix, group, threshold = FALSE, alpha = 0.01,
           adjust = c("none", "holm", "hochberg", "hommel",
           "bonferroni", "BH", "BY", "fdr"), mode = c("tested",
           "all"), diag = TRUE, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{matrix}{
String indicating the matrix to be extracted. For \code{\link{ri_clpm}} / \code{ri_clpn} models this may also be a block name such as \code{"omega_zeta_within"} or \code{"omega_zeta_between"}; see Details.
}
  \item{group}{
Integer indicating the group for the matrix to be extracted.
}
\item{threshold}{Logical. Should the matrix be thresholded (non-significant values set to zero? Can also be a value with an absolute threshold below which parameters are set to zero.)}
  \item{alpha}{
Significance level to use.
}
  \item{adjust}{
p-value adjustment method to use. See \code{p.adjust}.
}
\item{mode}{
Mode for adjusting for multiple comparisons. Should all parameters be considered as the total number of tests or only the tested parameters (parameters of interest)?
}
\item{diag}{Set to FALSE to set diagonal elements to zero.}
\item{time}{
Only for \code{\link{ri_clpm}} / \code{ri_clpn} models: the wave whose block should be returned, either an integer or a column name of the design matrix. See Details.
}
}
\details{
\bold{RI-CLPM / RI-CLPN blocks.} A model made with \code{\link{ri_clpm}} or \code{ri_clpn} is built on the \code{\link{lvm}} framework, so its model matrices span all waves at once: the \eqn{m} within-person components of each of the \eqn{T} waves, followed by the \eqn{m} random intercepts. The quantities of interest are the \eqn{m \times m} blocks of those matrices, which can be addressed directly with the same matrix names that \code{\link{panelvar}} and \code{\link{panelgvar}} use, together with the \code{time} argument:

\itemize{
\item \code{getmatrix(x, "beta", time = t)} -- the temporal (cross-lagged) network of the transition from wave \eqn{t} to wave \eqn{t+1}, with rows the outcome and columns the predictor, exactly as in \code{panelvar}.
\item \code{getmatrix(x, "PDC", time = t)} -- the partial directed correlations of that transition (rows = from, columns = to), the form that feeds directly into \code{qgraph}. The RI-CLPM has no PDC parameters of its own, so these are derived from the temporal block and the innovation precision of the target wave \eqn{t+1}; for a fully stationary model this reproduces the \code{panelvar} PDC exactly.
\item \code{getmatrix(x, "omega_zeta_within", time = t)} -- the contemporaneous network of wave \eqn{t}. Note that at \code{time = 1} this is the \emph{exogenous} wave-1 block, which the RI-CLPM leaves free, and not an innovation structure. The other parameterizations (\code{sigma_}, \code{kappa_}, \code{lowertri_}, \code{rho_}, \code{SD_}, \code{delta_zeta_within}) work the same way, and the plain name plus \code{time} (e.g. \code{"omega_zeta"}) is accepted as a shorthand.
\item \code{getmatrix(x, "omega_zeta_between")} -- the between-person (random intercept) network. No \code{time} is given or accepted, since the random intercepts do not vary over waves.
\item \code{getmatrix(x, "nu", time = t)} (or \code{"nu_within"}) -- the wave-specific intercepts, and \code{getmatrix(x, "nu_between")} -- the random-intercept means that \code{ri_clpm_stationary(x, "intercepts")} frees.
}

Blocks are returned with the variable names of the design matrix as row and column names. Omitting \code{time} returns a named list with one block per wave; waves whose parameters were equated by \code{\link{ri_clpm_stationary}} are collapsed into one entry, so a stationary structure yields a single block rather than several identical copies. Supplying a plain matrix name without \code{time} (for example \code{getmatrix(x, "beta")}) returns the full latent-sized matrix as in earlier versions.

Thresholding is applied to the full matrix and the block is then cut out of the result, so a thresholded block is exactly the corresponding submatrix of the thresholded full matrix and the multiple-comparison adjustment is unaffected by which block is requested. \code{diag = FALSE}, in contrast, applies to the block, so that it removes the autoregressive effects of a temporal block (which do not lie on the diagonal of the full \code{beta} matrix).
}
\value{
A matrix of parameter estimates, or a list of such matrices for multiple group models.
}
\author{
Sacha Epskamp
}
\examples{# Load bfi data from psych package:
library("psychTools")
data(bfi)

# Also load dplyr for the pipe operator:
library("dplyr")

# Let's take the agreeableness items, and gender:
ConsData <- bfi \%>\% 
  select(A1:A5, gender) \%>\% 
  na.omit # Let's remove missingness (otherwise use Estimator = "FIML)

# Define variables:
vars <- names(ConsData)[1:5]

# Let's fit a full GGM:
mod <- ggm(ConsData, vars = vars, omega = "full")

# Run model:
mod <- mod \%>\% runmodel

# Obtain network:
mod \%>\% getmatrix("omega")
}