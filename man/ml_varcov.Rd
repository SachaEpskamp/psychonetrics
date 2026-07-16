\name{ml_varcov}
\alias{ml_varcov}
\alias{ml_ggm}
\alias{ml_corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-level variance-covariance family of psychonetrics models
}
\description{
This is the multi-level (two-level) analogue of the \code{\link{varcov}} family: it models the within-cluster and between-cluster variance-covariance matrices of clustered data directly, without a latent-variable layer. Each level can be parameterized independently through the \code{within} and \code{between} arguments: \code{"cov"} models a covariance matrix directly, \code{"chol"} a Cholesky decomposition, \code{"prec"} a precision matrix, \code{"ggm"} a Gaussian graphical model (the multi-level GGM; alias \code{ml_ggm()} sets both levels to \code{"ggm"}) and \code{"cor"} a correlation matrix (alias \code{ml_corr()} sets both levels to \code{"cor"}).

The model is estimated with the two-level sufficient-statistics Gaussian maximum-likelihood estimator (the same estimator used by \code{\link{ml_lvm}} with \code{estimator = "ML"}), so the computational cost does not depend on the number of units per cluster.

The multi-level GGM (\code{ml_ggm}) is also implemented, using a different (nodewise multilevel regression) estimator, in the \code{mlGGM} function of the \code{mlVAR} package. The \code{psychonetrics} implementation here estimates the full two-level covariance structure by maximum likelihood instead.

\strong{This framework is experimental} (a message is shown once per session).
}
\usage{
ml_varcov(data, clusters, vars,
          within  = c("cov", "chol", "prec", "ggm", "cor"),
          between = c("cov", "chol", "prec", "ggm", "cor"),
          groups, groupvar,
          sigma_within = "full", lowertri_within = "full",
          omega_within = "full", delta_within = "diag",
          kappa_within = "full", rho_within = "full", SD_within = "full",
          sigma_between = "full", lowertri_between = "full",
          omega_between = "full", delta_between = "diag",
          kappa_between = "full", rho_between = "full", SD_between = "full",
          mu, equal = "none", baseline_saturated = TRUE,
          estimator = "ML", optimizer, storedata = FALSE,
          standardize = c("none", "z", "quantile"), sampleStats,
          verbose = FALSE)

ml_ggm(...)

ml_corr(...)
}
\arguments{
  \item{data}{A data frame in long format: one row per unit, with a column identifying the cluster each unit belongs to.}
  \item{clusters}{String indicating the name of the cluster (grouping) variable in \code{data}. May not be missing.}
  \item{vars}{Character vector of the variables to include. If missing, all columns except \code{clusters} and the grouping variable are used.}
  \item{within}{The parameterization of the within-cluster covariance matrix. One of \code{"cov"}, \code{"chol"}, \code{"prec"}, \code{"ggm"} or \code{"cor"}.}
  \item{between}{The parameterization of the between-cluster covariance matrix. One of \code{"cov"}, \code{"chol"}, \code{"prec"}, \code{"ggm"} or \code{"cor"}.}
  \item{groups}{Deprecated alias for \code{groupvar}.}
  \item{groupvar}{An optional string indicating the name of a variable for a multi-group analysis.}
  \item{sigma_within}{Only used when \code{within = "cov"}. Either \code{"full"}, \code{"diag"} or a matrix specifying the covariance structure of the within block (0 = fixed to 0, 1 = free, higher integers = equality constraints).}
  \item{lowertri_within}{Only used when \code{within = "chol"}. As \code{sigma_within}, for the lower-triangular Cholesky factor.}
  \item{omega_within}{Only used when \code{within = "ggm"}. As \code{sigma_within}, for the within partial-correlation network (only the off-diagonal is used).}
  \item{delta_within}{Only used when \code{within = "ggm"}. As \code{sigma_within}, for the within diagonal scaling matrix.}
  \item{kappa_within}{Only used when \code{within = "prec"}. As \code{sigma_within}, for the within precision matrix.}
  \item{rho_within}{Only used when \code{within = "cor"}. As \code{sigma_within}, for the within correlation matrix.}
  \item{SD_within}{Only used when \code{within = "cor"}. As \code{sigma_within}, for the within standard-deviation scaling matrix.}
  \item{sigma_between}{As \code{sigma_within}, for the between-cluster block (used when \code{between = "cov"}).}
  \item{lowertri_between}{As \code{lowertri_within}, for the between-cluster block.}
  \item{omega_between}{As \code{omega_within}, for the between-cluster block.}
  \item{delta_between}{As \code{delta_within}, for the between-cluster block.}
  \item{kappa_between}{As \code{kappa_within}, for the between-cluster block.}
  \item{rho_between}{As \code{rho_within}, for the between-cluster block.}
  \item{SD_between}{As \code{SD_within}, for the between-cluster block.}
  \item{mu}{Optional vector encoding the mean structure. Set elements to 0 to indicate fixed zeroes.}
  \item{equal}{A character vector indicating which matrices should be constrained equal across groups.}
  \item{baseline_saturated}{A logical indicating whether the baseline and saturated models should be included. Mostly used internally and should be kept to \code{TRUE}.}
  \item{estimator}{Only \code{"ML"} (the two-level sufficient-statistics ML estimator) is supported.}
  \item{optimizer}{The optimizer to be used. Usually left missing to use the default.}
  \item{storedata}{Logical, should the raw data be stored? Needed for bootstrapping.}
  \item{standardize}{Optional grand standardization of the raw variables before estimation: \code{"none"}, \code{"z"} (z-scores) or \code{"quantile"}.}
  \item{sampleStats}{An optional sample statistics object. Mostly used internally.}
  \item{verbose}{Logical, should progress be printed to the console?}
  \item{\dots}{Arguments passed to \code{ml_varcov} by \code{ml_ggm} and \code{ml_corr}.}
}
\details{
For unit \eqn{j} in cluster \eqn{i} the observed response is decomposed into a within-cluster and a between-cluster component,
\deqn{y_{ij} = \mu + \eta^{(W)}_{ij} + \eta^{(B)}_{i},}{y_ij = mu + eta^W_ij + eta^B_i,}
with \eqn{\eta^{(W)}_{ij} \sim N(0, \Sigma_W)}{eta^W_ij ~ N(0, Sigma_W)} and \eqn{\eta^{(B)}_{i} \sim N(0, \Sigma_B)}{eta^B_i ~ N(0, Sigma_B)} independent, so that \eqn{\mathrm{var}(y_{ij}) = \Sigma_W + \Sigma_B}{var(y_ij) = Sigma_W + Sigma_B} and \eqn{\mathrm{cov}(y_{ij}, y_{ij'}) = \Sigma_B}{cov(y_ij, y_ij') = Sigma_B} for two units in the same cluster. \eqn{\Sigma_W}{Sigma_W} (\code{sigma_within}) and \eqn{\Sigma_B}{Sigma_B} (\code{sigma_between}) are each parameterized according to the \code{within} and \code{between} arguments; for the \code{"ggm"} parameterization \eqn{\Sigma = \Delta (I - \Omega)^{-1} \Delta}{Sigma = Delta (I - Omega)^-1 Delta} with partial-correlation network \eqn{\Omega}{Omega} (\code{omega_within} / \code{omega_between}) and scaling \eqn{\Delta}{Delta} (\code{delta_within} / \code{delta_between}).

This is the special case of \code{\link{ml_lvm}} in which the loadings matrix equals the identity and there is no residual variance, but it is implemented as its own framework with its own matrix names (\code{sigma_within}, \code{omega_within}, \code{sigma_between}, \code{omega_between}, ...) and its own (simpler) analytic derivatives.
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class}).
}
\references{
Epskamp, S., Rhemtulla, M., & Borsboom, D. (2017). Generalized network psychometrics: Combining network and latent variable models. \emph{Psychometrika, 82}(4), 904-927. \doi{10.1007/s11336-017-9557-x}

Hamaker, E. L., & Muthen, B. (2020). The fixed versus random effects debate and how it relates to centering in multilevel modeling. \emph{Psychological Methods, 25}(3), 365-379.

McDonald, R. P., & Goldstein, H. (1989). Balanced versus unbalanced designs for linear structural relations in two-level data. \emph{British Journal of Mathematical and Statistical Psychology, 42}(2), 215-232.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
\seealso{
\code{\link{varcov}} for the single-level variance-covariance family, \code{\link{ggm}} for the single-level GGM, and \code{\link{ml_lvm}} for the multi-level latent variable model that supplies the estimator.
}
\examples{
\dontrun{
library("dplyr")

# Simulate two-level data (200 clusters of 8 units, 4 variables):
set.seed(1)
nCluster <- 200; perCluster <- 8; p <- 4
cl <- rep(seq_len(nCluster), each = perCluster)
# True within and between covariance matrices:
SigmaW <- diag(p); SigmaW[1,2] <- SigmaW[2,1] <- 0.4
SigmaB <- diag(0.5, p); SigmaB[3,4] <- SigmaB[4,3] <- 0.2
b <- MASS::mvrnorm(nCluster, rep(0, p), SigmaB)
Y <- b[cl, ] + MASS::mvrnorm(nCluster * perCluster, rep(0, p), SigmaW)
Data <- as.data.frame(Y)
names(Data) <- paste0("V", 1:p)
Data$cluster <- cl

# Fit a multi-level GGM (both levels a Gaussian graphical model):
mod <- ml_ggm(Data, vars = paste0("V", 1:p), clusters = "cluster")
mod <- mod \%>\% runmodel

# Inspect the estimated within- and between-cluster networks:
getmatrix(mod, "omega_within")
getmatrix(mod, "omega_between")

# Prune non-significant edges and re-fit:
mod <- mod \%>\% prune(alpha = 0.05)
}
}
