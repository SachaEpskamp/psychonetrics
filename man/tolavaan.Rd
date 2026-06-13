\name{tolavaan}
\alias{tolavaan}
\title{
Convert a psychonetrics lvm to a lavaan model
}
\description{
Converts a psychonetrics latent variable model (\code{\link{lvm}}) with \code{latent = "cov"} and \code{residual = "cov"} into a fitted \code{lavaan} object, or into a lavaan parameter table. Only the public lavaan API (\code{lavaan::lavaan}) is used; no lavaan source code is reused.
}
\usage{
tolavaan(x, fit = TRUE, verbose = FALSE, \dots)
}
\arguments{
  \item{x}{
A \code{psychonetrics} object created with \code{\link{lvm}} (or a wrapper such as \code{lnm}) using \code{latent = "cov"} and \code{residual = "cov"}.
}
  \item{fit}{
Logical. If \code{TRUE} (default), the model is fitted with \code{lavaan::lavaan} and the resulting \code{lavaan} object is returned. If \code{FALSE}, the lavaan parameter table (a \code{data.frame}) is returned without fitting.
}
  \item{verbose}{
Logical. If \code{TRUE}, an experimental-feature note is shown.
}
  \item{\dots}{
Further arguments passed to \code{lavaan::lavaan}. For example, pass \code{do.fit = FALSE} to obtain an unfitted lavaan object that uses the psychonetrics estimates as starting values.
}
}
\details{
Only models with \code{latent = "cov"} and \code{residual = "cov"} can be converted, because network (\code{ggm}/\code{prec}) and Cholesky (\code{chol}) parameterizations have no direct lavaan equivalent. The function emits one lavaan parameter-table row for every element of the psychonetrics matrices \code{lambda}, \code{sigma_zeta}, \code{sigma_epsilon}, \code{beta}, \code{nu} and \code{nu_eta} (including elements that are fixed to zero, so that lavaan adds no default parameters). Free parameters are given sequential \code{free} indices and fixed parameters their fixed value via \code{ustart}. Equality constraints (shared psychonetrics \code{par} indices) are encoded with per-row \code{plabel}s and explicit \code{==} rows, which is the only encoding lavaan honors for cross-row equality.

\bold{Warm start.} When the psychonetrics model has been estimated (\code{x@computed} is \code{TRUE}), the estimates are also written as starting values (\code{ustart}) for the free parameters. lavaan then typically polishes this warm start by a few optimizer steps, so the refitted log-likelihood matches to within optimizer tolerance (about \eqn{10^{-5}}) rather than exactly.

\bold{Data.} If the psychonetrics model stored its raw data (fitted with \code{storedata = TRUE}), that raw data is passed to lavaan (with the grouping variable for multi-group models, and \code{missing = "ml"} for FIML). Otherwise the sample covariance matrix (and the sample means, when a mean structure is present) is passed with \code{sample.cov.rescale = FALSE} and \code{likelihood = "normal"}, since psychonetrics covariances already use the \eqn{n} (maximum-likelihood) normalization. FIML conversion therefore requires the raw data: re-fit the lvm with \code{storedata = TRUE} first.

\bold{Unsupported.} The following raise an informative error: non-\code{lvm} models, \code{lvm} models with \code{latent} or \code{residual} not equal to \code{"cov"}, correlation-input or ordinal models, estimators other than \code{"ML"} or \code{"FIML"}, and FIML models fitted without \code{storedata = TRUE}.
}
\value{
A fitted \code{lavaan} object when \code{fit = TRUE}, or a lavaan parameter table (\code{data.frame}) when \code{fit = FALSE}.
}
\references{
Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. \doi{10.18637/jss.v048.i02}
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{fromlavaan}}, \code{\link{lvm}}
}
\examples{
\donttest{
if (requireNamespace("lavaan", quietly = TRUE)){
  library("lavaan")
  library("dplyr")

  # A psychonetrics CFA:
  data(HolzingerSwineford1939, package = "lavaan")
  HS <- lavaan::HolzingerSwineford1939
  Lambda <- matrix(0, 9, 3)
  Lambda[1:3, 1] <- Lambda[4:6, 2] <- Lambda[7:9, 3] <- 1

  mod <- lvm(HS, lambda = Lambda, latent = "cov", residual = "cov",
             vars = paste0("x", 1:9),
             latents = c("visual", "textual", "speed"),
             identification = "loadings") \%>\%
    runmodel

  # Convert to a fitted lavaan object:
  lavfit <- tolavaan(mod)
  c(psychonetrics = mod@fitmeasures$logl,
    lavaan = as.numeric(logLik(lavfit)))

  # Or obtain just the parameter table:
  PT <- tolavaan(mod, fit = FALSE)
  head(PT)
}
}
}
