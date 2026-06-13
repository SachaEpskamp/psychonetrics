\name{fromlavaan}
\alias{fromlavaan}
\title{
Convert a lavaan model to a psychonetrics lvm
}
\description{
Converts a \code{lavaan} model into an equivalent psychonetrics latent variable model (\code{\link{lvm}}) with \code{latent = "cov"} and \code{residual = "cov"}. The input can be a fitted \code{lavaan} object (e.g., from \code{lavaan::cfa} or \code{lavaan::sem}), or lavaan model syntax together with a data frame. Only the public lavaan API is used; no lavaan source code is reused.
}
\usage{
fromlavaan(x, data = NULL, run = FALSE, estimator = "default",
           baseline_saturated = TRUE, verbose = FALSE, \dots)
}
\arguments{
  \item{x}{
A fitted \code{lavaan} object, or a character string with lavaan model syntax. When syntax is supplied, the model is set up internally with \code{lavaan::sem(model = x, data = data, meanstructure = TRUE, fixed.x = FALSE, do.fit = FALSE)}.
}
  \item{data}{
A data frame. Required when \code{x} is model syntax; ignored when \code{x} is a fitted lavaan object.
}
  \item{run}{
Logical. If \code{TRUE}, the resulting psychonetrics model is estimated with \code{\link{runmodel}} before being returned.
}
  \item{estimator}{
The estimator passed to \code{\link{lvm}}. The default (\code{"default"}) uses \code{"ML"}, or \code{"FIML"} when lavaan was fitted with full-information maximum likelihood for missing data.
}
  \item{baseline_saturated}{
Logical, passed to \code{\link{lvm}}; should the baseline and saturated models be included? Mostly used internally.
}
  \item{verbose}{
Logical. If \code{TRUE}, an experimental-feature note is shown and progress is printed.
}
  \item{\dots}{
Further arguments passed to \code{\link{lvm}}.
}
}
\details{
The converted model reproduces the lavaan specification using the LISREL-style model matrices that lavaan exposes through \code{lavInspect(fit, "free")} and \code{lavInspect(fit, "est")}: factor loadings (\code{lambda}), the latent variance--covariance matrix (\code{sigma_zeta}, lavaan's psi), the residual variance--covariance matrix (\code{sigma_epsilon}, lavaan's theta), structural regressions (\code{beta}), observed intercepts (\code{nu}) and latent intercepts (\code{nu_eta}, lavaan's alpha). Observed variables that appear as structural (regression) variables are carried as phantom latent columns; any latent whose name clashes with an observed variable is renamed to \code{eta_<name>}. Equality constraints are detected via shared parameter labels and simple \code{==} rows and mapped onto shared psychonetrics parameter (\code{par}) indices.

The following are guaranteed to match lavaan (for a standard ML fit, see the caveats below): the parameter estimates, the (expected-information) standard errors, the chi-square test statistic and degrees of freedom, the number of free parameters, and the log-likelihood.

\bold{Caveats.} (1) If the lavaan model has no meanstructure, psychonetrics adds a saturated mean structure (free \code{nu}, \code{nu_eta} fixed to 0): estimates, standard errors, chi-square and df still match, but the log-likelihood differs by a constant; a warning is issued. (2) With \code{fixed.x = TRUE} and exogenous covariates, psychonetrics treats the covariates as random (endogenous): estimates, standard errors and chi-square match, but df-based fit measures and the log-likelihood do not (lavaan conditions on the covariates); a warning is issued. (3) If lavaan used \code{likelihood = "wishart"} (the \eqn{n-1} normalization), psychonetrics uses the normal (\eqn{n}) likelihood, so the log-likelihood differs; a warning is issued. (4) Robust standard errors and robust/scaled test statistics are not transferred; the model is converted as plain ML, with a warning. (5) For FIML, psychonetrics uses expected information for standard errors, whereas lavaan defaults to observed information, so FIML standard errors may differ slightly even though the log-likelihood matches.

\bold{Unsupported.} The following lavaan features cannot be converted and raise an informative error: two-level / multilevel models (\code{cluster=}, \code{level:}), defined parameters (\code{:=}), inequality constraints (\code{<}, \code{>}), general (non-simple) equality constraints, formative indicators (\code{<~}), \code{conditional.x = TRUE}, and categorical (ordinal) endogenous variables.
}
\value{
An object of class \code{psychonetrics} (see \code{\link{psychonetrics-class}}).
}
\references{
Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. \doi{10.18637/jss.v048.i02}
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{tolavaan}}, \code{\link{lvm}}
}
\examples{
\donttest{
if (requireNamespace("lavaan", quietly = TRUE)){
  library("lavaan")

  # A confirmatory factor analysis fitted in lavaan:
  data(HolzingerSwineford1939)
  HS <- HolzingerSwineford1939
  fit <- cfa(' visual  =~ x1 + x2 + x3
               textual =~ x4 + x5 + x6
               speed   =~ x7 + x8 + x9 ',
             data = HS, meanstructure = TRUE)

  # Convert to a psychonetrics lvm and run it:
  mod <- fromlavaan(fit, run = TRUE)
  fit(mod)

  # The log-likelihood matches lavaan:
  c(lavaan = as.numeric(logLik(fit)), psychonetrics = mod@fitmeasures$logl)

  # Equivalent starting from syntax + data (without fitting in lavaan first):
  mod2 <- fromlavaan(' visual  =~ x1 + x2 + x3
                       textual =~ x4 + x5 + x6
                       speed   =~ x7 + x8 + x9 ',
                     data = HS, run = TRUE)
}
}
}
