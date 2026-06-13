\name{setestimator}
\alias{setestimator}
\alias{setoptimizer}
\alias{usecpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convenience functions
}
\description{
These functions can be used to change some estimator options.
}
\usage{
setestimator(x, estimator)

setoptimizer(x, optimizer = c("default", "nlminb", "ucminf",
                 "nloptr_TNEWTON", "LBFGS++"), optim.args)


usecpp(x, use = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{estimator}{
A string indicating the estimator to be used. Maximum likelihood (\code{"ML"}),
full-information maximum likelihood (\code{"FIML"}), and the least-squares
estimators (\code{"WLS"}, \code{"DWLS"}, \code{"ULS"}; \code{"WLSMV"} is a
synonym for \code{"DWLS"} with a mean-and-variance adjusted scaled test) are
supported. In addition, the following \emph{robust maximum likelihood}
estimators are available for complete continuous data (raw data required).
They use plain ML for the point estimates but report robust (sandwich) standard
errors, a scaled test statistic and robust RMSEA/CFI/TLI:
\describe{
  \item{\code{"MLM"}}{Browne (1984) "robust.sem" sandwich standard errors and the Satorra-Bentler (1994) mean-adjusted scaled chi-square.}
  \item{\code{"MLMV"}}{Robust.sem standard errors and the scaled-and-shifted chi-square (Satorra 2000; Asparouhov & Muthen 2010).}
  \item{\code{"MLMVS"}}{Robust.sem standard errors and the mean-and-variance adjusted (Satterthwaite) scaled chi-square with fractional degrees of freedom.}
  \item{\code{"MLR"}}{Huber-White (sandwich) standard errors and the Yuan-Bentler (2000) scaled chi-square (the Mplus variant). Requires the raw data to be stored; the model constructors set \code{storedata = TRUE} automatically when \code{estimator = "MLR"}. \code{"MLR"} also supports \emph{missing} data: when the raw data contain missing values it estimates the point estimates by full-information maximum likelihood (FIML) and computes the Huber-White standard errors, the Yuan-Bentler-Mplus scaled chi-square and the FIML-corrected (Savalei 2010) robust RMSEA/CFI/TLI under missingness.}
}
These robust estimators map internally to \code{estimator = "ML"} (or, for
\code{"MLR"} with missing data, \code{estimator = "FIML"}) plus robust
standard-error and test flags, so the point estimates are identical to plain ML
(or FIML). The scaled statistics are stored in the fit measures
(\code{chisq.scaled}, \code{chisq.scaling.factor}, \code{df.scaled},
\code{chisq.shift.parameters}, \code{rmsea.robust}, \code{cfi.robust},
\code{tli.robust}). The Satorra-Bentler family
(\code{"MLM"}/\code{"MLMV"}/\code{"MLMVS"}) requires complete data; for missing
data use \code{"MLR"}.
}
  \item{optimizer}{
The optimizer to be used. The following optimizers are available:

\bold{R-based optimizers:}
\describe{
  \item{\code{"nlminb"}}{Port trust-region Newton-like optimizer from base R. Uses analytic gradients and a trust-region method that adaptively controls step sizes, making it very stable for SEM problems. This is the default and recommended optimizer for most models.}
  \item{\code{"ucminf"}}{Unconstrained minimization using a quasi-Newton method (via the \code{optimr} package). Does not support box constraints. Fast for unconstrained problems.}
}

\bold{C++ based optimizer:}
\describe{
  \item{\code{"LBFGS++"}}{L-BFGS-B optimizer from the LBFGSpp library (Yixuan Qiu). This is a pure C++ implementation that computes the objective function and gradient in a single combined call, avoiding redundant computation via internal caching. Supports box constraints. Recommended when speed is important.}
}

\bold{NLopt-based optimizer (via nloptr):}
\describe{
  \item{\code{"nloptr_TNEWTON"}}{Preconditioned truncated Newton with restarts from the NLopt library (via the \code{nloptr} package). Uses the \code{NLOPT_LD_TNEWTON_PRECOND_RESTART} algorithm, which builds a local quadratic model of the objective and solves the Newton system approximately using a preconditioned conjugate-gradient method. Supports box constraints and uses analytic gradients. Can handle large parameter spaces efficiently. See Dembo & Steihaug (1982) for the underlying method.}
}

Defaults to \code{"nlminb"}.
}
\item{use}{
Logical indicating if C++ should be used (currently only used in FIML)
}
\item{optim.args}{List of arguments to sent to the optimizer.}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\details{The default optimizer is nlminb with the following arguments:
\itemize{
\item eval.max=20000L
\item iter.max=10000L
\item trace=0L
\item abs.tol=sqrt(.Machine$double.eps)
\item rel.tol=sqrt(.Machine$double.eps)
\item step.min=1.0
\item step.max=1.0
\item x.tol=1.5e-8
\item xf.tol=2.2e-14
}

}
\references{
Browne, M. W. (1984). Asymptotically distribution-free methods for the analysis
of covariance structures. \emph{British Journal of Mathematical and Statistical
Psychology, 37}(1), 62-83.

Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
standard errors in covariance structure analysis. In A. von Eye & C. C. Clogg
(Eds.), \emph{Latent variables analysis: Applications for developmental
research} (pp. 399-419). Sage.

Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for mean
and covariance structure analysis with nonnormal missing data. \emph{Sociological
Methodology, 30}(1), 165-200.

Savalei, V. (2010). Expected versus observed information in SEM with incomplete
normal and nonnormal data. \emph{Psychological Methods, 15}(4), 352-367.

Brosseau-Liard, P. E., & Savalei, V. (2014). Adjusting incremental fit indices
for nonnormality. \emph{Multivariate Behavioral Research, 49}(5), 460-470.

Asparouhov, T., & Muthen, B. (2010). \emph{Simple second order chi-square
correction.} Mplus Technical Appendix.
}
\author{
Sacha Epskamp
}
\examples{
\donttest{
# Robust maximum likelihood (MLR) on the Holzinger-Swineford data:
library("lavaan")
data("HolzingerSwineford1939")
vars <- paste0("x", 1:9)

# Three-factor CFA loadings matrix:
Lambda <- matrix(0, 9, 3)
Lambda[1:3, 1] <- Lambda[4:6, 2] <- Lambda[7:9, 3] <- 1

# Robust ML (Huber-White SEs + Yuan-Bentler scaled chi-square). The point
# estimates are identical to plain ML; only the SEs, the scaled test and the
# robust fit indices differ:
mod <- lvm(HolzingerSwineford1939, lambda = Lambda, vars = vars,
           identification = "variance", estimator = "MLR")
mod <- runmodel(mod)

# Robust fit measures (chisq.scaled, rmsea.robust, cfi.robust, ...):
fit(mod)

# Equivalent via setestimator() (the model must store the raw data):
mod2 <- lvm(HolzingerSwineford1939, lambda = Lambda, vars = vars,
            identification = "variance", estimator = "ML", storedata = TRUE)
mod2 <- runmodel(setestimator(mod2, "MLM"))

# MLR also handles missing data via FIML. With missing values in the raw data,
# estimator = "MLR" estimates the point estimates by FIML and reports
# Huber-White standard errors, the Yuan-Bentler-Mplus scaled chi-square and the
# FIML-corrected robust RMSEA/CFI/TLI:
HSmis <- HolzingerSwineford1939
set.seed(1)
for (v in vars) HSmis[[v]][runif(nrow(HSmis)) < 0.1] <- NA
modM <- lvm(HSmis, lambda = Lambda, vars = vars,
            identification = "variance", estimator = "MLR")
modM <- runmodel(modM)   # estimator is "FIML" internally, with MLR robust output
fit(modM)
}
}

