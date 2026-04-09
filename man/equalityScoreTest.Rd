\name{equalityScoreTest}
\alias{equalityScoreTest}
\title{
Score (Lagrange Multiplier) test for cross-group equality constraints
}
\description{
Computes the score test (also known as the Lagrange Multiplier test) for releasing
cross-group equality constraints in a multi-group \code{psychonetrics} model. Both
the per-group univariate score statistics and the joint multivariate score
statistic (one statistic per equality-constrained parameter, with df = G - 1)
are returned.

This function is inspired by \code{lavaan::lavTestScore} from the \pkg{lavaan}
package and is intended to produce comparable output without depending on
\pkg{lavaan}. The underlying machinery reuses the Fisher-information
infrastructure already used by \code{addMIs} for ordinary modification
indices, so no refitting of the model is required.

The reported statistics are scaled (by a factor of 2) to be on the same
\eqn{-2\log L} scale as \code{lavaan::lavTestScore}, so for balanced
multi-group models the values match \code{lavaan} exactly. Two minor
sources of discrepancy can remain: (i) \pkg{psychonetrics} can use either
the maximum-likelihood (\eqn{1/n}) or unbiased (\eqn{1/(n-1)}) sample
covariance via the \code{covtype} argument of e.g. \code{\link{varcov}},
while \pkg{lavaan} always uses the unbiased version, and (ii)
\pkg{psychonetrics} constructs the score test from the Schur complement
of the augmented Fisher information (the same construction used for
ordinary modification indices), while \pkg{lavaan} uses a
constraint-Jacobian construction. The two formulations are asymptotically
equivalent and follow the same chi-square reference distribution.
}
\usage{
equalityScoreTest(x, matrices, top = 10, verbose = TRUE)
}
\arguments{
  \item{x}{A fitted (run) multi-group \code{psychonetrics} model.}
  \item{matrices}{Optional vector of matrix names to restrict the test to. By
    default all matrices are considered.}
  \item{top}{How many rows to print in each table.}
  \item{verbose}{Logical, should the formatted tables be printed?}
}
\value{
Invisibly returns a list with two data frames:
\describe{
  \item{\code{uni}}{Univariate score tests, one row per released parameter
    (i.e. one row per group beyond the reference group), sorted from largest
    to smallest \eqn{X^2}.}
  \item{\code{total}}{Joint multivariate score tests, one row per
    equality-constrained parameter (i.e. one row per (matrix, row, column) triple),
    sorted from largest to smallest \eqn{X^2}. The joint test has
    \code{df = G - 1} when the parameter is constrained equal across all
    \eqn{G} groups.}
}
}
\details{
For two groups (\eqn{G = 2}) the joint test reduces to the univariate test
(both have \code{df = 1}). The two tests differ for \eqn{G \ge 3}, where the
joint test correctly accounts for the covariance between the released
parameters via the (Schur-complemented) Fisher information, instead of summing
or otherwise aggregating univariate statistics.
}
\seealso{
\code{\link{MIs}}, \code{\link{partialprune}}
}
\references{
Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling.
Journal of Statistical Software, 48(2), 1-36. \doi{10.18637/jss.v048.i02}
}
\examples{
\donttest{
library("dplyr")
library("psychTools")
data(bfi)

ExData <- bfi \%>\% select(E1:E5, gender) \%>\% na.omit
vars <- names(ExData)[1:5]

mod <- ggm(ExData, vars = vars, groups = "gender") \%>\%
  runmodel \%>\%
  groupequal("omega") \%>\%
  runmodel

equalityScoreTest(mod)
}
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
