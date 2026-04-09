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

The algorithm used by \code{method = "jacobian"} (the default) is a direct
port of the score-test logic in \code{lavaan::lavTestScore} (Rosseel, 2012)
and reproduces \pkg{lavaan}'s output to numerical precision. The derivation
closely follows the implementation in \code{R/lav_test_score.R} of the
\pkg{lavaan} source code (\url{https://github.com/yrosseel/lavaan}), reused
here with thanks to Yves Rosseel under the GPL (\eqn{\geq} 2) license shared
by both packages. The underlying gradient and Fisher-information machinery
is provided by \pkg{psychonetrics}, so no refitting and no run-time
dependence on \pkg{lavaan} is required.

With the default \code{method = "jacobian"} the algorithm reproduces
\code{lavaan::lavTestScore} exactly: the model is augmented with G-1
extra parameters per equality-constrained tuple (one per non-reference
group), the gradient and expected information are computed at the
constrained MLE, and for each constraint the test statistic
\deqn{T = N \cdot s' \cdot \tilde{J}^{-1} \cdot s}
is formed, where \eqn{s} is the gradient of the log-likelihood and
\eqn{\tilde{J}^{-1}} is the upper-left block of the generalised inverse
of the bordered information matrix \eqn{[J\ R'^T;\ R\ 0]}, with \eqn{R}
the rows of the constraint Jacobian for the constraints that are NOT
being tested. This is the same construction used by
\code{lavaan::lavTestScore}.

\code{method = "schur"} uses a faster Schur-complement construction
(the same used by \pkg{psychonetrics} for ordinary modification
indices). It is much faster for complex models because it inverts only
the (smaller) block of the information matrix corresponding to the
currently free parameters, but it is not exactly equivalent to
\code{lavaan} in unbalanced multi-group samples (it agrees in
balanced samples). Both methods follow the same asymptotic
chi-square reference distribution.
}
\usage{
equalityScoreTest(x, matrices, top = 10, verbose = TRUE,
                  method = c("jacobian","schur"))
}
\arguments{
  \item{x}{A fitted (run) multi-group \code{psychonetrics} model.}
  \item{matrices}{Optional vector of matrix names to restrict the test to. By
    default all matrices are considered.}
  \item{top}{How many rows to print in each table.}
  \item{verbose}{Logical, should the formatted tables be printed?}
  \item{method}{Construction of the score test. \code{"jacobian"}
    (default) reproduces \code{lavaan::lavTestScore} exactly.
    \code{"schur"} uses a Schur-complement construction that is faster
    for complex models but only matches \code{lavaan} for balanced
    multi-group samples.}
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

lavaan source code: \url{https://github.com/yrosseel/lavaan} (in particular
\code{R/lav_test_score.R}).
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
Sacha Epskamp <mail@sachaepskamp.com>. The \code{method = "jacobian"}
implementation is a port of the score-test machinery in
\code{lavaan::lavTestScore} by Yves Rosseel.
}
