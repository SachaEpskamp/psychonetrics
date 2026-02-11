\name{loop_psychonetrics}
\alias{loop_psychonetrics}
\title{
Parallel loop for psychonetrics
}
\description{
A convenience wrapper for running repeated psychonetrics analyses in parallel. Replaces the typical boilerplate of \code{makeCluster} / \code{clusterExport} / \code{parLapply} / \code{stopCluster} with a single function call. Particularly useful for bootstrapping, simulation studies, and cross-validation.

Variables referenced in the expression that exist in the calling environment are automatically exported to parallel workers, so explicit \code{clusterExport} calls are not needed in most cases.
}
\usage{
loop_psychonetrics(expr, reps = 1000, nCores = 1,
                   export = character(0),
                   packages = c("psychonetrics", "dplyr"),
                   verbose = TRUE, envir = parent.frame())
}
\arguments{
  \item{expr}{An expression (code block) to evaluate on each iteration. Typically enclosed in curly braces \code{\{...\}}. The expression should return a single object (e.g., a psychonetrics model).}
  \item{reps}{Number of replications (default: 1000).}
  \item{nCores}{Number of CPU cores to use. Set to 1 for sequential execution (default), or higher for parallel execution via \code{makePSOCKcluster}.}
  \item{export}{Character vector of additional variable names to export to parallel workers. Usually not needed because variables referenced in \code{expr} are auto-detected and exported.}
  \item{packages}{Character vector of packages to load on parallel workers (default: \code{c("psychonetrics", "dplyr")}).}
  \item{verbose}{Logical; if \code{TRUE} (default), shows a progress bar via \code{pbapply}.}
  \item{envir}{The environment in which to look up variables for auto-export and to evaluate the expression in sequential mode. Defaults to the calling environment.}
}
\value{
A list of length \code{reps}, where each element is the result of evaluating \code{expr} once.
}
\details{
When \code{nCores > 1}, a PSOCK cluster is created and automatically cleaned up when the function exits (even on error). The specified \code{packages} are loaded on each worker, and variables referenced in \code{expr} are automatically detected and exported.

For bootstrapping, use this function together with \code{bootstrap = "nonparametric"} in the model constructor and \code{\link{aggregate_bootstraps}} to summarize the results.
}
\seealso{
\code{\link{aggregate_bootstraps}}, \code{\link{bootstrap}}, \code{\link{CIplot}}
}
\examples{
\donttest{
library(dplyr)
data(NA2020)

# Fit the original model:
model <- ggm(NA2020, estimator = "FIML") \%>\% runmodel

# Bootstrap with pruning (parallel, 100 reps):
bootstraps <- loop_psychonetrics({
  ggm(NA2020, bootstrap = "nonparametric",
      estimator = "FIML") \%>\%
    runmodel \%>\%
    prune(alpha = 0.05)
}, reps = 100, nCores = 2)

# Aggregate and inspect:
boot_agg <- aggregate_bootstraps(
  sample = model \%>\% prune(alpha = 0.05),
  bootstraps = bootstraps
)

# View results:
parameters(boot_agg)
CIplot(boot_agg, "omega", split0 = TRUE)
}
}
\keyword{utilities}
