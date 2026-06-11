\name{allequal}
\alias{allequal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constrain all elements of a matrix equal within groups
}
\description{
The \code{allequal} function constrains all free elements of a model matrix to be equal to one another \emph{within} each group (a single shared parameter per group). In contrast to \code{\link{groupequal}}, which constrains corresponding elements equal \emph{across} groups, \code{allequal} shares one parameter across the elements of the matrix and does so separately for each group, so the constraint is not imposed across groups. This is convenient for, for example, constraining all thresholds (\code{tau}) or all Blume-Capel quadratic potentials (\code{delta}) to a common value, or fitting an equal-edge-weight network (\code{omega}).
}
\usage{
allequal(x, matrix, row, col, group, verbose, log = TRUE,
         runmodel = FALSE, identify = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{
A \code{psychonetrics} model.
}
  \item{matrix}{
String indicating the matrix whose elements should be constrained equal.
}
  \item{row}{
Optional integer or string vector indicating which rows of the matrix to include. Defaults to all rows.
}
  \item{col}{
Optional integer or string vector indicating which columns of the matrix to include. Defaults to all columns.
}
  \item{group}{
Optional vector of group ids to which the constraint is applied. Defaults to all groups (each constrained separately).
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{log}{
Logical, should the log be updated?
}
  \item{runmodel}{
Logical, should the model be updated?
}
  \item{identify}{
Logical, should the model be identified?
}
 \item{\dots}{
Arguments sent to \code{runmodel}
}
}
\details{
Only free parameters are affected; elements that are fixed (e.g. fixed to zero) are left unchanged. For symmetric matrices all free elements, including a free diagonal, are constrained equal; if this is not desired (for example to keep the diagonal separate from the off-diagonal elements), use the \code{row} and \code{col} arguments to restrict the constraint to the relevant elements. Within each group the selected free elements are given the same parameter number and a common starting value (the mean of their current values).
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{groupequal}}, \code{\link{parequal}}, \code{\link{fixpar}}
}
\examples{
\dontrun{
library("dplyr")

# Simulate a small Blume-Capel data set (requires the development version of
# IsingSampler that supports the 'delta' argument):
library("IsingSampler")
nNode  <- 5
graph  <- 0.4 * (matrix(1, nNode, nNode) - diag(nNode)) *
            (matrix(c(0,1,0,0,1, 1,0,1,0,0, 0,1,0,1,0, 0,0,1,0,1, 1,0,0,1,0), 5, 5))
data <- IsingSampler(1000, graph = graph, thresholds = rep(0, nNode),
                     beta = 1, delta = 1, responses = c(-1, 0, 1))

# Blume-Capel model with all quadratic (delta) potentials constrained equal,
# and all thresholds (tau) constrained equal, within the group:
mod <- BlumeCapel(data, responses = c(-1, 0, 1)) \%>\%
  allequal("delta") \%>\%
  allequal("tau") \%>\%
  runmodel
}
}
