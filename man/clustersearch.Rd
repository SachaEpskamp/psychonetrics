\name{clustersearch}
\alias{clustersearch}
\alias{print.clustersearch}
\alias{plot.clustersearch}
\title{
Search for clusters of variables
}
\description{
Agglomerative search for clusters of variables, in which every cluster is modeled as a composite (\code{\link{cvm}} family) or as a latent variable (\code{\link{lvm}} family). The search starts with the saturated model, in which all variables are freely related. In every step, each pair of clusters (or variables that are not yet clustered) is combined into one cluster, and the combination with the lowest value of the criterion (by default the BIC) becomes the new model. This is repeated until few clusters are left. The result contains the model of every step and the model with the lowest value of the criterion, and can be plotted as a dendrogram next to the criterion. This function is new in version 0.18 and experimental (a message is shown once per session): arguments and defaults may still change.
}
\usage{
clustersearch(data, vars, cluster = c("composite", "latent"),
              type = c("cov", "ggm"), criterion = "bic", min_clusters,
              improper = FALSE, nCores = 1, verbose = TRUE, ...)

\method{print}{clustersearch}(x, digits = 2, ...)

\method{plot}{clustersearch}(x, labels, cex = 1, nodelabels = TRUE,
     col = "black", grid = "grey90", ...)
}
\arguments{
  \item{data}{
A data frame encoding the data used in the analysis. Can be missing if \code{covs} and \code{nobs} are supplied.
}
  \item{vars}{
An optional character vector with the variables to cluster. Defaults to all variables in \code{data}.
}
  \item{cluster}{
How clusters are modeled: \code{"composite"} (default) to model every cluster as a composite of its variables (using \code{\link{cvm}}), or \code{"latent"} to model every cluster as a latent variable that is measured by its variables (using \code{\link{lvm}}). All clusters in a search are of the same kind.
}
  \item{type}{
The type of the models that are returned: \code{"cov"} (default) for variance--covariance matrices, or \code{"ggm"} for Gaussian graphical models. With \code{"ggm"}, a search for composites returns composite network models (see \code{\link{cnm}}: a network of the variables within every composite, \code{omega_c}, and a network of the clusters, \code{omega_zeta}), and a search for latent variables returns latent network models (see \code{\link{lnm}}: a network of the clusters, \code{omega_zeta}). The two types give equivalent models, so the type does not change which clusters are found.
}
  \item{criterion}{
The fit measure that is minimized: \code{"bic"} (default), \code{"aic"}, or one of \code{"aic.ll"}, \code{"aic.ll2"}, \code{"bic2"}, \code{"ebic.25"}, \code{"ebic.5"}, \code{"ebic.75"} and \code{"ebic1"} (see \code{\link{fit}}).
}
  \item{min_clusters}{
The number of clusters at which the search stops (variables that are not part of a cluster count as a cluster). Defaults to 2 in a search for composites and to 1 in a search for latent variables. See Details.
}
  \item{improper}{
Logical, can a model with an improper solution be selected in a search for latent variables? Defaults to \code{FALSE}. See Details.
}
  \item{nCores}{
The number of cores to use. With more than one core, the candidate models of a step are estimated in parallel.
}
  \item{verbose}{
Logical, should progress be printed to the console?
}
  \item{\dots}{
For \code{clustersearch}: arguments sent to \code{\link{cvm}} or \code{\link{lvm}} on the data (\code{groupvar}, \code{covs}, \code{means}, \code{nobs}, \code{missing}, \code{covtype}, \code{standardize}, \code{storedata}), the estimator (\code{estimator}, which can be \code{"ML"} or \code{"FIML"}), the optimizer (\code{optimizer}), the identification (\code{identification}) and equality constraints across groups (\code{equal}). In a search for composites also \code{bad_marker_crit} and \code{sumscores}. For the \code{print} and \code{plot} methods: not used.
}
  \item{x}{
The result of \code{clustersearch}.
}
  \item{digits}{
Number of digits to print.
}
  \item{labels}{
Optional labels for the variables in the plot.
}
  \item{cex}{
Scaling of the text and symbols in the plot.
}
  \item{nodelabels}{
Logical, should the names of the clusters be shown in the plot?
}
  \item{col}{
Color of the dendrogram and the criterion.
}
  \item{grid}{
Color of the grid lines.
}
}
\details{
\strong{The search.} In the starting model every variable is its own construct, and all constructs are freely related (the saturated model). A step considers every pair of constructs. Combining two constructs gives one cluster of all their variables: combining the cluster of the variables A and B with the variable C gives a cluster with A, B and C as its variables. The models in which the clusters are not combined remain as they are: all clusters and remaining variables are freely related to each other. The combination with the lowest value of the criterion is kept, also when the criterion is higher than in the previous step, and the search continues until \code{min_clusters} clusters are left. The selected model is the model with the lowest value of the criterion over all steps. With \eqn{p}{p} variables the first step estimates \eqn{p(p-1)/2}{p * (p - 1) / 2} candidate models, and the number of candidate models over all steps grows with the third power of the number of variables. The search is greedy: a combination is never undone, and a cluster of three or more variables can only be found by first combining two of its variables.

\strong{Composites.} With \code{cluster = "composite"} a cluster is a composite: a weighted sum through which its variables relate to the other clusters and variables, while the covariances between the variables of the same cluster are not restricted (see \code{\link{cvm}}). A variable that is not part of a cluster is a composite with a single indicator, which represents that variable itself. Because the covariances within a composite are not restricted, combining clusters into a larger composite can increase the number of parameters, and a single composite of all variables is the saturated model again. The search for composites therefore stops at two clusters by default. With \code{sumscores = TRUE} every cluster is the sum of its variables.

\strong{Latent variables.} With \code{cluster = "latent"} a cluster is a latent variable that is measured by its variables with uncorrelated residuals (see \code{\link{lvm}}). A variable that is not part of a cluster is a latent variable with a single indicator and no residual variance. The latent variables are by default identified by fixing their variances to 1 (\code{identification = "variance"}), such that a cluster does not depend on the loading of its first variable. A candidate model with an improper solution (a residual variance on its lower bound of zero, or a variance--covariance matrix of the clusters that is not positive definite) is not selected, unless \code{improper = TRUE}. The search continues to a single latent variable by default.

\strong{Two variables.} For a cluster of two variables a composite and a latent variable (with a proper solution) give the same fit, so that both kinds of search evaluate pairs of variables in the same way. They differ for clusters of three or more variables, for which a latent variable is the more restrictive model.

\strong{Candidate models that are skipped.} A candidate model is not selected when it gives an error, when the optimizer does not report convergence, or (in a search for latent variables) when its solution is improper. The \code{candidates} element of the output lists all candidate models with the reason. If no candidate model can be selected in a step, the search stops with a warning.

\strong{The models.} All models are regular models of the \code{\link{cvm}} or \code{\link{lvm}} family, estimated with \code{\link{runmodel}}. The cluster that is formed in step \eqn{s}{s} is named \code{C}\eqn{s}{s} in the models, and a variable that is not part of a cluster is a construct with the name of that variable. Because such a construct has the same name as its indicator, the models cannot be converted with \code{\link{tolavaan}}.
}
\value{
An object of the class \code{clustersearch}, which is a list with the following elements:
\item{model}{The model with the lowest value of the criterion (an object of the class psychonetrics, \link{psychonetrics-class}).}
\item{models}{A list with the model of every step, named by the number of clusters (for example, \code{models[["3"]]} is the model with three clusters).}
\item{path}{A data frame with one row per step: the number of clusters, the name of the new cluster (\code{cluster}), the two constructs that are combined (\code{cluster1} and \code{cluster2}), the number of variables in the new cluster (\code{size}), the chi-square, degrees of freedom and number of parameters of the model, and the criterion.}
\item{candidates}{A data frame with all candidate models: the step, the two constructs that are combined, the chi-square, degrees of freedom, number of parameters and criterion, whether the optimizer reported convergence, whether the solution is improper, the error message (if any), whether the candidate could be selected (\code{eligible}) and whether it was selected.}
\item{clusters}{A named list with the variables of every cluster (and of every variable that is not part of a cluster) in the selected model.}
\item{partitions}{A list with the clusters of the model of every step, named by the number of clusters.}
\item{warnings}{A list with the warnings given when the model of every step was estimated.}
\item{best}{The step of the selected model.}
\item{criterion, cluster, type, vars}{The settings of the search.}

The \code{plot} method draws the dendrogram of the search above the criterion, with the steps of the search on the shared horizontal axis. A dashed line marks the selected model. It invisibly returns the positions of the nodes in the dendrogram.
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{cvm}}, \code{\link{cnm}}, \code{\link{lvm}}, \code{\link{lnm}}
}
\examples{
\donttest{
# Load data:
data(NA2020)

# Search for composites (the default):
res <- clustersearch(NA2020)
res

# Dendrogram and BIC:
plot(res)

# The selected model and its clusters:
res$clusters
res$model

# Another model, for example with four clusters:
res$models[["4"]]

# Search for latent variables:
res_latent <- clustersearch(NA2020, cluster = "latent")
res_latent
plot(res_latent)
}
}
