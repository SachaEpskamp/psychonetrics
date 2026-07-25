\name{plot.psychonetrics}
\alias{plot.psychonetrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the networks of a psychonetrics model
}
\description{
A wrapper around \code{\link[qgraph]{qgraph}} that draws one or more model matrices of a \code{psychonetrics} model as networks, in the spirit of the plot method for \code{estimateNetwork} objects in the \code{bootnet} package: the same overwritten defaults (\code{cut = 0}, \code{theme = "colorblind"}, \code{layout = "spring"}, \code{parallelEdge = TRUE}) and everything else passed on to \code{qgraph} through \code{\dots}. Which matrices are plotted, and how, is derived from the model, so that a plain \code{plot(mod)} gives a sensible picture for every model family.
}
\usage{
\method{plot}{psychonetrics}(x, matrix, \dots, time, group,
     threshold = FALSE, alpha = 0.01,
     adjust = c("none", "holm", "hochberg", "hommel",
                "bonferroni", "BH", "BY", "fdr"),
     mode = c("tested", "all"),
     covariances = FALSE, labels, layout = "spring",
     cut = 0, theme = "colorblind", parallelEdge = TRUE,
     titles = TRUE, equalLayout = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{matrix}{
A character vector of matrices to plot; each element becomes a panel, in the given order. Unlike \code{\link{getmatrix}}, more than one matrix may be given. When omitted, a default is chosen from the model family (see Details).
}
  \item{\dots}{
Arguments passed to \code{\link[qgraph]{qgraph}} (or to \code{\link[qgraph]{qgraph.loadings}} for a factor loading matrix, or to \code{semPaths} for a latent variable model). Anything given here overrides the defaults set by this method.
}
  \item{time}{
Only for \code{\link{ri_clpm}} / \code{ri_clpn} models: the wave whose block should be plotted, either an integer or a column name of the design matrix. When omitted, every wave gets its own panel. See \code{\link{getmatrix}}.
}
  \item{group}{
The group(s) to plot, as labels or numbers. Defaults to all groups, one row of panels per group.
}
\item{threshold}{Logical (or a numeric absolute threshold), as in \code{\link{getmatrix}}: should non-significant edges be set to zero?}
  \item{alpha}{
Significance level used when \code{threshold = TRUE}.
}
  \item{adjust}{
p-value adjustment method used when \code{threshold = TRUE}. See \code{p.adjust}.
}
\item{mode}{
Mode for adjusting for multiple comparisons, as in \code{\link{getmatrix}}.
}
\item{covariances}{
Logical. By default, covariance matrices are drawn as correlation networks and precision matrices as partial correlation networks, so that the edge widths are comparable across differently scaled variables. Set to \code{TRUE} to draw the unstandardized matrix instead, with self-loops showing the diagonal (the variances, respectively the diagonal of the precision matrix).
}
\item{labels}{
Node labels. Defaults to the variable names of the model. May be a character vector (used for every panel) or a list with one element per panel.
}
\item{layout}{
The \code{qgraph} layout. Defaults to \code{"spring"}.
}
\item{cut}{
The \code{qgraph} \code{cut} value. Defaults to \code{0} (no cut), as in \code{bootnet}.
}
\item{theme}{
The \code{qgraph} theme. Defaults to \code{"colorblind"}.
}
\item{parallelEdge}{
Should reciprocal directed edges be drawn as parallel lines? Defaults to \code{TRUE}.
}
\item{titles}{
\code{TRUE} (default) to title each panel with what it shows, \code{FALSE} for no titles, or a character vector with one title per panel.
}
\item{equalLayout}{
Logical. When several panels share the same nodes, should they be constrained to a common layout (through \code{\link[qgraph]{averageLayout}}) so that they are visually comparable? Defaults to \code{TRUE}. Ignored when \code{layout} is given as a matrix.
}
\item{verbose}{
Logical, should messages about how a matrix is being drawn be printed?
}
}
\details{
\bold{Default matrices.} When \code{matrix} is omitted, the matrices are chosen from the model family, in the parameterization that was estimated:

\itemize{
\item \code{\link{varcov}} family (\code{\link{ggm}}, \code{precision}, \code{corr}, ...): the (co)variance matrix itself, one panel.
\item \code{\link{Ising}} and \code{BlumeCapel}: the network \code{omega}, one panel.
\item \code{\link{var1}} / \code{gvar}: two panels, the temporal and the contemporaneous network.
\item \code{\link{panelvar}} / \code{\link{panelgvar}}, \code{\link{dlvm1}}, \code{ml_var1} and \code{\link{ml_lvm}}: three panels, the temporal, the contemporaneous (within-person) and the between-person network.
\item \code{\link{tsdlvm1}}: the temporal and the contemporaneous network.
\item \code{ml_varcov}: the within-person and the between-person network.
\item \code{meta_varcov}: the pooled (co)variance matrix.
\item \code{\link{ri_clpm}} / \code{ri_clpn}: the temporal, contemporaneous and between-person \emph{blocks} (see \code{\link{getmatrix}}); one panel per wave, with waves that were equated by \code{\link{ri_clpm_stationary}} collapsed into one.
\item \code{\link{lnm}}, \code{rnm} and \code{lrnm}: the latent and/or residual network.
\item Any other \code{\link{lvm}} model has no natural network representation, so the full path diagram is drawn with \code{semPaths} from the \code{semPlot} package (via \code{\link{tolavaan}}). Give a \code{matrix} to plot a model matrix as a network instead.
}

If the model uses the \code{"PDC"} temporal parameterization, the temporal panel shows \code{PDC} rather than \code{beta}.

\bold{How each matrix is drawn.} The type of graph follows from the matrix:

\itemize{
\item \code{omega} (Gaussian graphical model): an undirected network of partial correlations.
\item \code{kappa} (precision): an undirected network. Standardized to partial correlations unless \code{covariances = TRUE}.
\item \code{sigma} (covariance): a network with two arrowheads per edge. Standardized to correlations unless \code{covariances = TRUE}, in which case the variances are added as self-loops.
\item \code{rho} (correlation): a network with two arrowheads per edge.
\item \code{beta} and \code{PDC} (temporal): a directed network, with the autoregressive effects drawn as self-loops.
\item \code{lambda} (factor loadings): a reflective measurement diagram, drawn with \code{\link[qgraph]{qgraph.loadings}}.
\item \code{lowertri} and \code{exo_cholesky} (Cholesky decompositions) have no network interpretation and raise an error; plot the corresponding (co)variance or partial correlation matrix instead.
\item Diagonal scaling matrices (\code{delta}, \code{SD}), mean vectors (\code{mu}, \code{nu}, \code{tau}) and the Ising inverse temperature (\code{beta}) are not networks and raise an error.
}
}
\value{
Invisibly, the \code{qgraph} object of the plot, or a named list of such objects when more than one panel is drawn.
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{getmatrix}}, \code{\link{CIplot}}, \code{\link[qgraph]{qgraph}}
}
\examples{
\donttest{
# Simulate data from a chain graph:
set.seed(1)
omega <- matrix(0, 5, 5)
for (i in 1:4) omega[i, i + 1] <- omega[i + 1, i] <- 0.35
Sigma <- cov2cor(solve(diag(5) - omega))
data  <- as.data.frame(MASS::mvrnorm(500, rep(0, 5), Sigma))
names(data) <- paste0("V", 1:5)

# Fit a GGM and plot it:
mod <- runmodel(ggm(data))
plot(mod)

# Only the significant edges:
plot(mod, threshold = TRUE)

# The same data as a covariance model: drawn as a correlation network ...
mod2 <- runmodel(varcov(data, type = "cov"))
plot(mod2)
# ... or as covariances, with the variances as self-loops:
plot(mod2, covariances = TRUE)

# qgraph arguments are passed through:
plot(mod, layout = "circle", vsize = 12, labels = paste0("x", 1:5))
}
}
