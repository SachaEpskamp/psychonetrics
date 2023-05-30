\name{esa}
\alias{esa}
\alias{esa_manual}
\alias{plot.esa}
\alias{plot.esa_manual}
\alias{print.esa}
\alias{print.esa_manual}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ergodic Subspace Analysis
}
\description{
These functions implement Ergodic Subspace Analysis by von Oertzen, Schmiedek and Voelkle (2020). The functions can be used on the output of a \code{\link{dlvm1}} model, or manually by supplying a within persons and between persons variance-covariance matrix.
}
\usage{
esa(x, cutoff = 0.1,
    between = c("crosssection", "between"))
esa_manual(sigma_wp, sigma_bp, cutoff = 0.1)
\method{print}{esa}(x, printref = TRUE, ...)
\method{print}{esa_manual}(x, printref = TRUE, ...)
\method{plot}{esa}(x, plot = c("observed", "latent"), ...)
\method{plot}{esa_manual}(x,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Output of a \code{\link{dlvm1}} model
}
\item{sigma_wp}{
Manual within-person variance-covariance matrix
}
\item{sigma_bp}{
Manual between-person variance-covariance matrix
}
  \item{cutoff}{
Cutoff used to determine ergodicity
}
 \item{printref}{
Logical, should the reference be printed?
}
\item{plot}{
Should ergodicity of observed or latent variables be plotted?
}
\item{between}{
Should the between-persons variance-covariance matrix be based on exected cross-sectional or between-person relations
}
 \item{\dots}{
Not used
}
}

\value{
For each group a \code{esa_manual} object with the following elements:
\item{ergodicity }{Ergodicity values of each component}
\item{Q_esa}{Component loadings}
\item{V_bp}{Between persons subspace}
\item{V_ergodic}{Ergodic subspace}
\item{V_wp}{Within person subspace}
\item{cutoff}{Cutoff value used}
%% ...
}
\references{
von Oertzen, T., Schmiedek, F., and Voelkle, M. C. (2020). Ergodic Subspace Analysis. Journal of Intelligence, 8(1), 3.
}
\author{
Sacha Epskamp <mail@sachaepskamp.com>
}

