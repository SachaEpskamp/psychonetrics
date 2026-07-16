\name{psychonetrics-class}
\Rdversion{1.1}
\docType{class}
\alias{psychonetrics-class}
\alias{resid,psychonetrics-method}
\alias{residuals,psychonetrics-method}
\alias{show,psychonetrics-method}
\alias{summary,psychonetrics-method}
\alias{psychonetrics_samplestats-class}
\alias{sessionInfo-class}


\title{Class \code{"psychonetrics"}}
\description{
Main class for psychonetrics results.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("psychonetrics", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"character"}: the model framework (e.g. \code{"varcov"}, \code{"lvm"}, \code{"var1"}).}
    \item{\code{submodel}:}{Object of class \code{"character"}: the submodel/type within the framework (e.g. \code{"ggm"}).}
    \item{\code{parameters}:}{Object of class \code{"data.frame"}: the parameter table (one row per model parameter, with estimates, standard errors, matrix/row/col indices and parameter numbers).}
    \item{\code{matrices}:}{Object of class \code{"data.frame"}: bookkeeping table of the model matrices and their dimensions.}
    \item{\code{meanstructure}:}{Object of class \code{"logical"}: whether a mean structure is included in the model.}
    \item{\code{computed}:}{Object of class \code{"logical"}: whether the model has been estimated (run).}
    \item{\code{sample}:}{Object of class \code{"psychonetrics_samplestats"}: the sample statistics used in estimation.}
    \item{\code{modelmatrices}:}{Object of class \code{"list"}: the model-implied matrices in list form (per group).}
    \item{\code{log}:}{Object of class \code{"psychonetrics_log"}: the logbook of actions performed on the model.}
    \item{\code{optim}:}{Object of class \code{"list"}: the raw output of the optimizer.}
    \item{\code{fitmeasures}:}{Object of class \code{"list"}: computed fit measures.}
    \item{\code{baseline_saturated}:}{Object of class \code{"list"}: the baseline and saturated models used for fit computation.}
    \item{\code{equal}:}{Object of class \code{"character"}: names of matrices constrained equal across groups (legacy representation; the parameters table is authoritative).}
    \item{\code{objective}:}{Object of class \code{"numeric"}: the value of the fit function at the solution.}
    \item{\code{information}:}{Object of class \code{"matrix"}: the Fisher information matrix.}
    \item{\code{identification}:}{Object of class \code{"character"}: the identification method used (e.g. \code{"loadings"} or \code{"variance"}).}
    \item{\code{optimizer}:}{Object of class \code{"character"}: the optimizer used.}
    \item{\code{optim.args}:}{Object of class \code{"list"}: arguments passed to the optimizer.}
    \item{\code{estimator}:}{Object of class \code{"character"}: the estimator used (e.g. \code{"ML"}, \code{"FIML"}).}
    \item{\code{distribution}:}{Object of class \code{"character"}: the assumed distribution (e.g. \code{"Gaussian"}).}
    \item{\code{extramatrices}:}{Object of class \code{"list"}: additional helper matrices used internally.}
    \item{\code{rawts}:}{Object of class \code{"logical"}: whether raw time-series data are used.}
    \item{\code{Drawts}:}{Object of class \code{"list"}: duplication/derivative matrices for raw time-series estimation.}
    \item{\code{types}:}{Object of class \code{"list"}: the type of each model matrix per group.}
    \item{\code{cpp}:}{Object of class \code{"logical"}: whether the C++ backend is used.}
    \item{\code{verbose}:}{Object of class \code{"logical"}: whether progress messages are printed.}
    \item{\code{penalty}:}{Object of class \code{"list"}: the penalty configuration for penalized (PML/PFIML) estimation, holding the penalty strength (\code{lambda}) and elastic-net mixing parameter (\code{alpha}).}
  }
}
\section{Methods}{
  \describe{
    \item{resid}{\code{signature(object = "psychonetrics")}: ... }
    \item{residuals}{\code{signature(object = "psychonetrics")}: ... }
    \item{show}{\code{signature(object = "psychonetrics")}: ... }
    \item{summary}{\code{signature(object = "psychonetrics")}: The print output extended with the main fit measures (log-likelihood, number of parameters, AIC, BIC, RMSEA, CFI, TLI). Use \code{\link{fit}} for the complete list of fit measures. }
	 }
}
\author{
Sacha Epskamp
}

\examples{
showClass("psychonetrics")
}
\keyword{classes}
