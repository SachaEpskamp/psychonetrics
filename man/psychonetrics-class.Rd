\name{psychonetrics-class}
\Rdversion{1.1}
\docType{class}
\alias{psychonetrics-class}
\alias{resid,psychonetrics-method}
\alias{residuals,psychonetrics-method}
\alias{show,psychonetrics-method}

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
    \item{\code{model}:}{Object of class \code{"character"} ~~ }
    \item{\code{submodel}:}{Object of class \code{"character"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{matrices}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{meanstructure}:}{Object of class \code{"logical"} ~~ }
    \item{\code{computed}:}{Object of class \code{"logical"} ~~ }
    \item{\code{sample}:}{Object of class \code{"psychonetrics_samplestats"} ~~ }
    \item{\code{modelmatrices}:}{Object of class \code{"list"} ~~ }
    \item{\code{log}:}{Object of class \code{"psychonetrics_log"} ~~ }
    \item{\code{optim}:}{Object of class \code{"list"} ~~ }
    \item{\code{fitmeasures}:}{Object of class \code{"list"} ~~ }
    \item{\code{baseline_saturated}:}{Object of class \code{"list"} ~~ }
    \item{\code{equal}:}{Object of class \code{"character"} ~~ }
    \item{\code{objective}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{information}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{identification}:}{Object of class \code{"character"} ~~ }
    \item{\code{optimizer}:}{Object of class \code{"character"} ~~ }
    \item{\code{estimator}:}{Object of class \code{"character"} ~~ }
    \item{\code{distribution}:}{Object of class \code{"character"} ~~ }
    \item{\code{extramatrices}:}{Object of class \code{"list"} ~~ }
    \item{\code{rawts}:}{Object of class \code{"logical"} ~~ }
    \item{\code{Drawts}:}{Object of class \code{"list"} ~~ }
    \item{\code{types}:}{Object of class \code{"list"} ~~ }
    \item{\code{cpp}:}{Object of class \code{"logical"} ~~ }
    \item{\code{verbose}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{resid}{\code{signature(object = "psychonetrics")}: ... }
    \item{residuals}{\code{signature(object = "psychonetrics")}: ... }
    \item{show}{\code{signature(object = "psychonetrics")}: ... }
	 }
}
\author{
Sacha Epskamp
}

\examples{
showClass("psychonetrics")
}
\keyword{classes}
