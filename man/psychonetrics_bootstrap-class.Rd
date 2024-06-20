\name{psychonetrics_bootstrap-class}
\Rdversion{1.1}
\docType{class}
\alias{psychonetrics_bootstrap-class}
\alias{show,psychonetrics_bootstrap-method}

\title{Class \code{"psychonetrics_bootstrap"}}
\description{
Class for aggregated bootstrap results.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("psychonetrics_bootstrap", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"character"} ~~ }
    \item{\code{submodel}:}{Object of class \code{"character"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{models}:}{Object of class \code{"list"} ~~ }
    \item{\code{matrices}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{fitmeasures}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{distribution}:}{Object of class \code{"character"} ~~ }
    \item{\code{verbose}:}{Object of class \code{"logical"} ~~ }
    \item{\code{boot_sub}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{boot_resample}:}{Object of class \code{"logical"} ~~ }
    \item{\code{n_fail}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{n_success}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{types}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "psychonetrics_bootstrap")}: ... }
	 }
}
\author{
Sacha Epskamp
}

\examples{
showClass("psychonetrics_bootstrap")
}
\keyword{classes}
