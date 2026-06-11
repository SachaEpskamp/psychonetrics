\name{psychonetrics_log-class}
\Rdversion{1.1}
\docType{class}
\alias{psychonetrics_log-class}
\alias{show,psychonetrics_log-method}

\title{Class \code{"psychonetrics_log"}}
\description{
A logbook entry in the psychonetrics logbook
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("psychonetrics_log", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{event}:}{Object of class \code{"character"}: a description of the logged event.}
    \item{\code{time}:}{Object of class \code{"POSIXct"}: the time at which the event was logged.}
    \item{\code{sessionInfo}:}{Object of class \code{"sessionInfo"}: the R session information at the time of the event.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "psychonetrics_log")}: ... }
	 }
}
\author{
Sacha Epskamp
}

\examples{
showClass("psychonetrics_log")
}
\keyword{classes}
