\name{setestimator}
\alias{setestimator}
\alias{setoptimizer}
\alias{usecpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convenience functions
}
\description{
These functions can be used to change some estimator options.
}
\usage{
setestimator(x, estimator)

setoptimizer(x, optimizer = c("default", "nlminb", "ucminf",
                 "cpp_L-BFGS-B", "cpp_BFGS", "cpp_CG", "cpp_SANN",
                 "cpp_Nelder-Mead"), optim.args)


usecpp(x, use = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{estimator}{
A string indicating the estimator to be used
}
  \item{optimizer}{
The optimizer to be used. Can be one of \code{"nlminb"} (the default R \code{nlminb} function), \code{"ucminf"} (from the \code{optimr} package), and C++ based optimizers \code{"cpp_L-BFGS-B"}, \code{"cpp_BFGS"}, \code{"cpp_CG"}, \code{"cpp_SANN"}, and \code{"cpp_Nelder-Mead"}. The C++ optimizers are faster but slightly less stable. Defaults to \code{"nlminb"}.
}
\item{use}{
Logical indicating if C++ should be used (currently only used in FIML)
}
\item{optim.args}{List of arguments to sent to the optimizer.}
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\details{The default optimizer is nlminb with the following arguments:
\itemize{
\item eval.max=20000L
\item iter.max=10000L
\item trace=0L
\item abs.tol=sqrt(.Machine$double.eps)
\item rel.tol=sqrt(.Machine$double.eps)
\item step.min=1.0
\item step.max=1.0
\item x.tol=1.5e-8
\item xf.tol=2.2e-14
}

}
\author{
Sacha Epskamp
}

