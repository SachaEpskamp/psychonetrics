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

setoptimizer(x, optimizer = c("default","nlminb","ucminf","cpp_L-BFGS-B",
            "cpp_BFGS","cpp_CG","cpp_SANN","cpp_Nelder-Mead"))

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
}
\value{
An object of the class psychonetrics (\link{psychonetrics-class})
}
\author{
Sacha Epskamp
}

