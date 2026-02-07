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
                 "cpp_L-BFGS-B", "LBFGS++", "cpp_BFGS", "cpp_CG",
                 "cpp_SANN", "cpp_Nelder-Mead"), optim.args)


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
The optimizer to be used. The following optimizers are available:

\bold{R-based optimizers:}
\describe{
  \item{\code{"nlminb"}}{Port trust-region Newton-like optimizer from base R. Uses analytic gradients and a trust-region method that adaptively controls step sizes, making it very stable for SEM problems. This is the default and recommended optimizer for most models.}
  \item{\code{"ucminf"}}{Unconstrained minimization using a quasi-Newton method (via the \code{optimr} package). Does not support box constraints. Fast for unconstrained problems.}
}

\bold{C++ based optimizers (via roptim, wrapping R's \code{optim()} C routines):}
\describe{
  \item{\code{"cpp_L-BFGS-B"}}{Limited-memory BFGS with box constraints. Calls R's internal L-BFGS-B implementation from C++ via the \code{roptim} package.}
  \item{\code{"cpp_BFGS"}}{BFGS quasi-Newton method without box constraints.}
  \item{\code{"cpp_CG"}}{Conjugate gradient method.}
  \item{\code{"cpp_SANN"}}{Simulated annealing. Does not use gradients.}
  \item{\code{"cpp_Nelder-Mead"}}{Nelder-Mead simplex method. Does not use gradients.}
}

\bold{Pure C++ optimizer:}
\describe{
  \item{\code{"LBFGS++"}}{L-BFGS-B optimizer from the LBFGSpp library (Yixuan Qiu). This is a pure C++ implementation that computes the objective function and gradient in a single combined call, avoiding redundant computation via internal caching. Supports box constraints. Approximately 1.6x faster than \code{"cpp_L-BFGS-B"} on benchmark models. Recommended when speed is important.}
}

Defaults to \code{"nlminb"}.
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

