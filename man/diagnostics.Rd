\name{diagnostics}
\alias{checkJacobian}
\alias{checkFisher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic functions
}
\description{
The 'checkJacobian' function can be used to check if the analytic gradient / 
Jacobian is aligned with the numerically approximated gradient / Jacobian, and 
the 'checkFisher' function can be used to check if the analytic Hessian is 
aligned with the numerically approximated Hessian.
}
\usage{
checkJacobian(x, f = "default", jac = "default", transpose = FALSE, 
          plot = TRUE, perturbStart = FALSE, method = "Richardson")

checkFisher(x, f = "default", fis = "default", transpose = FALSE, 
          plot = TRUE,  perturbStart = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A 'psychonetrics' object
}
  \item{f}{
A custom fit function or the psychonetrics default fit function (default).
}
  \item{jac}{
A custom Jacobian function or the psychonetrics default Jacobian function 
(default).
}
  \item{fis}{
A custom Fischer information function or the psychonetrics default Fischer 
information function (default).
}
  \item{transpose}{
Should the numeric Jacobian be transposed?
}
  \item{plot}{
Should a diagnostic plot be produced?
}
  \item{perturbStart}{
Should start values be perturbed (only used in development)
}
  \item{method}{
Numeric derivative method (default: Richardson)
}
}
\author{
Sacha Epskamp
}
