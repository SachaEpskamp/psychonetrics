\name{penalize}
\alias{penalize}
\alias{unpenalize}
\alias{penaltyVector}
\alias{refit}
\title{
Penalized Maximum Likelihood Functions
}
\description{
Functions for penalized maximum likelihood (PML) and penalized full-information maximum likelihood (FIPML) estimation. \code{penalize} adds an L1 (LASSO) or elastic net penalty to specified model parameters. \code{unpenalize} removes the penalty. \code{penaltyVector} extracts the per-free-parameter penalty strengths. \code{refit} performs post-selection inference by fixing penalized zeros and re-running with standard ML (or FIML for FIPML) estimation to obtain standard errors.
}
\usage{
penalize(x, matrix, row, col, lambda, group, verbose, log = TRUE)
unpenalize(x, matrix, row, col, group, verbose, log = TRUE)
penaltyVector(x)
refit(x, threshold = 1e-8, verbose, log = TRUE, ...)
}
\arguments{
  \item{x}{
A \code{psychonetrics} model.
}
  \item{matrix}{
Character vector of matrix name(s) to penalize/unpenalize. If missing, defaults to the matrices returned by \code{defaultPenalizeMatrices(x)} for \code{penalize}, or all currently penalized matrices for \code{unpenalize}.
}
  \item{row}{
Integer or character indicating specific row(s) of the matrix. If missing, the entire matrix is penalized/unpenalized.
}
  \item{col}{
Integer or character indicating specific column(s) of the matrix. If missing, the entire matrix is penalized/unpenalized.
}
  \item{lambda}{
Numeric, the penalty strength per parameter. If missing, uses the global \code{x@penalty$lambda}.
}
  \item{group}{
Integer indicating specific group(s). If missing, all groups are included.
}
  \item{verbose}{
Logical, should messages be printed?
}
  \item{log}{
Logical, should the log be updated?
}
  \item{threshold}{
For \code{refit}: threshold below which penalized parameters are considered zero and fixed.
}
  \item{\dots}{
Additional arguments passed to \code{runmodel} in \code{refit}.
}
}
\details{
Penalized ML estimation uses the \code{"PML"} estimator, which wraps the standard ML fit function and gradient with an additional penalty term. For models with missing data, the \code{"FIPML"} estimator wraps the FIML fit function and gradient with the same penalty. The ISTA (Iterative Shrinkage-Thresholding Algorithm) proximal gradient optimizer is used automatically for both PML and FIPML models.

For symmetric matrices (e.g., omega), only off-diagonal elements are penalized by default. For non-symmetric matrices (e.g., beta), all elements are penalized.

The \code{refit} function is used for post-selection inference: after a penalized model is estimated, \code{refit} fixes parameters that were shrunk to zero and re-estimates the model with standard ML (for PML) or FIML (for FIPML) to obtain valid standard errors and fit indices.

Model constructors (e.g., \code{ggm}, \code{lvm}, \code{var1}) accept \code{penalty_lambda}, \code{penalty_alpha}, and \code{penalize_matrices} arguments for convenient PML/FIPML setup. When \code{missing = "auto"} (the default), specifying \code{estimator = "PML"} with missing data will automatically switch to \code{"FIPML"}.
}
\value{
\code{penalize}, \code{unpenalize}, and \code{refit} return an object of class \code{psychonetrics} (\link{psychonetrics-class}).
\code{penaltyVector} returns a numeric vector of penalty strengths per free parameter.
}
\author{
Sacha Epskamp
}
\examples{
\dontrun{
# Penalized GGM estimation:
mod <- ggm(data, estimator = "PML", penalty_lambda = 0.1)
mod <- runmodel(mod)

# Post-selection refit for SEs:
mod_refit <- refit(mod)

# Manual penalize/unpenalize:
mod2 <- ggm(data, estimator = "PML", penalty_lambda = 0.1)
mod2 <- unpenalize(mod2, "omega", row = 1, col = 2)  # unpenalize specific edge
mod2 <- runmodel(mod2)
}
}
