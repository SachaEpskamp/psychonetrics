\name{penalize}
\alias{penalize}
\alias{unpenalize}
\alias{penaltyVector}
\alias{refit}
\title{
Penalized Maximum Likelihood Functions
}
\description{
Functions for penalized maximum likelihood (PML) and penalized full-information maximum likelihood (PFIML) estimation. \code{penalize} adds an L1 (LASSO) or elastic net penalty to specified model parameters. \code{unpenalize} removes the penalty. \code{penaltyVector} extracts the per-free-parameter penalty strengths. \code{refit} performs post-selection inference by fixing penalized zeros and re-running with standard ML (or FIML for PFIML) estimation to obtain standard errors.
}
\usage{
penalize(x, matrix, row, col, lambda, group,
         what = c("free", "fixed"), verbose, log = TRUE)
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
The penalty strength. Can be a scalar (applied to all selected parameters) or a matrix with dimensions matching the model matrix. When a matrix is provided: \code{0} means do not penalize, \code{NA} means auto-select via EBIC, and a positive number means use that fixed penalty strength. If missing, uses the global \code{x@penalty$lambda}.
}
  \item{group}{
Integer indicating specific group(s). If missing, all groups are included.
}
  \item{what}{
Character, either \code{"free"} (default) or \code{"fixed"}. Controls which parameters are targeted: \code{"free"} penalizes parameters that are already free; \code{"fixed"} first frees parameters that are currently fixed to zero, then penalizes them (starting value set to 0). Parameters fixed to non-zero values (e.g., marker loadings fixed to 1) are never affected by \code{what = "fixed"}. When using a matrix \code{lambda}, mismatched parameters are silently skipped. When using a scalar \code{lambda} with specific \code{row}/\code{col}, a warning is issued for mismatched parameters.
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
Penalized ML estimation uses the \code{"PML"} estimator, which wraps the standard ML fit function and gradient with an additional penalty term. For models with missing data, the \code{"PFIML"} estimator wraps the FIML fit function and gradient with the same penalty. The ISTA (Iterative Shrinkage-Thresholding Algorithm) proximal gradient optimizer is used automatically for both PML and PFIML models.

For symmetric matrices (e.g., omega), only off-diagonal elements are penalized by default. For non-symmetric matrices (e.g., beta), all elements are penalized.

The \code{what} argument is particularly useful for regularized confirmatory factor analysis (CFA). For example, starting from a simple-structure CFA model, \code{what = "fixed"} can free and penalize all cross-loadings, allowing LASSO to select which cross-loadings to retain.

The \code{refit} function is used for post-selection inference: after a penalized model is estimated, \code{refit} fixes parameters that were shrunk to zero and re-estimates the model with standard ML (for PML) or FIML (for PFIML) to obtain valid standard errors and fit indices.

Model constructors (e.g., \code{ggm}, \code{lvm}, \code{var1}) accept \code{penalty_lambda}, \code{penalty_alpha}, and \code{penalize_matrices} arguments for convenient PML/PFIML setup. When \code{missing = "auto"} (the default), specifying \code{estimator = "PML"} with missing data will automatically switch to \code{"PFIML"}.
}
\value{
\code{penalize}, \code{unpenalize}, and \code{refit} return an object of class \code{psychonetrics} (\link{psychonetrics-class}).
\code{penaltyVector} returns a numeric vector of penalty strengths per free parameter.
}
\author{
Sacha Epskamp
}
\seealso{
\code{\link{find_penalized_lambda}} for automatic lambda selection via EBIC grid search,
\code{\link{runmodel}} for running models (calls \code{find_penalized_lambda} automatically when \code{penalty_lambda = NA}),
\code{\link{ggm}} and other model constructors for the \code{penalty_lambda}, \code{penalty_alpha}, and \code{penalize_matrices} arguments.
}
\examples{
\donttest{
# Load bfi data from psych package:
library("psychTools")
library("dplyr")
data(bfi)

ConsData <- bfi \%>\%
  select(A1:A5) \%>\%
  na.omit

vars <- names(ConsData)

# Penalized GGM estimation with manual lambda:
mod <- ggm(ConsData, vars = vars, estimator = "PML", penalty_lambda = 0.1)
mod <- mod \%>\% runmodel

# Post-selection refit for standard errors:
mod_refit <- mod \%>\% refit
mod_refit \%>\% parameters

# Manual unpenalize: free a specific edge from the penalty:
mod2 <- ggm(ConsData, vars = vars, estimator = "PML", penalty_lambda = 0.1)
mod2 <- unpenalize(mod2, "omega", row = 1, col = 2)
mod2 <- mod2 \%>\% runmodel
}
}
